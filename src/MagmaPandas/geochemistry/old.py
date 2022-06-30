def Fe_restore(
    inclusion: Union[Melt_inclusion, Melt],
    FeO_target: Union[int, float],
    P_bar: float,
    **kwargs,
):
    """
    Reverse melt inclusion Fe loss (or gain) through Fe-Mg exchange until a set initial melt FeO is reached.
    Isothermal and isobaric.
    """
    # Grab model parameters
    converge = kwargs.get(
        "converge", getattr(PEC_configuration, "FeO_converge")
    )  # In wt. %
    temperature_converge = kwargs.get(
        "temperature_converge", getattr(PEC_configuration, "temperature_converge")
    )  # In degrees
    stepsize = kwargs.get(
        "stepsize", getattr(PEC_configuration, "stepsize")
    )  # In molar fraction
    # Select Fe loss or gain
    if inclusion["FeO"] > FeO_target:
        stepsize = -stepsize
    forsterite = kwargs.get("forsterite", 0.8)
    QFM_logshift = kwargs.get("QFM_logshift", configuration().QFMlogshift)
    # Parameters for the while loop
    olivine_corrected = 0
    olivine_stepsize_reduction = 4
    decrease_factor = getattr(PEC_configuration, "decrease_factor")
    # Collect configured models
    Fe3Fe2_model = getattr(Fe_redox, configuration().Fe3Fe2_model)
    Kd_model = getattr(Kd_FeMg, configuration().Kd_model)

    # Normalise inclusion composition
    inclusion = inclusion[inclusion.elements]
    inclusion = inclusion.normalise()
    # Calculate temperature and fO2
    temperature = inclusion.melt_temperature(P_bar=P_bar)
    fO2 = fO2_QFM(QFM_logshift, temperature, P_bar)

    # Set up initial data
    moles = Melt_inclusion(
        columns=inclusion.elements, units="mol fraction", datatype="oxide"
    )
    moles.loc[0] = inclusion.moles[inclusion.elements].values
    # Fe-Mg exchange vector
    FeMg_exchange = pd.Series(0, index=moles.columns)
    FeMg_exchange.loc[["FeO", "MgO"]] = stepsize, -stepsize
    # Inclusion starting FeO
    FeO = inclusion["FeO"]

    while not np.isclose(FeO, FeO_target, atol=converge):
        # Exchange Fe and Mg
        idx = moles.index[-1] + stepsize
        moles.loc[idx] = (moles.iloc[-1] + FeMg_exchange).values
        # New model temperature after Fe-Mg exchange
        temperature_new = moles.iloc[-1].convert_moles_wtPercent.melt_temperature(
            P_bar=P_bar
        )

        # Calculate equilibrium olivine Fo#
        # melt Fe3+/Fe2+
        Fe3Fe2 = Fe3Fe2_model(moles.iloc[-1], temperature, fO2)
        # FeMg Kd
        Kd = Kd_model(moles.iloc[-1], forsterite, temperature, P_bar, Fe3Fe2)
        Fe2_FeTotal = 1 / (1 + Fe3Fe2)
        Fe2Mg = moles.loc[idx, "FeO"] * Fe2_FeTotal / moles.loc[idx, "MgO"]
        # Equilibrium forsterite content
        forsterite_EQ = 1 / (1 + Kd * Fe2Mg)
        # Equilibrium olivine composition in oxide mol fraction
        olivine = MagmaSeries(
            {"MgO": forsterite_EQ * 2, "FeO": (1 - forsterite_EQ) * 2, "SiO2": 1}
        )
        olivine = olivine.normalise().reindex(moles.columns, fill_value=0.0)

        # Crystallise or melt olivine to remain isothermal
        add_olivine = moles.iloc[-1]
        # Set stepsize for olivine addition/removal
        olivine_stepsize = stepsize / olivine_stepsize_reduction
        while not np.isclose(temperature_new, temperature, atol=temperature_converge):
            # crystallise olivine until calculated temperature is back at initial
            add_olivine = add_olivine + olivine * olivine_stepsize
            temperature_new = add_olivine.convert_moles_wtPercent.melt_temperature(
                P_bar=P_bar
            )
            # Record added/removed olivine amount. Negative values for crystallisation
            olivine_corrected += olivine_stepsize
            T_overstepped = np.sign(temperature - temperature_new) != np.sign(stepsize)
            # Reverse one iteration and reduce stepsize if temperature was
            # overstepped by more than the convergence value
            if T_overstepped and not np.isclose(
                temperature_new, temperature, atol=temperature_converge
            ):
                add_olivine = add_olivine - olivine * olivine_stepsize
                olivine_corrected -= olivine_stepsize
                olivine_stepsize = olivine_stepsize / decrease_factor
                continue

        # Copy olivine corrected composition
        moles.iloc[-1] = add_olivine.values

        # New inclusion FeO
        FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]

        overstepped = np.sign(FeO_target - FeO) != np.sign(stepsize)
        # Reverse one iteration and reduce stepsize if FeO content
        # gets oversteppend by more than the convergence value
        if overstepped and not np.isclose(FeO, FeO_target, atol=converge):
            moles.drop([idx], inplace=True)
            FeO = moles.iloc[-1].convert_moles_wtPercent["FeO"]
            stepsize = stepsize / decrease_factor
            FeMg_exchange = FeMg_exchange.div(decrease_factor)

    # Recalculate compositions to oxide wt. %
    wtPercent = moles.convert_moles_wtPercent

    return wtPercent, temperature, temperature_new, olivine_corrected