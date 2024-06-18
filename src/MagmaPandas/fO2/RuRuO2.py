def fO2_RuRuO2(logshift, T_K, P_bar):
    """
    Armstrong et al. (2018) Supplementary materials equation S3.

    calibrated between 1 bar and 25 GPa and 773 - 2773 K

    O'Neill and Nell (1997) have an alternative formulation
    """

    P_GPa = P_bar / 1e4
    offset = 10**logshift

    log10fO2 = (
        (7.782 - 9.96e-3 * 1.932e-3 * P_GPa**2)
        + (-13763 + 593 * P_GPa - 3.955 * P_GPa**2) / T_K
        + (-1.05e6 - 4622 * P_bar) / T_K**2
    )

    fO2 = 10**log10fO2 * offset

    return fO2
