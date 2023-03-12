from MagmaPandas.Kd.Ol_melt import calculate_FeMg_Kd
from MagmaPandas.Fe_redox import Fe_redox
from MagmaPandas.fO2 import fO2_QFM


from MagmaPandas.configuration import configuration


def _root_temperature(olivine_amount, melt_x_moles, olivine_x_moles, T_K, P_bar):

    melt_x_new = melt_x_moles + olivine_x_moles.mul(olivine_amount)
    melt_x_new = melt_x_new.normalise()
    temperature_new = melt_x_new.convert_moles_wtPercent.melt_temperature(P_bar=P_bar)

    return T_K - temperature_new


def _root_Kd(exchange_amount, melt_x_moles, exchange_vector, forsterite, P_bar, kwargs):

    melt_x_new = melt_x_moles + exchange_vector.mul(exchange_amount)
    melt_x_new = melt_x_new.normalise()
    Kd_equilibrium, Kd_real = calculate_Kds(melt_x_new, P_bar, forsterite, **kwargs)

    return Kd_equilibrium - Kd_real


def calculate_Kds(melt_x_moles, P_bar, forsterite, **kwargs):

    Fe3Fe2_model = getattr(Fe_redox, configuration.Fe3Fe2_model)
    Kd_model = calculate_FeMg_Kd
    dQFM = kwargs.get("dQFM", configuration.dQFM)

    T_K = melt_x_moles.convert_moles_wtPercent.melt_temperature(P_bar)
    fO2 = fO2_QFM(dQFM, T_K, P_bar)
    Fe3Fe2 = Fe3Fe2_model(melt_x_moles, T_K, fO2)

    Kd_observed = calculate_observed_Kd(melt_x_moles, Fe3Fe2, forsterite)

    Kd_eq = Kd_model(melt_x_moles, forsterite, T_K, Fe3Fe2, P_bar=P_bar)

    return Kd_eq, Kd_observed


def calculate_observed_Kd(melt_x_moles, Fe3Fe2, forsterite):

    melt_x_moles = melt_x_moles.normalise()
    Fe2_FeTotal = 1 / (1 + Fe3Fe2)
    melt_MgFe = melt_x_moles["MgO"] / (melt_x_moles["FeO"] * Fe2_FeTotal)
    olivine_MgFe = forsterite / (1 - forsterite)

    return melt_MgFe / olivine_MgFe
