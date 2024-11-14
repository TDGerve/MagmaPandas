import elementMass as e
import pandas as pd
from typing_extensions import Self

from MagmaPandas.enums import Datatype, Unit
from MagmaPandas.MagmaFrames.magmaFrame import MagmaFrame
from MagmaPandas.parse_io import check_components

_endmember_required_elements = ["MgO", "MnO", "FeO", "Fe2O3", "TiO2", "Al2O3"]


class Magnetite(MagmaFrame):
    """
    Subclass of :py:class:`~MagmaPandas.MagmaFrames.magmaFrame.MagmaFrame` extended with magnetite specific methods.
    """

    def endmembers(self):
        """
        Calculations according to Lindsley, as implemented in QUILF
        """

        composition = check_components(
            composition=self, components=_endmember_required_elements
        )

        cations = composition.cations(norm_to=3)

        magnetite = (
            cations[["Mg", "Mn", "Fe"]].sum(axis=1)
            - 2 * cations["Ti"]
            - cations["Al"] / 2
        ) / 3
        mg_ulvospinel = cations["Mg"] / 2
        mn_ulvospinel = cations["Mn"] / 2
        ulvospinel = cations["Ti"] - mn_ulvospinel - mg_ulvospinel

        total = magnetite + ulvospinel + mg_ulvospinel
        ulvospinel = (ulvospinel + mg_ulvospinel) / total
        magnetite = 1 - ulvospinel

        return pd.DataFrame(
            {"magnetite": magnetite, "ulvospinel": ulvospinel}, index=self.index
        )

    def Fe_speciation(self) -> Self:
        """
        Calculations according to Lindsley, as implemented in QUILF
        """

        # TODO add code to catch when total iron is given as Fe2O3

        old_datatype = self._datatype
        elements = {Datatype.OXIDE: ("FeO", "Fe2O3"), Datatype.CATION: ("Fe", "Fe3")}[
            old_datatype
        ]

        if all(c in self.elements for c in elements):
            return self

        composition = check_components(
            composition=self, components=_endmember_required_elements
        )

        cations = composition.cations(normalise=False)

        totals = cations["total"]

        cations_norm = cations.div(totals, axis=0) * 3  # normalised to 3 cations total

        magnetite = (
            cations_norm[["Mg", "Mn", "Fe"]].sum(axis=1)
            - 2 * cations_norm["Ti"]
            - cations_norm["Al"] / 2
        ) / 3

        Fe3 = 2 * magnetite
        Fe2 = (
            magnetite
            + 2 * cations_norm["Ti"]
            + cations_norm["Al"] / 2
            - cations_norm["Mg"]
            - cations_norm["Mn"]
        )

        # cations_norm = cations_norm.drop(columns=["Fe"])
        cations_norm["Fe3"] = Fe3
        cations_norm["Fe"] = Fe2

        cations_Fe_calc = (cations_norm / 3).mul(totals, axis=0).recalculate()

        if old_datatype == Datatype.CATION:
            return cations_Fe_calc

        oxide_wt_pc = cations_Fe_calc.oxides(normalise=False).wt_pc(normalise=False)
        oxide_wt_pc["total"] = oxide_wt_pc[oxide_wt_pc.elements].sum(axis=1)

        # oxide_names = e.get_oxide_names(cations_Fe_calc.elements)
        # oxide_weights = e.compound_weights(oxide_names)
        # cation_amounts = e.cation_numbers(oxide_names)

        # oxide_wt_pc = cations_Fe_calc.rename(
        #     columns={
        #         cation: oxide
        #         for cation, oxide in zip(cations_Fe_calc.elements, oxide_names)
        #     }
        # ).recalculate()

        # oxide_wt_pc[oxide_wt_pc.elements] = (
        #     oxide_wt_pc[oxide_wt_pc.elements]
        #     .mul(oxide_weights, axis=1)
        #     .div(cation_amounts, axis=1)
        # )

        # oxide_wt_pc._units = Unit.WT_PERCENT
        # oxide_wt_pc._datatype = Datatype.OXIDE

        return oxide_wt_pc
