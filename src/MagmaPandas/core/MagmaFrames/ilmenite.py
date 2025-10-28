# import elementMass as e
import pandas as pd
from typing_extensions import Self

from MagmaPandas.core.enums import Datatype, Unit
from MagmaPandas.core.MagmaFrames.magmaFrame import MagmaFrame
from MagmaPandas.parse_io import check_components

_endmember_required_elements = ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "Cr2O3"]


class Ilmenite(MagmaFrame):
    """
    Subclass of :py:class:`~MagmaPandas.MagmaFrames.magmaFrame.MagmaFrame` extended with ilmenite specific methods.
    """

    def endmembers(self):
        """
        Calculations according to Andersen et al. (1993), as implemented in QUILF
        """
        composition = check_components(
            composition=self, components=_endmember_required_elements
        )

        cations = composition.cations(norm_to=2)

        hematite_moles = (cations[["Fe", "Mg", "Mn"]].sum(axis=1) - cations["Ti"]) / 2

        Fe2 = cations["Ti"] - cations["Mg"] - cations["Mn"]

        ilmenite_moles = Fe2 + cations["Al"] / 2

        Fe3 = 2 * hematite_moles

        total_components = (
            hematite_moles + ilmenite_moles + cations[["Mn", "Mg"]].sum(axis=1)
        )

        hematite = hematite_moles / total_components
        ilmenite = ilmenite_moles / total_components
        geikielite = cations["Mg"] / total_components
        pyrophanite = cations["Mn"] / total_components

        return pd.DataFrame(
            {
                "hematite": hematite,
                "ilmenite": ilmenite,
                "geikielite": geikielite,
                "pyrophanite": pyrophanite,
            }
        )

    def Fe_speciation(self, normalise=False) -> Self:
        """
        Calculations according to Andersen et al. (1993), as implemented in QUILF
        """

        units = self._units
        datatype = self._datatype

        elements = {Datatype.OXIDE: ("FeO", "Fe2O3"), Datatype.CATION: ("Fe", "Fe3")}[
            datatype
        ]

        if all(c in self.elements for c in elements):
            return self

        composition = check_components(
            composition=self, components=_endmember_required_elements
        )

        cations = composition.cations(
            normalise=False
        )  # convert to cation mol fractions
        totals = cations["total"]
        cations = cations.div(totals, axis=0).mul(2, axis=0)

        hematite_moles = (cations[["Fe", "Mg", "Mn"]].sum(axis=1) - cations["Ti"]) / 2

        Fe2 = cations["Ti"] - cations["Mg"] - cations["Mn"]
        Fe3 = 2 * hematite_moles

        cations["Fe"] = Fe2
        cations["Fe3"] = Fe3
        cations = cations.recalculate()
        cations = cations.div(2, axis=0).mul(totals, axis=0)  # revert normalisation

        if units == Unit.MOL_FRACTIONS:
            if datatype == Datatype.CATION:
                return cations
            return cations.oxides(normalise=normalise)

        cations_wt_pc = cations.wt_pc(normalise=normalise)
        if datatype == Datatype.CATION:
            return cations_wt_pc

        return cations_wt_pc.oxides(normalise=normalise)
