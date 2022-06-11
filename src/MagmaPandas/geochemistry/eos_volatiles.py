import numpy as np
from ..parse.validate import _check_argument

"""
Re-written visual basic code from Holloway and Blank, 1994 (Reviews
in Mineralogy and Geochemistry, vol. 30). 
Original python code from the supplementary material of:

C. M. Allison, K. Roggensack & A. B. Clarke (2022) 
MafiCH: a general model for H2O–CO2 solubility in mafic magmas. 
Contributions to Mineralogy and Petrology 177: 40
"""


class hollowayBlank:
    @staticmethod
    @_check_argument("species", [None, "H2O", "CO2"])
    def fugacity(T_K, P_bar, species):
        """
        Calculates CO2 fugacity given temperature in Kelvin and pressure
        in bar
        """
        P0 = 4000

        # MRK below 4kb and Saxena above 4kb
        if P_bar > 4000 and species == "CO2":
            # Calculate FP at T and 4Kb
            iPUREG = hollowayBlank._RKCALC(T_K, P0, species)
            XLNF = hollowayBlank._Saxena(T_K, P_bar)
            PUREG = iPUREG + XLNF
        else:
            PUREG = hollowayBlank._RKCALC(T_K, P_bar, species)

        # Convert from LN fugacity to fugacity
        stdf = np.exp(PUREG)
        return stdf

    @staticmethod
    def _RKCALC(T_K, P_bar, species):
        """
        Calculation of pure gas MRK (modified Redlick-Kwong) properties following Holloway 1981, 1987
        P is in bar, T in Kelvin. Returns ln fugacity
        """

        R = 82.05736
        PBLN = np.log(P_bar)
        T_C = T_K - 273.15
        RXT = R * T_K
        RT = R * T_K ** 1.5 * 0.000001  # 1e-6

        if species == "CO2":
            # Calculate T dependent MRK A parameter CO2
            ACO2M = 73.03 - 0.0714 * T_C + 2.157e-05 * T_C * T_C
            # Define MRK B parameter for CO2
            BSUM = 29.7
            ASUM = ACO2M / (BSUM * RT)

        elif species == "H2O":
            AH2OM = 115.98 - 0.0016295 * T_K - 1.4984e-05 * T_K ** 2
            BSUM = 14.5
            ASUM = AH2OM / (BSUM * RT)
        else:
            raise ValueError(
                f"species: {species} not recognised, please choose 'CO2' or 'H2O'"
            )

        BSUM = P_bar * BSUM / RXT
        XLNFP = hollowayBlank._REDKW(BSUM, ASUM)
        # Convert to LN fugacity
        PUREG = XLNFP + PBLN
        return PUREG

    @staticmethod
    def _REDKW(BP, A2B):
        """
        The RK routine
        A routine to calc compressibility factor and fugacity coefficient with the
        Redlick-Kwong equation following Edmister (1968). This solution for
        supercritical fluid returns a value of 1 for FP if arguments are out of range.
        """

        TH = 1 / 3

        if A2B < 1e-10:
            A2B = 0.001

        RR = -A2B * BP * BP
        QQ = BP * (A2B - BP - 1)
        XN = QQ * TH + RR - 0.074074
        XM = QQ - TH
        XNN = XN * XN / 4
        XMM = XM * XM * XM / 27
        ARG = XNN + XMM

        if ARG > 0:
            X = np.sqrt(ARG)
            F = 1
            XN2 = -XN / 2
            iXMM = XN2 + X
            if iXMM < 0:
                F = -1
            XMM = F * ((F * iXMM) ** TH)
            F = 1
            iXNN = XN2 - X
            if iXNN < 0:
                F = -1
            XNN = F * ((F * iXNN) ** TH)
            Z = XMM + XNN + TH
            ZBP = Z - BP
            if ZBP < 0.000001:
                ZBP = 0.000001
            BPZ = 1 + BP / Z
            FP = Z - 1 - np.log(ZBP) - A2B * np.log(BPZ)
            if FP < -37 or FP > 37:
                FP = 0.000001

        elif ARG < 0:
            COSPHI = np.sqrt(-XNN / XMM)
            if XN > 0:
                COSPHI = -COSPHI
            TANPHI = np.sqrt(1 - COSPHI * COSPHI) / COSPHI
            PHI = np.arctan(TANPHI) * TH
            FAC = 2 * np.sqrt(-XM * TH)
            # sort for largest root
            R1 = np.cos(PHI)
            R2 = np.cos(PHI + 2.0944)
            R3 = np.cos(PHI + 4.18879)
            RH = R2
            if R1 > R2:
                RH = R1
            if R3 > RH:
                RH = R3
            Z = RH * FAC + TH
            ZBP = Z - BP
            if ZBP < 0.000001:
                ZBP = 0.000001
            BPZ = 1 + BP / Z
            FP = Z - 1 - np.log(ZBP) - A2B * np.log(BPZ)
            if FP < -37 or FP > 37:
                FP = 0.000001

        else:
            FP = 1
            Z = 1
        XLNFP = FP
        return XLNFP

    @staticmethod
    def _Saxena(T_K, P_bar):
        """
        High pressure corresponding states routines from Saxena and Fei (1987) GCA
        Vol. 51, 783-791
        Returns natural np.log of the ratio F(P)/F(4000 bar) as XLNF array
        """

        # Define integration limit
        PO = 4000

        # Critical temperatures and pressures for CO2
        TR = T_K / 304.2
        PR = P_bar / 73.9
        PC = 73.9

        # Virial coefficients
        A = 2.0614 - 2.2351 / TR ** 2 - 0.39411 * np.log(TR)
        B = 0.055125 / TR + 0.039344 / TR ** 2
        C = -1.8935e-06 / TR - 1.1092e-05 / TR ** 2 - 2.1892e-05 / TR ** 3
        D = 5.0527e-11 / TR - 6.3033e-21 / TR ** 3

        # Calculate molar volume
        Z = A + B * PR + C * PR ** 2 + D * PR ** 3
        V = Z * 83.0117 * T_K / P_bar

        # integrate from PO (4000 bars) to P to calculate ln fugacity
        LNF = (
            A * np.log(P_bar / PO)
            + (B / PC) * (P_bar - PO)
            + (C / (2 * PC ** 2)) * (P_bar ** 2 - PO ** 2)
        )
        LNF = LNF + (D / (3 * PC ** 3)) * (P_bar ** 3 - PO ** 3)
        XLNF = LNF
        return XLNF
