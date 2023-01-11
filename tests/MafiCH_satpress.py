#This script calculates saturation pressures for a series of samples from a csv
#input file using the MafiCH model of Allison et al. (Contributions to 
#Mineralogy and Petrology)

#The only section of the code that should be edited by users is the USER INPUTS
#section that directly follows this header

#For further information about using this script, including the template for the
#required input file, refer to the supplementary material for this publication

#All calculations using this script or other aspects of the MafiCH model should
#cite: Allison et al., MafiCH: a general model for H2O-CO2 solubility in mafic
#magmas. Contributions to Mineralogy and Petrology

#Contact author Dr. Chelsea Allison at cmalliso@asu.edu 

################################################################################
#USER INPUTS
################################################################################

#This is the only section of the code that should be changed

#SET THE PATH TO THE INPUT FILE:
#USE THE PROVIDED TEMPLATE TO CONSTRUCT THE CSV FILE
satpressinput = '/Users/Chelsea/Desktop/MafiCH_satpress_input.csv'

#SET THE PATH TO THE OUTPUT FILE:
#output file does not need to exist; it will be generated during the calculations
satpressoutput = '/Users/Chelsea/Desktop/MafiCH_satpress_output.csv'

#End of user input section. No further changes needed to run the code

################################################################################
#define function to calculate TD parameters with model
################################################################################

import math
import csv

def TDparam(wtSiO2, wtTiO2, wtAl2O3, wtFeO, wtMgO, wtCaO, wtNa2O, wtK2O):

    #formula weights
    FWSiO2 = 60.08
    FWTiO2 = 79.866
    FWAl2O3 = 101.96
    FWFeO = 71.844
    FWMgO = 40.3044
    FWCaO = 56.0774
    FWNa2O = 61.9789
    FWK2O = 94.2
    
    #divide wt.% by formula weight
    SiO2 = wtSiO2/FWSiO2
    TiO2 = wtTiO2/FWTiO2
    Al2O3 = wtAl2O3/FWAl2O3
    FeO = wtFeO/FWFeO
    MgO = wtMgO/FWMgO
    CaO = wtCaO/FWCaO
    Na2O = wtNa2O/FWNa2O
    K2O = wtK2O/FWK2O
    
    #Oxygen Factor
    OSiO2 = 2
    OTiO2 = 2
    OAl2O3 = 3
    OFeO = 1
    OMgO = 1
    OCaO = 1
    ONa2O = 1
    OK2O = 1
    
    #Cations per Oxygen
    CSiO2 = 1/2
    CTiO2 = 1/2
    CAl2O3 = 2/3
    CFeO = 1
    CMgO = 1
    CCaO = 1
    CNa2O = 2/1
    CK2O = 2/1
    
    #Cation
    CatSi = SiO2*OSiO2*CSiO2
    CatTi = TiO2*OTiO2*CTiO2
    CatAl = Al2O3*OAl2O3*CAl2O3
    CatFe2 = FeO*OFeO*CFeO
    CatMg = MgO*OMgO*CMgO
    CatCa = CaO*OCaO*CCaO
    CatNa = Na2O*ONa2O*CNa2O
    CatK = K2O*OK2O*CK2O
    
    #add cations
    CatSum = CatSi+CatTi+CatAl+CatFe2+CatMg+CatCa+CatNa+CatK
    
    #normalize to cation fraction
    Si = CatSi/CatSum
    Ti = CatTi/CatSum
    Al = CatAl/CatSum
    Fe2 = CatFe2/CatSum
    Mg = CatMg/CatSum
    Ca = CatCa/CatSum
    Na = CatNa/CatSum
    K = CatK/CatSum
    
    #round to thousandths
    Si = round(Si,3)
    Ti = round(Ti,3)
    Al = round(Al,3)
    Fe2 = round(Fe2,3)
    Mg = round(Mg,3)
    Ca = round(Ca,3)
    Na = round(Na,3)
    K = round(K,3)
    
    #calculate alkali ratio
    NaK = Na/(Na+K);
    NaK = round(NaK,3)
    
    #calculate deltaV and lnK0 using the model equations
    deltaV = -3350.65 + 2625.385*Ti + 3105.426*Al + 47.0037*NaK + 3375.552*(Si+Na) + 3795.115*K + 3628.018*Fe2 + 3323.32*(Mg+Ca)
    deltaV = round(deltaV,2)
    
    lnK0 = -128.365 + 122.644*(Fe2+Na+Ca) + 92.263*(Ti+Al) + 114.098*Si + 111.549*Mg + 138.855*K + 2.239*NaK;
    lnK0 = round(lnK0,2)
    
    return deltaV, lnK0

################################################################################
#define CO2 fugacity calculation function
################################################################################

#Re-written from visual basic code from Holloway and Blank, 1994 (Reviews
#in Mineralogy and Geochemistry, vol. 30)

#set up 3 functions needed for this calculation
#first function: RKCALC

def RKCALC(T, P):
    #Calculation of pure gas MRK properties following Holloway 1981, 1987
    #P is in atm, T in Kelvin. Returns ln fugacity

    R = 82.05736
    RR = 6732.2
    pb = 1.013*P
    PBLN = math.log(pb)
    TCEL = T-273.15
    RXT = R*T
    RT = R*T**1.5*.000001

    #Calculate T dependent MRK A parameter CO2
    ACO2M = 73.03-.0714*TCEL+2.157e-05*TCEL*TCEL

    #Define MRK B parameter for CO2
    BSUM = 29.7
    ASUM = ACO2M/(BSUM*RT)
    BSUM = P*BSUM/RXT
    XLNFP = REDKW(BSUM, ASUM)
    #Convert to LN fugacity
    PUREG = XLNFP+PBLN
    return PUREG

################################################################################

#second function: REDKW

def REDKW(BP, A2B):
    #The RK routine
    #A routine to calc compressibility factor and fugacity coefficient with the
    #Redlick-Kwong equation following Edmister (1968). This solution for 
    #supercritical fluid
    #Returns a value of 1 for FP if arguments out of range

    TH = 1/3

    if A2B < 1E-10:
        A2B = .001

    RR = -A2B*BP*BP
    QQ = BP*(A2B-BP-1)
    XN = QQ*TH+RR-.074074
    XM = QQ-TH
    XNN = XN*XN/4
    XMM = XM*XM*XM/27
    ARG = XNN+XMM

    if ARG > 0:
        X = math.sqrt(ARG)
        F = 1
        XN2 = -XN/2
        iXMM = XN2+X
        if iXMM < 0:
            F = -1
        XMM = F*((F*iXMM)**TH)
        F = 1
        iXNN = XN2-X
        if iXNN < 0:
            F= -1
        XNN = F*((F*iXNN)**TH)
        Z = XMM+XNN+TH
        ZBP = Z-BP
        if ZBP < .000001:
            ZBP = .000001
        BPZ = 1+BP/Z
        FP = Z-1-math.log(ZBP)-A2B*math.log(BPZ)
        if FP < -37 or FP > 37:
            FP = .000001
        XLNFP = FP
    elif ARG < 0:
        COSPHI = math.sqrt(-XNN/XMM)
        if XN > 0:
            COSPHI = -COSPHI
        TANPHI = math.sqrt(1-COSPHI*COSPHI)/COSPHI
        PHI = math.atan(TANPHI)*TH
        FAC = 2*math.sqrt(-XM*TH)
        #sort for largest root
        R1 = math.cos(PHI)
        R2 = math.cos(PHI+2.0944)
        R3 = math.cos(PHI+4.18879)
        RH = R2
        if R1 > R2:
            RH = R1
        if R3 > RH:
            RH = R3
        Z = RH*FAC+TH
        ZBP = Z-BP
        if ZBP < .000001:
            ZBP = .000001
        BPZ = 1+BP/Z
        FP = Z-1-math.log(ZBP)-A2B*math.log(BPZ)
        if FP < -37 or FP > 37:
            FP = .000001
        XLNFP = FP
    else:
        FP = 1
        Z = 1
        XLNFP = FP
    return XLNFP

################################################################################

#third function: Saxena

def Saxena(TK, pb):
    #High pressure corresponding states routines from Saxena and Fei (1987) GCA
    #Vol. 51, 783-791
    #Returns natural math.log of the ratio F(P)/F(4000 bar) as XLNF array

    #Define integration limit
    PO = 4000

    #Critical temperatures and pressures for CO2
    TR = TK/304.2
    PR = pb/73.9
    PC = 73.9

    #Virial coefficients
    A = 2.0614-2.2351/TR**2-.39411*math.log(TR)
    B = .055125/TR+.039344/TR**2
    C = -1.8935e-06/TR-1.1092e-05/TR**2-2.1892e-05/TR**3
    D = 5.0527e-11/TR-6.3033e-21/TR**3

    #Calculate molar volume
    Z = A+B*PR+C*PR**2+D*PR**3
    V = Z*83.0117*TK/pb

    #integrate from PO (4000 bars) to P to calculate ln fugacity
    LNF = A*math.log(pb/PO)+(B/PC)*(pb-PO)+(C/(2*PC**2))*(pb**2-PO**2)
    LNF = LNF+(D/(3*PC**3))*(pb**3-PO**3)
    XLNF = LNF
    return XLNF

################################################################################

#now use all three functions to determine CO2 fugacity

def fco2calc(TC, pb):
    #Calculates CO2 fugacity given temperature in degrees Celsius and pressure
    #in bar

    #Convert to atmospheres and Kelvin
    PO = 4000/1.013
    PA = pb/1.013
    TK = TC+273.15

    #MRK below 4kb and Saxena above 4kb
    if pb > 4000:
            #Calculate FP at T and 4Kb
            iPUREG = RKCALC(TK, PO)
            XLNF = Saxena(TK, pb)
            PUREG = iPUREG+XLNF
    else:
            PUREG = RKCALC(TK, PA)

    #Convert from LN fugacity to fugacity
    stdfco2 = math.exp(PUREG)
    return stdfco2

################################################################################
#define H2O fugacity calculation function
################################################################################

def fh2ocalc(TC, pb):
    #calculate H2O fugacity given temperature in degrees Celsius and pressure
    #in bar

    #rewritten from visual basic code from Holloway and Blank, 1994 (Reviews in
    #Mineramath.logy and Geochemistry vol. 30)

    #Convert to Kelvin
    TK = TC+273.15

    R = 83.144621
    b = 14.5

    RXT = R*TK
    RT = R*TK**1.5*.000001

    ah2o = 115.98-.0016295*TK-1.4984e-05*TK**2

    A2B = ah2o/(b*RT)
    BP = pb*b/RXT

    TH = 1/3

    if A2B < 1e-10:
        A2B = .001

    RR = -A2B*BP*BP
    QQ = BP*(A2B-BP-1)
    XN = QQ*TH+RR-.074074
    xm = QQ-TH
    XNN = XN*XN/4
    XMM = xm*xm*xm/27
    ARG = XNN+XMM

    if ARG > 0:
        X = math.sqrt(ARG)
        F = 1
        XN2 = -XN/2
        iXMM = XN2+X
        if iXMM < 0:
            F = -1
        XMM = F*((F*iXMM)**TH)
        F = 1
        iXNN = XN2-X
        if iXNN < 0:
            F= -1
        XNN = F*((F*iXNN)**TH)
        Z = XMM+XNN+TH
        ZBP = Z-BP
        if ZBP < .000001:
            ZBP = .000001
        BPZ = 1+BP/Z
        FP = Z-1-math.log(ZBP)-A2B*math.log(BPZ)
        if FP < -37 or FP > 37:
            FP = .000001
    elif ARG < 0:
        COSPHI = math.sqrt(-XNN/XMM)
        if XN > 0:
            COSPHI = -COSPHI
        TANPHI = math.sqrt(1-COSPHI*COSPHI)/COSPHI
        PHI = math.atan(TANPHI)*TH
        FAC = 2*math.sqrt(-xm*TH)
        #sort for largest root
        R1 = math.cos(PHI)
        R2 = math.cos(PHI+2.0944)
        R3 = math.cos(PHI+4.18879)
        RH = R2
        if R1 > R2:
            RH = R1
        if R3 > RH:
            RH = R3
        Z = RH*FAC+TH
        ZBP = Z-BP
        if ZBP < .000001:
            ZBP = .000001
        BPZ = 1+BP/Z
        FP = Z-1-math.log(ZBP)-A2B*math.log(BPZ)
        if FP < -37 or FP > 37:
            FP = .000001
    else:
        FP = 1
        Z = 1

    PUREG = FP + np.log(P_bar)
    stdfh2o = math.exp(PUREG)
    return stdfh2o

################################################################################
#define saturation pressure function
################################################################################

def SatP(H2O,CO2,lnK0,deltaV):
    #calculates saturation pressure and fluid composition
    #requires dissolved H2O & CO2 and thermodynamic parameters lnK0 & deltaV
    #needs access to fco2calc and fh2ocalc

    #temperature 1200 C
    TC = 1200

    #Determine partial pressure of CO2 for this CO2 concentration using solubility model
    #work backwards from concentration: determine XmCO3, then Kf
    XCO3 = CO2*(1/44.01)/(10**4*(100/36.594)-(CO2/36.594))
    Kf = XCO3/(1+XCO3)
    #converting Kf to K requires CO2 fugacity, which requires pressure
    #start at low pressure
    CO2pp = 1
    #calculate K from the CO2 fugacity (K = Kf/fugacity)
    KbyFug = Kf/fco2calc(TC,CO2pp)
    #calculate K from the partial pressure in the solubility equation
    KbyPress = math.exp(lnK0)*math.exp(-deltaV*(CO2pp-1000)/(83.144621*(TC+273.15)))
    #determine difference between K calculated by each method
    Kdiff = KbyFug - KbyPress
    #loop until the two methods of K value determinations converge
    while Kdiff > 0:
        Kdiffold = Kdiff
        CO2ppold = CO2pp
        CO2pp = CO2pp + 100    #large step size until difference is negative
        KbyFug = Kf/fco2calc(TC,CO2pp)
        KbyPress = math.exp(lnK0)*math.exp(-deltaV*(CO2pp-1000)/(83.144621*(TC+273.15)))
        Kdiff = KbyFug - KbyPress

    while Kdiff < 0:
        Kdiffold = Kdiff
        CO2ppold = CO2pp
        CO2pp = CO2pp - 10     #smaller steps until difference becomes positive again
        KbyFug = Kf/fco2calc(TC,CO2pp)
        KbyPress = math.exp(lnK0)*math.exp(-deltaV*(CO2pp-1000)/(83.144621*(TC+273.15)))
        Kdiff = KbyFug - KbyPress

    while Kdiff > 0:
        Kdiffold = Kdiff
        CO2ppold = CO2pp
        CO2pp = CO2pp + 1      #1 bar increments as difference converges to 0
        KbyFug = Kf/fco2calc(TC,CO2pp)
        KbyPress = math.exp(lnK0)*math.exp(-deltaV*(CO2pp-1000)/(83.144621*(TC+273.15)))
        Kdiff = KbyFug - KbyPress

    #if previous pressure step yielded a difference closer to 0, use that pressure 
    if Kdiffold < abs(Kdiff):
        CO2pp = CO2ppold

    #calculate CO2 fugacity for this partial pressure
    fCO2 = fco2calc(TC,CO2pp)

    #determine H2O fugacity for this H2O concentration
    #Lesne Etna power equation
    fH2O = 104.98*H2O**1.83

    #determine pressure and fluid composition at which these fugacities can coexist
    #start at ppCO2 value (has to be higher than that for mixed fluid compositions)
    P = CO2pp
    #determine fluid composition based on H2O fugacity compared to pure H2O fugacity at this pressure
    fH2Opure = fh2ocalc(TC,P)
    XfH2O = fH2O/fH2Opure
    #determine CO2 fugacity at this pressure and the fluid composition according to H2O
    fCO2pure = fco2calc(TC,P)
    fCO2trial = fCO2pure*(1-XfH2O)
    #compare this CO2 fugacity to the known CO2 fugacity
    fdiff = fCO2 - fCO2trial
    #loop until these two determinations of CO2 fugacity converge
    while fdiff > 0:
        Pold = P
        fdiffold = fdiff
        XfH2Oold = XfH2O
        P = P + 100            #large step size until difference becomes negative
        #calculate new fugacity difference at this new P and fluid composition
        fH2Opure = fh2ocalc(TC,P)
        XfH2O = fH2O/fH2Opure
        fCO2pure = fco2calc(TC,P)
        fCO2trial = fCO2pure*(1-XfH2O)
        fdiff = fCO2 - fCO2trial
    
    while fdiff < 0:
        Pold = P
        fdiffold = fdiff
        XfH2Oold = XfH2O
        P = P - 10             #smaller steps until difference becomes positive again
        #calculate new fugacity difference at this new P and fluid composition
        fH2Opure = fh2ocalc(TC,P)
        XfH2O = fH2O/fH2Opure
        fCO2pure = fco2calc(TC,P)
        fCO2trial = fCO2pure*(1-XfH2O)
        fdiff = fCO2 - fCO2trial
    
    while fdiff > 0:
        Pold = P
        fdiffold = fdiff
        XfH2Oold = XfH2O
        P = P + 1              #1 bar steps as the difference converges to 0
        #calculate new fugacity difference at this new P and fluid composition
        fH2Opure = fh2ocalc(TC,P)
        XfH2O = fH2O/fH2Opure
        fCO2pure = fco2calc(TC,P)
        fCO2trial = fCO2pure*(1-XfH2O)
        fdiff = fCO2 - fCO2trial
    
    #if previous pressure step yielded a fugacity difference closer to 0, use that pressure 
    if fdiffold < abs(fdiff):
        P = Pold
        XfH2O = XfH2Oold

    Press = P
    Fraction = round(XfH2O,3)
    return Press, Fraction

################################################################################
#read variables from csv file & perform calculations
################################################################################

#read inputs from a file
with open(satpressinput, 'r') as inputsatpressfile:
    csvinputfile = csv.reader(inputsatpressfile)
    n = -1
    #set up a list for outputs
    satpressinputlist = []
    for row in csvinputfile:
        n = n + 1
        if n > 1:
            satpressinputlist.append(row)

numsamples = len(satpressinputlist)
satpressoutputlist =[]
templist =[]
for m in range(numsamples):
    templist.clear()
    templist.append(satpressinputlist[m][0])
    for g in range(12):
        templist.append(float(satpressinputlist[m][g+1]))
    #calculate TD parameters
    deltaV_calc, lnK0_calc = TDparam(templist[3], templist[4], templist[5], templist[6], templist[8], templist[9], templist[10], templist[11])
    #calculate saturation pressure & fluid composition
    P_calc, F_calc = SatP(templist[1], templist[2], lnK0_calc, deltaV_calc)
    #add results to list for this composition
    templist.append(P_calc)
    templist.append(F_calc)
    #add full results to master list
    satpressoutputlist.append(templist[:])

################################################################################
#PRINT RESULTS
################################################################################

header = ['sample', 'dissolved H2O wt%', 'dissolved CO2 ppm', 'wt% SiO2','wt% TiO2','wt% Al2O3','wt% FeO','wt% MnO','wt% MgO','wt% CaO','wt% Na2O','wt% K2O','wt% P2O5','Pressure (bar)','XfH2O']
with open(satpressoutput, 'w', newline='') as outputcsv:
    csvwriter = csv.writer(outputcsv)
    csvwriter.writerow(header)
    csvwriter.writerows(satpressoutputlist)
