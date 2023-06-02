import CoolProp.CoolProp as CP
from sympy import *



#-----------------------------------------------------
# Revisar mM_aire
def Camara_Combustion(T_a2, P_a2, m_dot_aire, m_dot_CH4):
    """
    PC_CH4_mol = poder calirifico del combustible
    AFR_estq = corresponde a la razon aire combustible molar (moles-Aire/moles-Combustible)
    AFR_gr = corresponde al  la  razón aire combustibel grabimetrico (masa-Aire / masa.Combustible)

    """
    #------------------------Calculo de flujos de aire y combustible------------------------
    #--------------- variables de entrada------
    mM_aire=28.96
    mM_CH4=16.043
    mM_CO2=44.01
    mM_O2=16
    mM_H2O=18.01
    mM_N2=28.01
    PC_CH4_grav=55.050                       #[kJ/kg]
    AFR_estq=9.52                           #Molar (valor constante )
    #GASES IDEALES
    R=8314                     #constante de los gases [j/kmol*K]
    T_ref=25+273.15            #temperatura de referencia [K]
    P_ref=101325                #Presión de referencia [Pa]
    #----------


    m_dot_comb = m_dot_aire + m_dot_CH4        #Flujo masico total

    #AFR
    AFR_gr=(m_dot_aire/m_dot_CH4)           #grabimétrico
    AFR_mol=(m_dot_aire/m_dot_CH4)*(mM_CH4/mM_aire) #Molar
    #RELACION DE EQUIVALENCIA  DE AIRE(exceso)
    R_eq=(AFR_mol/AFR_estq)

    #NUMERO DE MOLES (reactivos y productos)-------
    #en esta seccion se calcula la cantidad de moles de los productos dependiendo de los reactivos

    n_CH4=1
    n_aire=R_eq*(2*(1+79/21))
    n_CO2=1
    n_H2O=2                    #a verificar (varia dependiendo de la humedad del aire).
    if (R_eq>1):                #si es mayor a la estequimétrica
        n_N2=(2*(79/21))+2*(R_eq-1)*(79/21)
        n_O2=2*(R_eq-1)
    else:
        n_N2=(2*(79/21))
        n_O2=0


    n_productos=n_CO2+n_H2O+n_N2+n_O2       #Moles totales productos
    Fr_CO2=n_CO2/n_productos                #fracciones molares de los productps
    Fr_H2O=n_H2O/n_productos                #para el calculo de propiedades
    Fr_N2=n_N2/n_productos
    Fr_O2=n_O2/n_productos

    mM_prod=Fr_CO2*mM_CO2+Fr_O2*mM_O2+Fr_H2O*mM_H2O+Fr_N2*mM_N2
    m_prod_total=mM_CO2*n_CO2+mM_O2*n_O2+mM_H2O*n_H2O+mM_N2*n_N2

    FrM_CO2=(mM_CO2*n_CO2)/m_prod_total    #fracciones molares de los productps
    FrM_H2O=(mM_H2O*n_H2O)/m_prod_total    #para el calculo de propiedades
    FrM_N2=(mM_N2*n_N2)/m_prod_total
    FrM_O2=(mM_O2*n_O2)/m_prod_total

    #PODER CALORIFICO DE COMBUSTIBLES (inferior)
    PC_CH4_mol=PC_CH4_grav*mM_CH4                                             #[kJ/kmol]
    PC_CH4_mol=PC_CH4_mol*1000                                                #[J/kmol]

    #CAPACIDAD CALORIFICA A PRESION CONSTANTE(Cp)[kJ/kmol*K]
    #aire
    A_aire=28.11
    B_aire=0.1967*10**-2
    C_aire=0.4802*10**-5
    D_aire=1.966*10**-9
    #N2
    A_N2=28.90
    B_N2=-0.1571*10**-2
    C_N2=0.8081*10**-5
    D_N2=-2.873*10**-9
    #O2
    A_O2=25.48
    B_O2=1.520*10**-2
    C_O2=-0.7155*10**-5
    D_O2=1.312*10**-9
    #CO2
    A_CO2=22.26
    B_CO2=5.981*10**-2
    C_CO2=-3.501*10**-5
    D_CO2=7.469*10**-9
    #H2O (vapor)
    A_H2O=32.24
    B_H2O=0.1923*10**-2
    C_H2O=1.055*10**-5
    D_H2O=-3.595*10**-9




    #---------------------------------------------------
    P_prod3=P_a2
    h_a2=CP.PropsSI('H','T',T_a2,'P',P_a2,'Air')
    h_a2_ref=h_a2-(CP.PropsSI('H','T',298.15,'P',P_a2,'Air'))
    h_a2_ref=(h_a2_ref/mM_aire)

    # baleance de enrgia en la camara de combustion

    IN_CC=(n_CH4)*(PC_CH4_mol)+(h_a2_ref)*(n_aire)   #[J/Kmol]+[J/kmol]
    T_prod3 = symbols("T_prod3")

    Cp_CO2=A_CO2+B_CO2*T_prod3+C_CO2*T_prod3**2+D_CO2*T_prod3**3
    Cp_H2O=A_H2O+B_H2O*T_prod3+C_H2O*T_prod3**2+D_H2O*T_prod3**3
    Cp_O2=A_O2+B_O2*T_prod3+C_O2*T_prod3**2+D_O2*T_prod3**3
    Cp_N2=A_N2+B_N2*T_prod3+C_N2*T_prod3**2+D_N2*T_prod3**3

    niCpi_prom3=Cp_H2O*n_H2O+Cp_CO2*n_CO2+Cp_O2*n_O2+Cp_N2*n_N2     #[kJ/kmol*K]
    OUT_CC = integrate(niCpi_prom3,(T_prod3, (298.5, T_prod3)))

    print(f"OUT_cc: {OUT_CC}")
    print(f"in_cc : {IN_CC}")

    BE_CC = OUT_CC - IN_CC
    BE_CC = Eq(BE_CC, 0)
    print(f"BE_CC ecuation : {BE_CC}")
    print(BE_CC, type(BE_CC))
    T_prod3 = solveset(BE_CC, T_prod3, domain=S.Reals)
    print(AFR_gr)
    return print(T_prod3.args, type(T_prod3.args))

