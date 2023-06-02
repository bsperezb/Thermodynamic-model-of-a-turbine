from modulos.compresor import compresor
from modulos.camara_combustion import Camara_Combustion
import CoolProp.CoolProp as CP



# variables de entrada y constantes (por ahora):

                        #DATOS DE ENTRADA


#"Aire"
T_a1 = 30+ 273.15    #"Temperatura del Aire a la entrada del Compresor"
P_a1 =101325         #Presión del Aire a la entrada del Compresor"
m_dot_aire =0.304  #0.3004 - 582645    #kg/s
#"METANO"
T_CH4 = 25 + 273.15   #Temperatura del Metano CH4"
P_CH4 = 370      #3.7    Presión del Meano CH4"
m_dot_CH4 =0.006        #0.002  - 13.645        #kg/s

#RELACION DE COMPRESIÓN

RP_comp = 3.6                       #20  Relación de Presión en el Compresor"
RP_turb =3.6                       # 20  Relación de Presión en la Turbina"

#EFICIENCIAS
N_c = 0.78                           #Eficiencia adiabática Compresor (eta_com)
N_turb = 0.83                      #Eficiencia adiabática Turbina
eta_ph = 0.70                        #Eficiencia Intercambiador de Calor Precalentador
eta_cc = 0.98                        #Eficiencia Cámara de Combustión
eta_HRSG = 0.80                      #Eficiencia HRSG

m_dot_comb = m_dot_aire + m_dot_CH4        #Flujo masico total

#Diferencia estado de referencia coolprop con ecuación. AIR
#Vref_h=125850                                                             #+-0.1(entalpia)

#---------------(Propiedades de entrada)--------------------------

h_a1=CP.PropsSI('H','P',P_a1,'T',T_a1,'AIR')
s_a1=CP.PropsSI('S','P',P_a1,'T',T_a1,'Air')

#-------------(Calcuos del estado dos [compresor])---------------

dic_comp = compresor(P_a1, T_a1, RP_comp, N_c)
e2 = dic_comp
# print(dic_comp)
# print(AFR_gr)
# print(f"El AFR molar : {AFR_mol}")
#-------(Calculos del estado tres [Camara de Combustion])-------

print(e2["T"], e2["P"] )
Camara_Combustion( e2["T"], e2["P"], m_dot_aire, m_dot_CH4)

#Camara_Combustion(T_a2, P_a2, mM_aire, n_CH4, n_aire, PC_CH4_mol)

