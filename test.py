import CoolProp
import numpy as np

RP_comp = 3.6 
# Propiedades del aire
air_fluid = 'Air'
T_a1 = 30+ 273.15
P_a1 = 101325
m_dot_aire = 0.304
air = CoolProp.AbstractState('HEOS', air_fluid)
air.update(CoolProp.PT_INPUTS, P_a1, T_a1)
h_a1 = air.hmolar()

# Propiedades del metano
methane_fluid = 'Methane'
T_CH4 = 25 + 273.15
P_CH4 = 370
m_dot_CH4 = 0.006
methane = CoolProp.AbstractState('HEOS', methane_fluid)
methane.update(CoolProp.PT_INPUTS, P_CH4, T_CH4)
h_CH4 = methane.hmass() # Entalpía del metano

# Calculo de la masa de aire y de combustible
alpha = 1.5*(1+1/0.068)
m_dot_fuel = m_dot_CH4 / alpha
m_dot_total = m_dot_aire + m_dot_fuel

# Calor de combustión del metano
LHV_CH4 = 50.1e6 # J/kg
Q_comb = LHV_CH4 * m_dot_fuel

# Propiedades del agua
water_fluid = 'Water'
water = CoolProp.AbstractState('HEOS', water_fluid)

# Cálculo del ciclo Brayton en cogeneración
# Proceso 1-2 (Compresión adiabática del aire)
P_a2 = RP_comp * P_a1
air.update(CoolProp.PSmolar_INPUTS, P_a2, h_a1)
T_a2 = air.T() # Temperatura del aire a la salida del compresor
s_a2 = air.smass() # Entropía del aire a la salida del compresor
W_c = m_dot_total * (h_a1 - air.hmass()) / N_c # Trabajo realizado por el compresor

# Proceso 2-3 (Precalentamiento)
T_PH_in = T_a2
T_PH_out = T_PH_in + (Q_comb/(m_dot_total*Cp_air)) * eta_ph # Temperatura de salida del precalentador
Q_PH = m_dot_total * Cp_air * (T_PH_out - T_PH_in) # Calor transferido en el precalentador
P_PH_in = P_a2
P_PH_out = RP_turb * P_PH_in
water.update(CoolProp.PT_INPUTS, P_PH_in, T_PH_in)
h_w_in = water.hmass() # Entalpía del agua en la entrada del precalentador
s_w_in = water.smass() # Entropía del agua en la entrada del precalentador
water.update(CoolProp.PT_INPUTS, P_PH_out, T_PH_out)
h_w_out = water.hmass() # Entalpía del agua en la salida del precalentador
s_w_out = water.smass() # Entropía del agua en la salida del precalentador

# Proceso 3-4 (Expansión adiabática en la turbina)
air.update(CoolProp.PS_INPUTS, P_PH_out, s_a2)
T_a4 = air.T() # Temperatura del aire a la salida de la turbina
s_a4 = air.smass() # Entropía del aire a la salida de la turb