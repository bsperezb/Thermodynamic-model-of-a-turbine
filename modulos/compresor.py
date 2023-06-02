import  CoolProp.CoolProp as CP


def compresor(P_a1, T_a1, RP_comp, N_c):
    """"
    Este modulo calcula las propiedades de salida de un compresor
    donde sus valores de entrada son:

    *P = presipon a la entrada.
    *T = Temperatura a la entrada.
    *RP_comp = razón de compresion entre la entrada y la salida.
    *N_c = eficiencia isentrópica.

    """

    h_a1=CP.PropsSI('H','P',P_a1,'T',T_a1,'AIR')
    s_a1=CP.PropsSI('S','P',P_a1,'T',T_a1,'Air')

    s_sa2 = s_a1
    P_a2= RP_comp*(P_a1)

    T_sa2=CP.PropsSI('T','S',s_sa2,'P',P_a2,'AIR')
    h_sa2=CP.PropsSI('H','T',T_sa2,'P',P_a2,'Air')

    h_a2=h_a1+((h_sa2-h_a1)/(N_c))  #Calculo de entalpia real usando la eficiencia del compresor
    T_a2=CP.PropsSI('T','H',h_a2,'P',P_a2,'Air')
    s_a2=CP.PropsSI('S','H',h_a2,'P',P_a2,'Air')

    return { "T":T_a2, "Ts":T_sa2, "P":P_a2, "h":h_a2, "h_s":h_sa2, "s":s_a2, "s_s":s_sa2 }
