P = 4000000
poles = 52
w_rotor = 22.5 #rpm
flux_nom = 15.826 #Wb
Ra = 14.821 #ohm
jXs = 5.573 #mH

f = w_rotor*poles/120
w = 3.1456*f

Ea = w*flux_nom
Ia = P/(3*Ea)
Va = Ea - Ia*(Ra)
Va_i = -Ia*jXs

print(Ea, Va, Va_i)

