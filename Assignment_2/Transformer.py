import math as m
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


# CONSTANTS

#Transformer
f = 50
omega = 2*m.pi*f
V1_rated = 690
Vpoc = 33000
RATIO = V1_rated/Vpoc
L1 = 20e-6
L2_p = 20e-6
R1 = 10e-3
R2_p = 10e-3
RC = 124
LM = 25e-3

P_real = 3865000  # Active power from generator, 3 phases
P_im = P_real*0.2
P = P_real + P_im*1j
P_cmp = m.sqrt(P_real**2 + P_im**2)
delta = m.acos(P_real/P_cmp)
PF = (P_real/P_cmp)

S = P/3

#Cable
lenght = 60
r_cable = 0.5*lenght
l_cable = 0.5e-3*lenght
c_cable = 0.1e-6*lenght

######
def bring_to_secondary(x):
    return x/RATIO**2

def bring_to_primary(x):
    return x*RATIO**2


z1 = R1 * L1*omega*1j
ze = 1 / (1/RC + 1 / (1j*omega*LM))
zcomp_p = R2_p + L2_p*omega*1j + bring_to_primary(r_cable) + bring_to_primary(l_cable*omega*1j)
zc_p = bring_to_primary(1/(c_cable*omega*1j))


Vpoc_p = Vpoc/m.sqrt(3)*RATIO
print(Vpoc_p)


Ic_p = Vpoc_p/zc_p
print(Ic_p)

def equations(vars, ):
    V1, I1, Ie, I2_p, Ipoc_p = vars

    eq1 = V1 - I1*z1 - Ie*ze
    eq2 = Ie*ze - I2_p*zcomp_p - Ic_p*zc_p
    eq3 = I1 - Ie - I2_p
    eq4 = I2_p - Ic_p - Ipoc_p
    eq5 = S - I1*V1

delta = atan(V1.im/v1.real)