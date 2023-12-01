import math as m
import numpy as np
from sympy.solvers import solve
import sympy as sym
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

P_real = 3915910  # Active power from generator, 3 phases
P_im = P_real*0.2
P = P_real + P_im*1j
# P_cmp = m.sqrt(P_real**2 + P_im**2)
# delta = m.acos(P_real/P_cmp)
# PF = (P_real/P_cmp)

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

def rotate(var, delta):
    var = complex(var)
    x = var.real
    y = var.imag
    if x == 0:
        angle = m.pi*m.sign(y)/2
    else:
        angle = m.atan(y/x)
    module = m.sqrt(x**2+y**2)
    x = module*m.cos(angle+delta)
    y = module*m.sin(angle+delta)

    return complex(x,y)

z1 = R1 + L1*omega*1j
ze = 1 / (1/RC + 1 / (1j*omega*LM))
zcomp_p = R2_p + L2_p*omega*1j + bring_to_primary(r_cable) + bring_to_primary(l_cable*omega*1j)
zc_p = bring_to_primary(1/(c_cable*omega*1j))

Vpoc_p = Vpoc/m.sqrt(3)*RATIO

Ic_p = Vpoc_p/zc_p

## Set up system of equations
V1, I1, Ie, I2_p, Ipoc_p = sym.symbols("V1, I1, Ie, I2_p, Ipoc_p")


def transformer(S):
    eq1 = sym.Eq(V1 - I1*z1 - Ie*ze,0)
    eq2 = sym.Eq(Ie*ze - I2_p*zcomp_p - Ic_p * zc_p,0)
    eq3 = sym.Eq(I1 - Ie - I2_p,0)
    eq4 = sym.Eq(I2_p - Ic_p - Ipoc_p,0)
    eq5 = sym.Eq(np.conj(S)-I1*V1,0)

    ## Solve system of equations
    solution = solve([eq1,eq2,eq3,eq4,eq5],(V1, I1, Ie, I2_p, Ipoc_p))
    V1, I1, Ie, I2_p, Ipoc_p = (solution[1])
    V1 = np.conj(complex(V1))

    ## Transfer everything so that V1 = 0 rad
    delta = m.atan(V1.imag/V1.real)

print(rotate(V1,-delta))
print(rotate(I1,-delta))
print(rotate(I2_p,-delta))
