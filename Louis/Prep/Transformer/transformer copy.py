import math as m
import numpy as np
import matplotlib.pyplot as plt
#CONSTANTS

f = 60
omega = 2*m.pi*f
V1_nominal = 10e3
V2_nominal = 0.69e3
RATIO = V1_nominal/V2_nominal

def bring_to_secondary(x): #works only for impedences
    return x/RATIO**2

def bring_to_primary(x): #works only for impedences
    return x*RATIO**2


def line_to_line(V): #brings single phase voltage to line to line
    return V*np.sqrt(3)
    # return I*np.sqrt(3)

def line_to_neutral(V):
    return V/np.sqrt(3)
    # return I/np.sqrt(3)

def print_list(x, y,x_label, y_label):
    print(x_label,'\t',y_label)
    for i in range(len(x)):
        print(x[i],'\t',y[i])
    plt.figure()
    plt.plot(x, y.real, label='real')
    if np.iscomplex(y[0]):
        plt.plot(x, y.imag, label='imag')
    plt.legend()
    plt.show()


complex_impedence_given = True
if complex_impedence_given:
    z1_p = (4.2 + 1j * 32) * 10**(-3)
    z2_s = (0.04 + 1j * 0.23) * 10**(-3)
    z2_p = bring_to_primary(z2_s)
else:
    R1_p = 4.3e-3
    R2_p = 0
    R2_s = 0.04e-3
    L1_p = 32e-3
    L2_p = 0
    L2_s = 0.23e-3
    ######
    z1_p = R1_p + L1_p * omega * 1j
    if L2_p == 0:
        z2_s = R2_s + 1j * omega * L2_s
        z2_p = bring_to_primary(z2_s)
    else:
        z2_p = R2_p + L2_p * omega  *1j

#CORE IMPEDENCE REFERRED TO PRIMARY
RC_p = 0
LM_p = 0
RC_s = 0    #core resistance
LM_s = 0    #core magnetizing inductance
ze = 786 + 1j * 258
if ze == 0:
    if RC_p != 0:
        ze = 1 / (1/RC_p + 1 / (1j*omega*(LM_p)))
    else:
        ze = 1 / (1/bring_to_primary(RC_s) + 1 / (1j*omega*bring_to_primary(LM_s))) #core impedence

print('z1 =', z1_p, '\nz2_p =', z2_p, '\nze =', ze)


V1 = line_to_neutral(10e3)

load = False
if load:
    Rl_s = 0.16
    # Ll_s = np.arange(500e-6,400e-5,1e-6)
    Ll_s = 460e-6
    zl_p = bring_to_primary(Rl_s+1j*omega*Ll_s)

    Ie = V1 / ze
    Il = V1 / (z1_p + z2_p + zl_p) - Ie
    Vl = Il * zl_p

    # V2_p = V1 - I1*z1 - I2_p*(R2_p+1j*omega*L2_p)

    power_in = V1 * (Il + Ie)
    power_out = zl_p * Il
    power_factor = power_out.real / np.abs(power_out)

    efficiency = power_out.real / power_in.real

    print(efficiency)
    # print_list(Ll_s,Vl,'zl_p','power factor')

else:
    I1_rms = np.arange(5,76)
    I1_module = I1_rms*np.sqrt(2)
    power_factor = 0.8

    inductive = True            ###########current lagging
    if inductive:
        theta = -np.arccos(power_factor)
    else:
        theta = np.arccos(power_factor)

    I1 = I1_module*(np.cos(theta)+1j*np.sin(theta))
    power_in = V1*I1

    Ve_p = V1 - I1*z1_p
    Ie = Ve_p / ze

    I2_p = I1-Ie
    V2_p = Ve_p - I2_p*z2_p
    power_out = V2_p*I2_p

    efficiency = power_out.real / power_in.real

    print_list(I1_rms, efficiency, 'i1','power out')


