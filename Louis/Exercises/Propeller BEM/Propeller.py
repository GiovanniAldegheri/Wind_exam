#%%
import numpy as np
import os
import matplotlib.pyplot as plt
import math as m

# os.chdir(r'G:\Other computers\Grote Laptop\Desktop\TU Delft\MSc EWEM 1\Q1-2 DTU\45300 Wind turbine technology and aerodynamics\Exercises\Propeller')
path = 'Louis/Exercises/Propeller BEM/'


# Blade Geometry
data = np.loadtxt(path+'data.txt', skiprows = 1)
r = data[:, 0] * 10e-4          #m
r = r.round(5)
twist = np.radians(data[:, 1])  #rad
chord = data[:, 2] * 10e-4      #m
B = 2
pitch = 0

#Constants
P = 22e3
RPM = 2600
w = 2 * m.pi * RPM / 60 #rad/s
R = 0.690         #m
rho = 1.225

Vo = np.arange(1,41,1)  #m/s

#BEM LOOP
def BEM(r, chord, twist, Vo):
    a = 0
    aprime = 0
    convergenceFactor = 1e-7
    delta = 1
    deltaPrime = 1
    relax = 0.1
    solidity = (B * chord) / (2 * m.pi * r)
    count = 0

    while(delta > convergenceFactor and deltaPrime > convergenceFactor):
        count = count + 1
        if (count > 10000):
            print("No convergence!")
            break

        flowAngle = np.arctan(Vo * (1 + a) / (w * r * (1 - aprime)))     #DIFFERENT FROM NORMAL BEM
        localalpha =  - flowAngle + pitch + twist                        #DIFFERENT FROM NORMAL BEM
        Cl = 0.1 * np.degrees(localalpha) + 0.4
        Cd = 0.008

        Ct = Cl*m.sin(flowAngle) + Cd * m.cos(flowAngle)            #DIFFERENT FROM NORMAL BEM
        Cn = Cl*m.cos(flowAngle) - Cd * m.sin(flowAngle)            #DIFFERENT FROM NORMAL BEM


        F = 2 / np.pi * np.arccos(np.exp(-B / 2 * (R - r) / (r * np.sin(np.abs(flowAngle)))))
        # CT = ((1 - a)**2 * Cn *solidity) / (np.sin(flowAngle))**2

        aold = a
        aStar = (solidity * Cn * (1 + aold)) / (4 * F * (np.sin(flowAngle))**2)        #DIFFERENT FROM NORMAL BEM
        a = relax * aStar + (1 - relax) * aold
        delta = abs(a - aold)

        aprimeOld  = aprime
        aprimeStar = (solidity * Ct * (1 - aprimeOld)) / (4 * F * np.sin(flowAngle) * np.cos(flowAngle))      #DIFFERENT FROM NORMAL BEM
        aprime = relax * aprimeStar + (1 - relax) * aprimeOld
        deltaPrime = abs(aprime - aprimeOld)


    Vrel = m.sqrt(Vo**2 + (w * r)**2)
    Pn = 0.5 * rho * Vrel**2 * chord * Cn       
    Pt = 0.5 * rho * Vrel**2 * chord * Ct       
    
    if (m.isnan(Pt)|(m.isnan(Pn))):
        Pt, Pn = 0, 0

    return Pn, Pt

P = np.zeros(len(Vo))
T = np.zeros(len(Vo))

Cp = np.zeros(len(Vo))
Ct = np.zeros(len(Vo))

for ii in range(len(Vo)):
    Pn_lst = np.zeros(len(r))
    Pt_lst = np.zeros(len(r))

    for jj in range (len(r)-1):

        Pn, Pt = BEM(r[jj], chord[jj], twist[jj], Vo[ii])
        Pn_lst[jj] = Pn
        Pt_lst[jj] = Pt

    T[ii] = np.trapz(Pn_lst, r) * B
    P[ii] = np.trapz(Pt_lst * r, r) * w * B

    Cp[ii] = P[ii] / (0.5 * rho * Vo[ii]**3 * m.pi * R**2)
    Ct[ii] = T[ii] / (0.5 * rho * Vo[ii]**2 * m.pi * R**2)


plt.figure()
plt.plot(Vo, P)
plt.show()
# plt.plot(Vo, T)