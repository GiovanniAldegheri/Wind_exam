import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os


# os.chdir(r'data')
# os.chdir(r'G:\Other computers\Grote Laptop\Desktop\TU Delft\MSc EWEM 1\Q1-2 DTU\45300 Wind turbine technology and aerodynamics\Shared Git\Wind_exam\Louis\Exercises\HAWT BEM\data')

files=['FFA-W3-301.txt','FFA-W3-301.txt','FFA-W3-360.txt','FFA-W3-480.txt','FFA-W3-600.txt','cylinder.txt']
#Initializing tables    
cl_tab=np.zeros([105,6])
cd_tab=np.zeros([105,6])
cm_tab=np.zeros([105,6])
aoa_tab=np.zeros([105,])
#Readin of tables. Only do this once at startup of simulation
for i in range(np.size(files)):
    aoa_tab[:],cl_tab[:,i],cd_tab[:,i],cm_tab[:,i] = np.loadtxt(files[i], skiprows=0).T

thick_prof=np.zeros(6)
thick_prof[0]=24.1;
thick_prof[1]=30.1;
thick_prof[2]=36;
thick_prof[3]=48;
thick_prof[4]=60;
thick_prof[5]=100;

bladedat = pd.read_csv('bladedat.txt',sep="\t", header=None)
r_ref = bladedat[0].tolist() #m
r_ref.pop()
r_ref = np.array(r_ref)
c_ref = bladedat[2].tolist() #m
beta_ref = bladedat[1].tolist() #deg
tc_ref = bladedat[3].tolist() #%



#Functions____________
def contourplots(pitch, TSR, Cp, Ct):
    # Create a figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns of subplots

    # Subplot 1
    axs[0].set_title(r'$C_p(\lambda,\theta_p)$ Contour Plot')
    [X, Y] = np.meshgrid(pitch, TSR)
    cont1 = axs[0].contourf(Y, X, Cp, 60, cmap="turbo")
    axs[0].set_ylabel('Pitch angle (deg)')
    axs[0].set_xlabel('Tip speed ratio (-)')
    cbar1 = plt.colorbar(cont1, ax=axs[0])
    cbar1.set_label(r'$C_p$')

    # Subplot 2
    axs[1].set_title(r'$C_t(\lambda,\theta_p)$ Contour Plot')
    [X, Y] = np.meshgrid(pitch, TSR)
    cont2 = axs[1].contourf(Y, X, Ct, 60, cmap="turbo")
    axs[1].set_ylabel('Pitch angle (deg)')
    axs[1].set_xlabel('Tip speed ratio (-)')
    cbar2 = plt.colorbar(cont2, ax=axs[1])
    cbar2.set_label(r'$C_t$')

    plt.tight_layout()
    plt.savefig('plots/Glauert.png')
    plt.close()

def force_coeffs(localalpha,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    cl_aoa=np.zeros([1,6])
    cd_aoa=np.zeros([1,6])
    cm_aoa=np.zeros([1,6])
    
    #Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0,i]=np.interp (localalpha,aoa_tab,cl_tab[:,i])
        cd_aoa[0,i]=np.interp (localalpha,aoa_tab,cd_tab[:,i])
        cm_aoa[0,i]=np.interp (localalpha,aoa_tab,cm_tab[:,i])
    
    #Interpolate to current thickness:
    Cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    Cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    Cm=np.interp (thick,thick_prof,cm_aoa[0,:])
    return Cl, Cd, Cm 

def BEM(omega,pitch,r,c,twist,thick,aoa_tab,cl_tab,cd_tab,cm_tab,Vo):
    TSR = omega * R / Vo
    a = 0
    aprime = 0
    convergenceFactor = 1e-10
    delta = 1
    deltaPrime = 1

    relax = 0.1
    solidity = (B*c)/(2*m.pi*r)
    count = 0

    while(delta > convergenceFactor and deltaPrime > convergenceFactor):
        count = count + 1
        if (count > 1e4):
            print("No convergence!")
            break

        flowAngle = m.atan(((1 - a) * R) / ((1 + aprime) * TSR * r))
        localalpha =  m.degrees(flowAngle) - (pitch + twist)

        Cl,Cd,Cm = force_coeffs(localalpha, thick, aoa_tab, cl_tab, cd_tab, cm_tab)
        Ct = Cl*m.sin(flowAngle) - Cd*m.cos(flowAngle)
        Cn = Cl*m.cos(flowAngle) + Cd*m.sin(flowAngle)

        F = 2 / m.pi * m.acos(m.exp(-B * (R-r) / (2 * r *m.sin(abs(flowAngle)))))
        CT = ((1 - a)**2 * Cn * solidity) / m.sin(flowAngle)**2

        aold = a
        if(aold < 0.33):
            a = (solidity * Cn * (1 - aold)) / (4 * F * m.sin(flowAngle)**2)
        else:
            aStar = CT / (4 *F*(1 - 1/4 * (5 - 3 * aold) * aold))
            a = relax * aStar + (1 - relax) * aold

        aprimeOld  = aprime
        aprimeStar = (solidity * Ct * (1 + aprimeOld)) / (4 * F * m.sin(flowAngle) * m.cos(flowAngle))
        aprime = relax * aprimeStar + (1 - relax) * aprimeOld

        delta = abs(aprime - aprimeOld)
        deltaPrime = abs(aprime - aprimeOld)

    Vrel = m.sqrt(Vo**2+ (omega * r)**2)

    Pn = 0.5 * rho * Vrel**2 * c * Cn
    Pt = 0.5 * rho * Vrel**2 * c * Ct

    if (m.isnan(Pt)|(m.isnan(Pn))):
        Pt, Pn = 0,0

    return Pn, Pt

def single_BEM_loop(Vo):
    for k in range(len(r_ref)):
        Pn, Pt = BEM(omega,pitch,r_ref[k],c_ref[k],beta_ref[k],tc_ref[k],aoa_tab,cl_tab,cd_tab,cm_tab,Vo)
        # print(r_ref[k], Pn)
        Pn_lst[k] = Pn
        Pt_lst[k] = Pt

    T = np.trapz(Pn_lst, r_ref) * B
    P = np.trapz(Pt_lst * r_ref, r_ref) * omega * B

    Cp = P/(0.5*rho*Vo**3*m.pi*R**2)
    Ct = T/(0.5*rho*Vo**2*m.pi*R**2)
    return P, T, Cp, Ct, Pn_lst, Pt_lst

#Constants______________

R = 89.17 #m
B = 3
rho = 1.225 #kg/m3
# Vo = 15
Vo = np.arange(7,20,0.1)
Vo = np.round(Vo,1)

#Interpolate over r, tip speed ratio and pitch

# RPM = 1.5
omega = 1.5
# TSR = omega * R / Vo
# TSR = np.arange(RPM*R/Vo,RPM*R/Vo+1)

# pitch = np.arange(0,1)
pitch = 0

#Blade characteristics
P_max = 0
Cp_max = 0
TSR_max = 0
pitch_max = 0

# Cp=np.zeros([len(TSR),len(pitch)])
# Ct=np.zeros([len(TSR),len(pitch)])

Pn_lst = np.zeros(len(r_ref))
Pt_lst = np.zeros(len(r_ref))

def iterative_BEM_loop_pitch():
    pitch = np.arange(0,15,0.1)
    pitch = np.round(pitch,1)
    P = 0

    for ii, pitch in enumerate(pitch):
        P, T, Cp, Ct, Pn_lst, Pt_lst = single_BEM_loop(pitch)
        print('Pitch = ',pitch,'[deg] \t P = ',P,'[W]')
        if P <= 10e6:
            break

def iterative_BEM_loop_EXAM2023(x):
    P = 1
    Vo_store = []
    P_store = []
    count = 0
    for ii, Vo in enumerate(x):
        P, T, Cp, Ct, Pn_lst, Pt_lst = single_BEM_loop(Vo)
        print('Vo = ',Vo,'[m/s] \t P = ',P,'[W]')
        Vo_store.append(Vo)
        P_store.append(P)

        if P >= 10e6 and count == 0:
            v_0, P_0 = Vo, P
            count += 1
            
    print('V_0', v_0, '[m/s]', 'P_0', P_0 ,'[W]')
    plt.plot(Vo_store,P_store)
    plt.xlabel('Vo [m/s]')
    plt.ylabel('P [W]')
    plt.show()

def get_loads():
    P, T, Cp, Ct, Pn_lst, Pt_lst = single_BEM_loop()

    #When asked to plot deflection at certain conditions, copy output to loads_custom.txt and run deflection.py
    print('# Vo=',Vo,'m/s, pitch=',pitch,' deg, omega=',omega,' rad/s')
    print('#    r [m]   pn [kN/m]  pt [kN/m]')
    for i in range(len(Pn_lst)):
        print('  ',round(r_ref[i],4),'  ',round(Pn_lst[i]/1000,4),'  ',round(Pt_lst[i]/1000,4))
    print('   89.1660         0         0')

iterative_BEM_loop_EXAM2023(Vo)

# P =  22175125.080722693 [W]