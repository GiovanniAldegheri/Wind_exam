import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

'''
Code currently based on exams 2020 & 2022
'''
#CONSTANTS
f = 50
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

def line_to_neutral(V):
    return V/np.sqrt(3)

def print_list(x, y,x_label, y_label):
    print(x_label,'\t',y_label)
    for i in range(len(x)):
        print(x[i],'\t',y[i])
    plt.figure()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid()
    plt.plot(x, y.real, label='real')
    if np.iscomplex(y[0]):
        plt.plot(x, y.imag, label='imag')
    plt.legend()
    plt.show()

def complex_equations(vars, V1, z1, z2, zm, zl):
    # Unpack the variables
    I1_real, I1_imag, Im_real, Im_imag, I2_real, I2_imag = vars
    
    # Convert to complex numbers
    I1 = I1_real + 1j*I1_imag
    Im = Im_real + 1j*Im_imag
    I2 = I2_real + 1j*I2_imag
    
    # Define the equations
    eq1 = I1 - Im - I2
    eq2 = V1 - I1 * z1 - Im * zm
    eq3 = Im * zm - I2 * z2 - I2 * zl

    # Split the equations into real and imaginary parts (fsolve doesn't work with complex numbers)
    return [
        eq1.real, eq1.imag, 
        eq2.real, eq2.imag, 
        eq3.real, eq3.imag, 
    ]

def solve_system(V1, z1, z2, zm, zl, initial_guess):

    # Use fsolve to find the roots
    solution = fsolve(complex_equations, initial_guess, args=(V1, z1, z2, zm, zl))
    
    # Extract the solutions
    I1 = solution[0] + 1j*solution[1]
    Im = solution[2] + 1j*solution[3]
    I2 = solution[4] + 1j*solution[5]
    
    initial_guess=[solution[0],solution[1],solution[2],solution[3],solution[4],solution[5]]    

    return I1, Im, I2, initial_guess
    
    # # Rotate the solution so that V1 has an angle of 0 radians
    # delta = m.atan(V1.imag/V1.real)

    # V1_rotated = V1 * np.exp(-1j*delta)
    # V2_p_rotated = V2_p * np.exp(-1j*delta)
    # Vpoc_p_rotated = Vpoc_p * np.exp(-1j*delta)
    # I1_rotated = I1 * np.exp(-1j*delta)
    # I2_p_rotated = I2_p * np.exp(-1j*delta)
    # Ie_rotated = Ie * np.exp(-1j*delta)
    # Ipoc_p_rotated = Ipoc_p * np.exp(-1j*delta)
    
    # return V1_rotated, V2_p_rotated, Vpoc_p_rotated, I1_rotated, I2_p_rotated, Ie_rotated, Ipoc_p_rotated, initial_guess


#Calculate z1 and z2 primary
complex_impedence_given = True
if complex_impedence_given:
    R1_p = 3.52       #Primary Resistance
    X1_p = 4.32        #Primary Leakage Reactance

    R2_s = 0.04      #Secondary Resistance (Measured from secondary side)
    X2_s = 0.042      #Secondary Leakage Reactance (Measured from secondary side)
    
    z1_p = R1_p + 1j * X1_p
    z1_s = bring_to_secondary(z1_p)
    z2_s = R2_s + 1j * X2_s
    z2_p = bring_to_primary(z2_s)

else:
    R1_p = 0        #Primary Resistance
    R2_p = 0        #Secondary Resistance
    R2_s = 0        #Secondary Resistance (Measured from secondary side)
    L1_p = 0        #Primary Inductance
    L2_p = 0        #Secondary inductance
    L2_s = 0        #Secondary Inductance (Measured from secondary side)

    z1_p = R1_p + L1_p* omega * 1j
    if L2_p == 0:
        z2_s = R2_s + 1j * omega * L2_s
        z2_p = bring_to_primary(z2_s)
    else:
        z2_p = R2_p + L2_p*omega*1j
        z2_s = bring_to_secondary(z2_p)

#Calculate zm primary
RC_p = 0            #Core Resistance
XM_p = 0            #Core Magnetizing Reactance 
RC_s = 20.5          #Core Resistance (Measured from secondary side)
XM_s = 0.05*omega          #Core Magnetizing Reactance (Measured from secondary side)

if RC_p == 0:
    zm_s = 1 / ((1 / RC_s) + (1 / (XM_s * 1j)))
    zm_p = bring_to_primary(zm_s)
else:
    zm_p = 1/((1/RC_p) + (1 / (XM_p * 1j)))
    zm_s = bring_to_secondary(zm_p)

V1_p = line_to_neutral(10e3)
V1_s = V1_p / RATIO

print('z1_p =', z1_p, '\nz2_p =', z2_p, '\nzm_p =', zm_p)
print('\nz1_s =', z1_s, '\nz2_s =', z2_s, '\nzm_s =', zm_s)

load = False
if load == True:
    #Calculate zl primary
    Rl_s = 0.16
    Ll_s = 460e-6

    power_factor = np.cos(np.arctan(Ll_s*omega/Rl_s))

    zl_s = Rl_s + 1j *omega * Ll_s
    zl_p = bring_to_primary(zl_s)

    #Solve system of equations
    initial_guess= [0,0,0,0,0,0]
    I1_p, Im_p, I2_p, initial_guess = solve_system(V1_p, z1_p, z2_p, zm_p, zl_p, initial_guess)

    initial_guess= [0,0,0,0,0,0]
    I1_s, Im_s, I2_s, initial_guess = solve_system(V1_s, z1_s, z2_s, zm_s, zl_s, initial_guess)
    #You could also do I_s = I_p * RATIO but it iz what it iz

    Exam2020 = False
    if Exam2020 == True:
        # What is the (line to line) voltage at the transformer secondary?
        V2_p = line_to_line(I2_p * zl_p)
        V2_s = line_to_line(I2_s * zl_s)
        print('Voltage at transformer secondary =',V2_s)
        if V2_s.round(4) == (V2_p / RATIO).round(4):
            print('\tSanity check pass - V LV = V HV * RATIO')
        else:
            print('\tSanity check fail - V LV != V HV * RATIO')

        # What are active and reactive power at HV of the transformer?
        S1_p = V1_p * np.conj(I1_p) * 3
        S1_s = V1_s * np.conj(I1_s) * 3
        print('\nComplex power at primary =',V2_p)
        if S1_s.round(4) == S1_p.round(4):
            print('\tSanity check pass - Power constant LV / HV side')
        else:
            print('\tSanity check fail - Power not constant LV / HV side')

        # What is power factor of the load? How can this be improved?
        PF = np.cos(np.angle(V2_s) - np.angle(I2_s))
        print('\n Power Factor =', PF.round(3))

        # What the efficiency of the transformer?
        S2_p = line_to_neutral(V2_p) * np.conj(I2_p) * 3
        S2_s = line_to_neutral(V2_s) * np.conj(I2_s) * 3
        eff = S2_p.real / S1_p.real
        print('\n Efficiency =', eff.round(3))

else:
    V1_p = line_to_neutral(10e3)

    S_nominal = 1.732e6 / 3

    I1_module = np.arange(1,300,0.1)
    # I1_module = I1_rms*np.sqrt(2)

    power_factor = 0.9
    inductive = True            ###########current lagging the voltage
    if inductive:
        theta = -np.arccos(power_factor)
    else:
        theta = np.arccos(power_factor)

    I1 = I1_module*(np.cos(theta)+1j*np.sin(theta))
    power_in = V1_p * I1

    Vm_p = V1_p - I1 * z1_p

    Im_p = Vm_p / zm_p

    I2_p = I1 - Im_p

    V2_p = Vm_p - I2_p * z2_p

    power_out = V2_p * I2_p

    efficiency = power_out.real / power_in.real
    
    count = 0
    for i in range(len(I1_module)):
        if np.abs(power_out[i]) >= S_nominal and count == 0:
            index = i
            print('1_1 =', I1_module[index])
            count += 1

    
    print_list(I1_module, efficiency, r'$|I_1|$ [A]','Efficiency')

