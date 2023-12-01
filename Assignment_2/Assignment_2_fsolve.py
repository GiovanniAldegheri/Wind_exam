import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math as m
#%%
np.seterr(divide='ignore', invalid='ignore')
"""
PLOTS
"""
def plot_efficiency(eff):
    a_lst = [78.92240667,79.98490763,80.95343249,81.8766464,82.76986694,83.63938648,84.4879822,85.3167367,86.12579758,86.91474514,87.6827903,88.42888978,89.15181695,89.85020669,90.52258332,91.16737593,91.78292271,92.36746436,92.91912482,93.43587656,93.91548561,94.35542952,94.75277793,95.10402088,95.40482214,95.64966358,95.83132717,95.94013082,95.96278129,95.88061784,95.66685848,95.28216566,94.66728276,93.73037011,0,0,0,0,0,0,0,0,0,0,0]
    plt.figure()
    plt.plot(RPM[:32], a_lst[:32], label = 'Case A')
    plt.plot(RPM[:32], eff[:32], label = 'Case B')
    plt.xlim(5, 25)
    plt.xlabel('RPM')
    plt.ylabel('Efficiency [%]')
    plt.legend()

def plots_1():
    fig, axs = plt.subplots(2, 2, figsize=(15, 10)) # 2 rows, 2 columns
    # Note, lists are reversed (22.5, 22, 21.5 ...) 32 is the index of start of 0 values (cut-in)
    # Plot 1: Current, Terminal-, Induced Voltage
    # axs[0, 0].set_title('Current, Terminal, Induced Voltage')
    axs[0, 0].grid(True, linestyle='--', alpha=0.5)
    ax2 = axs[0, 0].twinx()
    ax2.plot(RPM[:32], Ia_lst[:32], marker='o', markersize = 3, linestyle='--', label='Phase Current Amplitude (A)')
    axs[0, 0].plot([], [], marker='s', markersize = 3, linestyle='--', label='Phase Current Amplitude (A)')
    axs[0, 0].plot(RPM[:32], Va_lst[:32], marker='s', markersize = 3, linestyle='--', label='Terminal Voltage Amplitude (V)')
    axs[0, 0].plot(RPM[:32], Ea_lst[:32], marker='^', markersize = 3, linestyle='--', label='Induced Voltage Amplitude (V)')
    axs[0, 0].set_xlim(5, 25)
    axs[0, 0].set_xlabel('RPM')
    axs[0, 0].set_ylabel('Voltage [V]')
    ax2.set_ylabel('Current [A]')
    axs[0, 0].legend()

    # Plot 2: Power
    # axs[0, 1].set_title('Power')
    axs[0, 1].grid(True, linestyle='--', alpha=0.5)
    axs[0, 1].plot(RPM[:32], Sg.imag[:32], marker='o', markersize = 3, linestyle='--', label='Reactive Power')
    axs[0, 1].plot(RPM[:32], Sg.real[:32], marker='s', markersize = 3, linestyle='--', label='Active Power')
    axs[0, 1].plot(RPM[:32], Sg_lst[:32], marker='^', markersize = 3, linestyle='--', label='Mechanical Power')
    axs[0, 1].set_xlim(5, 25)
    axs[0, 1].set_xlabel('RPM')
    axs[0, 1].set_ylabel('Power (W)')
    axs[0, 1].legend()

    # Plot 3: Power Loss
    # axs[1, 0].set_title('Power Loss')
    axs[1, 0].plot(RPM, P_loss, marker='o', markersize = 3, linestyle='--', label='Power Loss')
    axs[1, 0].set_xlim(5, 25)
    axs[1, 0].set_xlabel('RPM')
    axs[1, 0].set_ylabel('Power (W)')
    axs[1, 0].grid(True, linestyle='--', alpha=0.5)
    axs[1, 0].legend()

    # Plot 4: Efficiency
    # axs[1, 1].set_title('Efficiency')
    axs[1, 1].plot(RPM, eff, marker='o', markersize = 3, linestyle='--', label='Generator Efficiency')
    axs[1, 1].set_xlim(5, 25)
    axs[1, 1].set_xlabel('RPM')
    axs[1, 1].set_ylabel('Efficiency (%)')
    axs[1, 1].grid(True, linestyle='--', alpha=0.5)
    axs[1, 1].legend()

    plt.tight_layout()

def plots_2():
    fig, axs = plt.subplots(2, 2, figsize=(15, 10)) # 2 rows, 2 columns

    # Plot 1: Voltage
    # axs[0, 0].set_title('Voltage')
    axs[0, 0].grid(True, linestyle='--', alpha=0.5)
    axs[0, 0].plot(RPM[:32], abs(V1_lst[:32]), marker='o', markersize = 3, linestyle='--', label='Transformer Primary Voltage')
    axs[0, 0].plot(RPM[:32], abs(V2_p_lst[:32]), marker='s', markersize = 3, linestyle='--', label='Transformer Secondary Voltage')
    axs[0, 0].plot(RPM[:32], abs(Vpoc_p_lst[:32]), marker='^', markersize = 3, linestyle='--', label='POC Voltage')
    axs[0, 0].set_xlim(5, 25)
    axs[0, 0].set_xlabel('RPM')
    axs[0, 0].set_ylabel('Voltage [V]')
    axs[0, 0].legend()

    # Plot 2: Power
    # axs[0, 1].set_title('Active Power')
    axs[0, 1].grid(True, linestyle='--', alpha=0.5)
    axs[0, 1].plot(RPM[:32], Sg.real[:32], marker='o', markersize = 3, linestyle='--', label='Transformer Primary Active Power')
    axs[0, 1].plot(RPM[:32], S2_lst.real[:32], marker='s', markersize = 3, linestyle='--', label='Transformer Secondary Active Power')
    axs[0, 1].plot(RPM[:32], Spoc_lst.real[:32], marker='^', markersize = 3, linestyle='--', label='POC Active Power')
    axs[0, 1].set_xlim(5, 25)
    axs[0, 1].set_xlabel('RPM')
    axs[0, 1].set_ylabel('Power (W)')
    axs[0, 1].legend()

    # Plot 3: Power Loss
    # axs[1, 0].set_title('Reactive Power')
    axs[1, 0].grid(True, linestyle='--', alpha=0.5)
    axs[1, 0].plot(RPM[:32], Sg.real[:32]*0.2, marker='o', markersize = 3, linestyle='--', label='Transformer Primary Reactive Power')
    axs[1, 0].plot(RPM[:32], S2_lst.imag[:32], marker='s', markersize = 3, linestyle='--', label='Transformer Secondary Reactive Power')
    axs[1, 0].plot(RPM[:32], Spoc_lst.imag[:32], marker='^', markersize = 3, linestyle='--', label='POC Reactive Power')
    axs[1, 0].set_xlim(5, 25)
    axs[1, 0].set_xlabel('RPM')
    axs[1, 0].set_ylabel('Power (W)')
    axs[1, 0].legend()

    # Plot 4: Efficiency
    # axs[1, 1].set_title('Efficiency')
    axs[1, 1].plot(RPM[:32], eff_transformer[:32], marker='o', markersize = 3, linestyle='--', label='Transformer Efficiency')
    axs[1, 1].plot(RPM[:32], eff_cable[:32], marker='s', markersize = 3, linestyle='--', label='Cable Efficiency')
    axs[1, 1].plot(RPM[:32], eff_transmission[:32], marker='^', markersize = 3, linestyle='--', label='Transmission Efficiency')
    axs[1, 1].set_xlabel('RPM')
    axs[1, 1].set_ylabel('Efficiency (%)')
    axs[1, 1].grid(True, linestyle='--', alpha=0.5)
    axs[1, 1].legend()

    plt.tight_layout()

# Function to update the plot for each frame (*start point, *vector direction)
def update(frame):
    plt.clf()
    plt.xlabel('Real Axis')
    plt.ylabel('Imaginary Axis')
    plt.quiver(*(0, 0), *(Va[frame].real, Va[frame].imag), angles='xy', scale_units='xy', scale=1, width = 0.005, color='r', label = 'Va')
    plt.quiver(*(Va[frame].real, Va[frame].imag), *(Ra*Ia[frame].real, Ra*Ia[frame].imag), angles='xy', scale_units='xy', scale=1, width = 0.005, color='b', label = 'Ia Ra')
    plt.quiver(*(0, 0), *(Ea[frame].real, Ea[frame].imag), angles='xy', scale_units='xy', scale=1, width = 0.005, color='g', label = 'Ea')
    plt.quiver(*(Va[frame].real+Ra*Ia[frame].real, Va[frame].imag+Ra*Ia[frame].imag), *(-Ls*omega[frame]*Ia[frame].imag, Ls*omega[frame]*Ia[frame].real), angles='xy', scale_units='xy', scale=1, width = 0.005, label = 'j Xs Ia')
    plt.xlim(-100, 1100)
    plt.ylim(-100, 1100)
    plt.legend()
    plt.grid()
    plt.title(r"Phase diagram for RPM = {}".format(RPM[frame]))
 
def convert_to_latex(RPM, P, Va, Ea, Ia, P_loss, Sg, eff, case):
    header = r"\begin{table}[H]"
    header += r"\centering"
    header += r"\begin{adjustbox}{width=1\textwidth}"
    header += r"\begin{tabular}{|c|c|c|c|c|c|c|c|} \hline "
    header += "RPM & P [kW] & Va [V] & Ea[V] & Ia [A] & P\\textsubscript{loss} [kW] & Sg [kW] & $\\eta$\\textsubscript{gen} [\%]  \\\\ \hline"
    rows = "\n".join(
        f"{r:.1f} & {p*10**(-3):.0f} & {vr:.1f} & {er:.1f} & {ir:.1f} & {pl*10**(-3):.1f} & {sr*10**(-3):.1f} & {eff:.1f} \\\\"
        for r, p, vr, er, ir, pl, sr, eff in zip(RPM, P, Va, Ea, Ia, P_loss, Sg, eff)
    )
    footer = r"\hline \end{tabular}"
    footer += r"\end{adjustbox}"
    footer += r"\caption{Caption}"
    footer += r"\label{tab:1"+case+"}"
    footer += r"\end{table}"
    return header + "\n" + rows + "\n" + footer

def convert_to_latex_2(RPM, V1, S1, V2, S2, Vpoc, Spoc, e_tf, e_c, e_tm, case):
    header = r"\begin{table}[H]"
    header += r"\centering"
    header += r"\begin{adjustbox}{width=1\textwidth}"
    header += r"\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|} \hline "
    header += r"RPM & V1 [V] & S1 [kW] & V2 [V] & S2 [kW] & Vpoc [V] & Spoc [kW] & $\nu$\textsubscript{tf} [\%] & $\nu$\textsubscript{c} [\%] & $\nu$\textsubscript{tm} [\%] \\ \hline"
    rows = "\n".join(
        f"{r:.1f} & {v1:.1f} & {s1*10**(-3):.1f} & {v2:.1f} & {s2*10**(-3):.1f} & {vpoc:.1f} & {spoc*10**(-3):.1f} & {etf:.1f} & {ec:.1f} & {etm:.1f} \\\\"
        for r, v1, s1, v2, s2, vpoc, spoc, etf, ec, etm  in zip(RPM, V1, S1, V2, S2, Vpoc, Spoc, e_tf, e_c, e_tm)
    )
    footer = r"\hline \end{tabular}"
    footer += r"\end{adjustbox}"
    footer += r"\caption{Caption}"
    footer += r"\label{tab:2"+case+"}"
    footer += r"\end{table}"
    return header + "\n" + rows + "\n" + footer
#%%
"""
Functions
"""
def input_power(RPM):
    P = np.zeros_like(RPM)
    for i in range(len(RPM)):
        if RPM[i] <= 5.6:
            P[i] = 0
        elif RPM[i] >= 22.5:
            P[i] = 4e6 / 3
        else:
            P[i] = 351.2 * RPM[i] ** 3 / 3
    return P

def equations_a(vars, Ea, omega, P_phase):
    Ia, Va, delta = vars
    eq1 = Ea**2 - (Va + Ra * Ia)**2 - (Ls * Ia * omega)**2
    eq2 = P_phase - Ea * np.cos(delta) * Ia
    eq3 = np.cos(delta) - (Va + Ra * Ia) / Ea
    return [eq1, eq2, eq3]
    
def equations_b(vars, Ea, omega, Ia):
    Va, delta = vars
    # Ia = P_phase / Ea
    eq1 = Va * np.cos(delta) - (Ea - Ra * Ia)
    eq2 = Va * np.sin(delta) - Ia * omega * Ls
    return [eq1, eq2]

def complex_sol(Sg_lst, Ia_lst, Ea_lst, delta_lst, Va_lst, case):
    if case == 'a':
        Ea = Ea_lst * np.cos(delta_lst) + 1j*Ea_lst * np.sin(delta_lst)
        Sg = Sg_lst + 1j * 0
        Ia = Ia_lst + 1j * 0
        Va = Va_lst + 1j * 0
        return Ea, Sg, Ia, Va

    if case == 'b':
        Ia = Ia_lst * np.cos(delta_lst) + 1j* Ia_lst * np.sin(delta_lst)
        Sg = Sg_lst * np.cos(-delta_lst) + 1j * Sg_lst * np.sin(-delta_lst)
        Ea = Ea_lst * np.cos(delta_lst) + 1j * Ea_lst * np.sin(delta_lst)
        Va = Va_lst + 1j * 0
        return Ea, Sg, Ia, Va

def bring_to_secondary(x):
    return x/RATIO**2

def bring_to_primary(x):
    return x*RATIO**2

def print_results(V1, V2, Vpoc, I1, I2, Ipoc):
    complex_numbers = {'V1': V1, 'V2': V2, 'Vpoc': Vpoc, 'I1': I1, 'I2': I2, 'Ipoc': Ipoc}
    for name, value in complex_numbers.items():
        # Format with sign handling
        formatted_value = f"{value.real:.2f} {'+' if value.imag >= 0 else '-'} {abs(value.imag):.2f}i"
        print(f"{name} = {formatted_value}")

#%%
"""
CONSTANTS
"""
case = 'b'   ####### HELLO LOOK AT ME 

Np = 52
phi_nom = 15.826 
Ra = 14.826e-3
Ls = 5.573e-3

RPM = np.arange(22.5, 0, -0.5)
V = RPM/5.6*3
f = RPM * Np / 120
omega = f * 2 * np.pi
P_phase = input_power(RPM)
P = P_phase*3

#%%
"""
SYNCHRONOUS GENERATOR   
"""
def generator(case):
    Ea_lst = omega * phi_nom  
    Va_lst = np.full_like(omega, np.nan) #Define empty lists
    Ia_lst = np.full_like(omega, np.nan)
    delta_lst = np.full_like(omega, np.nan)

    if case == 'a':
        print('Running case A')
        initial_guess = [1740.45, 740.29, 0.6597]

        for i, Ea in enumerate(Ea_lst):
            if P_phase[i] == 0:
                Va_lst[i], Ia_lst[i], delta_lst[i], Ea_lst[i] = 0, 0, 0, 0
            else:
                try:
                    # Use fsolve to find the roots
                    solution = fsolve(equations_a, initial_guess, args=(Ea, omega[i], P_phase[i]))
                    Ia_lst[i], Va_lst[i], delta_lst[i] = solution
                    initial_guess = [Ia_lst[i], Va_lst[i], 0.6597]
                except Exception as e:
                    print(f"Solution not found for omega index {i}: {e}")

    if case == 'b':
        print('Running case B')
        initial_guess = [1058.92, 0.4593] 
        
        for i, Ea in enumerate(Ea_lst):
            if P_phase[i] == 0:
                Va_lst[i], Ia_lst[i], delta_lst[i], Ea_lst[i] = 0, 0, 0, 0
            else:
                Ia_lst[i] = P_phase[i] / Ea_lst[i]
                try:
                    # Use fsolve to find the roots
                    solution = fsolve(equations_b, initial_guess, args=(Ea, omega[i], Ia_lst[i]))
                    Va_lst[i], delta_lst[i] = solution
                    initial_guess = [Va_lst[i], 0.4593]
                except Exception as e:
                    print(f"Solution not found for omega index {i}: {e}")

    return Va_lst, Ia_lst, delta_lst, Ea_lst

Va_lst, Ia_lst, delta_lst, Ea_lst = generator(case)
Sg_lst = np.zeros_like(Va_lst)
P_loss = []
eff = []

for i in range(len(Va_lst)):
    Sg_lst[i] = Va_lst[i]*Ia_lst[i]*3

Ea, Sg, Ia, Va = complex_sol(Sg_lst, Ia_lst, Ea_lst, delta_lst, Va_lst ,case)

for i in range(len(Va_lst)):
    P_loss.append(P[i]-Sg[i].real)
    eff.append(Sg[i].real/(P_phase[i]*3)*100)

plots_1()


# Creating the animation
fig, ax = plt.subplots()
ani = animation.FuncAnimation(fig, update, frames=len(Va_lst), repeat=True)

f = r"animation_{}.gif".format(case)
writergif = animation.PillowWriter(fps=7, codec='libx264', bitrate=2) 
ani.save(f, writer=writergif)

plt.show()


#%%
"""
TRANSFORMER
"""
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

#Cable
lenght = 60
r_cable = 0.5*lenght
l_cable = 0.5e-3*lenght
c_cable = 0.1e-6*lenght

######
z1 = R1 + L1*omega*1j
ze = 1 / (1/RC + 1 / (1j*omega*LM))
zcomp_p = R2_p + L2_p*omega*1j + bring_to_primary(r_cable) + bring_to_primary(l_cable*omega*1j)
zc_p = bring_to_primary(1/(c_cable*omega*1j))

Vpoc_p = Vpoc/m.sqrt(3)*RATIO

Ic_p = Vpoc_p/zc_p

## Set up system of equations
def complex_equations(vars, z1, ze, zcomp_p, zc_p, S, R2_p, omega):
    # Unpack the variables
    V1_real, V1_imag, I1_real, I1_imag, Ie_real, Ie_imag, I2_p_real, I2_p_imag, Ipoc_p_real, Ipoc_p_imag = vars
    
    # Convert to complex numbers
    V1 = V1_real + 1j*V1_imag
    I1 = I1_real + 1j*I1_imag
    Ie = Ie_real + 1j*Ie_imag
    I2_p = I2_p_real + 1j*I2_p_imag
    Ipoc_p = Ipoc_p_real + 1j*Ipoc_p_imag
    
    # Define the equations
    eq1 = V1 - I1*z1 - Ie*ze
    eq2 = Ie*ze - I2_p*zcomp_p - Ic_p * zc_p
    eq3 = I1 - Ie - I2_p
    eq4 = I2_p - Ic_p - Ipoc_p
    eq5 = I1 - np.conj(S)/np.conj(V1)
    
    # Split the equations into real and imaginary parts
    return [
        eq1.real, eq1.imag, 
        eq2.real, eq2.imag, 
        eq3.real, eq3.imag, 
        eq4.real, eq4.imag, 
        eq5.real, eq5.imag
    ]

def solve_system(S, z1, ze, zcomp_p, zc_p, R2_p, omega,initial_guess):

    # Use fsolve to find the roots
    solution = fsolve(complex_equations, initial_guess, args=(z1, ze, zcomp_p, zc_p, S, R2_p, omega))
    
    # Extract the solutions
    V1 = solution[0] + 1j*solution[1]
    I1 = solution[2] + 1j*solution[3]
    Ie = solution[4] + 1j*solution[5]
    I2_p = solution[6] + 1j*solution[7]
    Ipoc_p = solution[8] + 1j*solution[9]
    
    initial_guess=[solution[0],solution[1],solution[2],solution[3],solution[4],solution[5],solution[6],solution[7],solution[8],solution[9]]    
    # Calculate V2_p based on the solved values
    #V1 = np.conj(V1)
    V2_p = V1 - I1*z1 - I2_p*(R2_p+1j*omega*L2_p)
    
    # Rotate the solution so that V1 has an angle of 0 radians
    delta = m.atan(V1.imag/V1.real)

    V1_rotated = V1 * np.exp(-1j*delta)
    V2_p_rotated = V2_p * np.exp(-1j*delta)
    Vpoc_p_rotated = Vpoc_p * np.exp(-1j*delta)
    I1_rotated = I1 * np.exp(-1j*delta)
    I2_p_rotated = I2_p * np.exp(-1j*delta)
    Ie_rotated = Ie * np.exp(-1j*delta)
    Ipoc_p_rotated = Ipoc_p * np.exp(-1j*delta)
    
    return V1_rotated, V2_p_rotated, Vpoc_p_rotated, I1_rotated, I2_p_rotated, Ie_rotated, Ipoc_p_rotated, initial_guess

#Initialise solution lists
V1_lst = np.zeros(len(RPM),dtype=np.complex_)
V2_p_lst = np.zeros_like(V1_lst)
Vpoc_p_lst = np.zeros_like(V1_lst)

I1_lst = np.zeros_like(V1_lst)
Ie_lst = np.zeros_like(V1_lst)
I2_p_lst = np.zeros_like(V1_lst)
Ipoc_p_lst = np.zeros_like(V1_lst)
Ic_p_lst = np.zeros_like(V1_lst)

S_lst = np.zeros_like(V1_lst)

#Main loop for q2
initial_guess = [490.98, 34.04, 2641.5, -281.68, 6.3987, -58.7008, 2641.5, -1999.3, 0, 1717.6]  # 5 complex variables split into real and imaginary parts
for i in range(len(RPM)):
    S_lst[i] = complex(Sg.real[i]/3, Sg.real[i]*0.2/3)
    V1_lst[i], V2_p_lst[i],Vpoc_p_lst[i], I1_lst[i], I2_p_lst[i], Ie_lst[i], Ipoc_p_lst[i], initial_guess = solve_system(S_lst[i], z1, ze, zcomp_p, zc_p, R2_p, omega, initial_guess)

S1_lst = V1_lst * np.conj(I1_lst) * 3
S2_lst = V2_p_lst * np.conj(I2_p_lst) * 3
Spoc_lst = Vpoc_p_lst * np.conj(Ipoc_p_lst) * 3

# Call the function with the first element of each list
print_results(V1_lst[0], V2_p_lst[0], Vpoc_p_lst[0], I1_lst[0], I2_p_lst[0], Ipoc_p_lst[0])


# %%
"""
EFFICIENCY OF THE WT SYSTEM
"""
eff_transformer = S2_lst.real / Sg.real * 100
eff_cable = Spoc_lst.real / S2_lst.real * 100
eff_transmission = Spoc_lst.real / Sg.real *100
eff_system = Spoc_lst.real / (P_phase*3) * 100

print('EFFICIENCIES:')
print('\ttransformer',eff_transformer[0])
print('\tcable',eff_cable[0])
print('\ttransmission',eff_transmission[0])
print('\tsystem',eff_system[0])


latex_table = convert_to_latex(RPM, P, Va, Ea, Ia, P_loss, Sg, eff, case)
latex_table_2 = convert_to_latex_2(RPM, V1_lst, S1_lst, V2_p_lst, S2_lst, Vpoc_p_lst, Spoc_lst, eff_transformer, eff_cable, eff_transmission, case)
# Save the LaTeX table to a .txt file
with open('latex_table.txt', 'w') as file:
    file.write(latex_table)
with open('latex_table_2.txt', 'w') as file:
    file.write(latex_table_2)

plots_2()
plot_efficiency(eff_system)
plt.show()


# %%
