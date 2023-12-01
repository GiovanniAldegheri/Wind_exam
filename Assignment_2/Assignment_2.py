import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import math as m
from sympy.solvers import solve
import sympy as sym
#%%
"""
PLOTS
"""

def plots_1():
    fig, axs = plt.subplots(2, 2, figsize=(15, 10)) # 2 rows, 2 columns

    # Plot 1: Current, Terminal-, Induced Voltage
    axs[0, 0].set_title('Current, Terminal-, Induced Voltage')
    axs[0, 0].grid(True, linestyle='--', alpha=0.5)
    axs[0, 0].plot(RPM, Ia_lst, marker='o', markersize = 3, linestyle='--', label='Phase Current conjugate (A)')
    axs[0, 0].plot(RPM, Va_lst, marker='s', markersize = 3, linestyle='--', label='Terminal Voltage Conjugate (V)')
    axs[0, 0].plot(RPM, Ea_lst, marker='^', markersize = 3, linestyle='--', label='Induced Voltage Conjugate (V)')
    axs[0, 0].set_xlim(5, 25)
    axs[0, 0].set_xlabel('RPM')
    axs[0, 0].set_ylabel('Voltage/Current')
    axs[0, 0].legend()

    # Plot 2: Power
    axs[0, 1].set_title('Power')
    axs[0, 1].grid(True, linestyle='--', alpha=0.5)
    axs[0, 1].plot(RPM, Sg.imag, marker='o', markersize = 3, linestyle='--', label='Reactive Power')
    axs[0, 1].plot(RPM, Sg.real, marker='s', markersize = 3, linestyle='--', label='Active Power')
    axs[0, 1].plot(RPM, Sg_lst, marker='^', markersize = 3, linestyle='--', label='Mechanical Power')
    axs[0, 1].set_xlim(5, 25)
    axs[0, 1].set_xlabel('RPM')
    axs[0, 1].set_ylabel('Power (W)')
    axs[0, 1].legend()

    # Plot 3: Power Loss
    axs[1, 0].set_title('Power Loss')
    axs[1, 0].plot(RPM, P_loss, marker='o', markersize = 3, linestyle='--')
    axs[1, 0].set_xlim(5, 25)
    axs[1, 0].set_xlabel('RPM')
    axs[1, 0].set_ylabel('Power loss (W)')
    axs[1, 0].grid(True, linestyle='--', alpha=0.5)

    # Plot 4: Efficiency
    axs[1, 1].set_title('Efficiency')
    axs[1, 1].plot(RPM, eff, marker='o', markersize = 3, linestyle='--')
    axs[1, 1].set_xlim(5, 25)
    axs[1, 1].set_xlabel('RPM')
    axs[1, 1].set_ylabel('Efficiency (%)')
    axs[1, 1].grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()

def plots_2():
    fig, axs = plt.subplots(2, 2, figsize=(15, 10)) # 2 rows, 2 columns

    # Plot 1: Voltage
    axs[0, 0].set_title('Voltage')
    axs[0, 0].grid(True, linestyle='--', alpha=0.5)
    axs[0, 0].plot(RPM, abs(V1_lst), marker='o', markersize = 3, linestyle='--', label='Transformer Primary Voltage')
    axs[0, 0].plot(RPM, abs(V2_p_lst), marker='s', markersize = 3, linestyle='--', label='Transformer Secondary Voltage')
    axs[0, 0].plot(RPM, abs(Vpoc_lst), marker='^', markersize = 3, linestyle='--', label='POC Volotage')
    axs[0, 0].set_xlim(5, 25)
    axs[0, 0].set_xlabel('RPM')
    axs[0, 0].set_ylabel('Voltage')
    axs[0, 0].legend()

    # Plot 2: Power
    axs[0, 1].set_title('Power')
    axs[0, 1].grid(True, linestyle='--', alpha=0.5)
    axs[0, 1].plot(RPM, Sg.real, marker='o', markersize = 3, linestyle='--', label='Transformer Primary Active Power')
    axs[0, 1].plot(RPM, S2_lst.real, marker='s', markersize = 3, linestyle='--', label='Transformer Secondary Active Power')
    axs[0, 1].plot(RPM, Spoc_lst.real, marker='^', markersize = 3, linestyle='--', label='POC Active Power')
    axs[0, 1].set_xlim(5, 25)
    axs[0, 1].set_xlabel('RPM')
    axs[0, 1].set_ylabel('Power (W)')
    axs[0, 1].legend()

    # Plot 3: Power Loss
    axs[1, 0].set_title('Power')
    axs[1, 0].grid(True, linestyle='--', alpha=0.5)
    axs[1, 0].plot(RPM, Sg.real*0.2, marker='o', markersize = 3, linestyle='--', label='Transformer Primary Reactive Power')
    axs[1, 0].plot(RPM, S2_lst.imag, marker='s', markersize = 3, linestyle='--', label='Transformer Secondary Reactive Power')
    axs[1, 0].plot(RPM, Spoc_lst.imag, marker='^', markersize = 3, linestyle='--', label='POC Reactive Power')
    axs[1, 0].set_xlim(5, 25)
    axs[1, 0].set_xlabel('RPM')
    axs[1, 0].set_ylabel('Power (W)')
    axs[1, 0].legend()

    # Plot 4: Efficiency
    axs[1, 1].set_title('Efficiency')
    # axs[1, 1].plot(RPM, eff, marker='x', linestyle='--')
    axs[1, 1].set_xlabel('RPM')
    axs[1, 1].set_ylabel('Efficiency (%)')
    axs[1, 1].grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
 
def convert_to_latex(RPM, P, Va, Ea, Ia, P_loss, Sg,case):
    header = r"\begin{table}[H]"
    header += r"\centering"
    header += r"\begin{adjustbox}{width=1\textwidth}"
    header += r"\begin{tabular}{|c|c|c|c|c|c|c|} \hline "
    header += "RPM & P [kW] & Va [V] & Ea[V] & Ia [A] & P\\textsubscript{loss} [kW] & Sg [kW] \\\\ \hline"
    rows = "\n".join(
        f"{r:.1f} & {p*10**(-3):.0f} & {vr:.1f} & {er:.1f} & {ir:.1f} & {pl*10**(-3):.1f} & {sr*10**(-3):.1f} \\\\"
        for r, p, vr, er, ir, pl, sr in zip(RPM, P, Va, Ea, Ia, P_loss, Sg)
    )
    footer = r"\hline \end{tabular}"
    footer += r"\end{adjustbox}"
    footer += r"\caption{Caption}"
    footer += r"\label{tab:1"+case+"}"
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

def rotate(var, delta):
    var = complex(var)
    x = var.real
    y = var.imag
    if x == 0:
        angle = m.pi*m.sin(y)/2
    else:
        angle = m.atan(y/x)
    module = m.sqrt(x**2+y**2)
    x = module*m.cos(angle+delta)
    y = module*m.sin(angle+delta)
    return complex(x,y)

def comp_multiply(V,I):
    S = 3*(abs(V)*abs(I)*np.cos(np.tanh(V.imag/V.real)-(I.imag/I.real))+1j*abs(V)*abs(I)*np.sin(np.tanh(V.imag/V.real)-(I.imag/I.real)))
    # S_re = [S[i].real for i in range(len(S))]
    # S_im = [S[i].imag for i in range(len(S))]
    return S

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

latex_table = convert_to_latex(RPM, P, Va, Ea, Ia, P_loss, Sg,case)
# Save the LaTeX table to a .txt file
with open('latex_table.txt', 'w') as file:
    file.write(latex_table)

plots_1()
# plt.show()
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
def transformer(S):
    if S == 0:
        return 0, 0, 0, 0, 0, 0, 0, 0
    else:
        V1, I1, Ie, I2_p, Ipoc_p = sym.symbols("V1, I1, Ie, I2_p, Ipoc_p")

        eq1 = sym.Eq(V1 - I1*z1 - Ie*ze,0)
        eq2 = sym.Eq(Ie*ze - I2_p*zcomp_p - Ic_p * zc_p,0)
        eq3 = sym.Eq(I1 - Ie - I2_p,0)
        eq4 = sym.Eq(I2_p - Ic_p - Ipoc_p,0)
        eq5 = sym.Eq(np.conj(S)-I1*V1,0)


        ## Solve system of equations
        solution = solve([eq1,eq2,eq3,eq4,eq5],(V1, I1, Ie, I2_p, Ipoc_p))
        V1, I1, Ie, I2_p, Ipoc_p = (solution[1])
        V1 = np.conj(complex(V1))
        V2_p = V1 - I1*z1 - I2_p*(R2_p+1j*omega*L2_p)

        ## Transfer everything so that V1 = 0 rad
        delta = m.atan(V1.imag/V1.real)
        return rotate(V1,-delta), rotate(V2_p,-delta), rotate(I1,-delta), rotate(I2_p,-delta), rotate(Ie,-delta), rotate(Ipoc_p,-delta), rotate(Vpoc_p,-delta), rotate(Ic_p,-delta)

#Initialise solution lists
V1_lst = np.zeros(len(RPM),dtype=np.complex_)
V2_p_lst = np.zeros_like(V1_lst)
Vpoc_lst = np.zeros_like(V1_lst)

I1_lst = np.zeros_like(V1_lst)
Ie_lst = np.zeros_like(V1_lst)
I2_p_lst = np.zeros_like(V1_lst)
Ipoc_p_lst = np.zeros_like(V1_lst)
Ic_p_lst = np.zeros_like(V1_lst)

S_lst = np.zeros_like(V1_lst)

#Main loop for q2
for i in range(len(RPM)):
    S_lst[i] = complex(Sg.real[i]/3, Sg.real[i]*0.2/3)
    V1_lst[i], V2_p_lst[i], I1_lst[i], I2_p_lst[i], Ie_lst[i], Ipoc_p_lst[i], Vpoc_lst[i], Ic_p_lst[i] = transformer(S_lst[i])
    print(len(RPM)-i)

# Ipoc_p_lst[0] = 2523700 + 1j * -2183580
S2_lst = V2_p_lst * np.conj(I2_p_lst) * 3
Spoc_lst = Vpoc_lst * np.conj(Ipoc_p_lst) * 3

print(S2_lst[0])
print(Spoc_lst[0])

print(V1_lst[0], V2_p_lst[0], I1_lst[0], I2_p_lst[0], Ipoc_p_lst[0])

plots_2()

# %%
plt.show()