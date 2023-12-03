#%%
import numpy as np
import matplotlib.pyplot as plt
import time
import datetime


def propeller_solver(R, v_0, r_lst, twist_lst, c_lst, omega, theta, B, error=1e-5, itermax=750):
    rho = 1.225
    # initialise BEM Loop
    pn_lst, pt_lst = np.zeros(len(r_lst)), np.zeros(len(r_lst))
    lmbda = omega * R / v_0

    for i in range(len(r_lst) - 1):
        # Initialise a, a_prime and error for looping
        a = 0
        a_prime = 0
        epsilon = 1
        allowed_epsilon = error
        iterno = 0

        while epsilon >= allowed_epsilon:
            # Update iteration
            iterno += 1
            if iterno > itermax:
                print(
                    f'No Convergence reached within {iterno - 1} iterations for lambda = {round(lmbda, 3)}, '
                    f'Theta = {theta}')
                break
            # Compute local flow angle
            phi = np.arctan(v_0 * (1 + a) / (omega * r_lst[i] * (1 - a_prime)))
            # Compute effective AoA
            alpha = -np.degrees(phi) + theta + twist_lst[i]
            # Obtain aerodynamic coefficients
            Cl = 0.1*alpha+0.4
            Cd = 0.008
            Cn = Cl * np.cos(phi) - Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) + Cd * np.cos(phi)

            # Thrust Coefficient
            sigma = c_lst[i] * B / (2 * np.pi * r_lst[i])
            # CT = Cn * sigma * (1 - a) ** 2 / (np.sin(phi)) ** 2

            F = 2 / np.pi * np.arccos(
                np.exp(-B / 2 * (R - r_lst[i]) / (r_lst[i] * np.sin(np.abs(phi)))))

            a_star = (1 + a)*sigma*Cn/(4*F*(np.sin(phi))**2)
            a_new = 0.1 * a_star + (1 - 0.1) * a

            a_prime_star = (1-a_prime) * sigma*Ct/(4*F*np.sin(phi)*np.cos(phi))
            a_prime_new = 0.1*a_prime_star + (1 - 0.1)*a_prime

            # Check for convergence and update a, a_prime
            epsilon1 = np.abs(a_new - a)
            epsilon2 = np.abs(a_prime_new - a_prime)
            epsilon = max(epsilon1, epsilon2)
            a, a_prime = a_new, a_prime_new

        # Compute sectional forces
        # print(phi)
        V_r = v_0 * (1 + a) / np.sin(phi)
        pn_lst[i] = 0.5 * rho * c_lst[i] * Cn * V_r ** 2
        pt_lst[i] = 0.5 * rho * c_lst[i] * Ct * V_r ** 2

    # Compute global forces and characteristics
    T = B * np.trapz(pn_lst, r_lst[:])
    M_R = B * np.trapz(r_lst[:] * pt_lst, r_lst[:])
    P = omega * M_R
    CP = P / (0.5 * rho * np.pi * (R ** 2) * v_0 ** 3)
    CT = T / (0.5 * rho * np.pi * (R ** 2) * v_0 ** 2)

    return CP, P, T, M_R, np.array(pn_lst), np.array(pt_lst), omega, CT


# Propeller input data
r_lst = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.645, 0.690])
twist_lst = np.array([23, 20, 16, 14.5, 13, 12.3, 11.7])
c_lst = np.array([0.106, 0.117, 0.112, 0.103, 0.088, 0.082, 0])

R = 0.69                        # m
omega = 2600                    # RPM
omega = 2*np.pi*omega/60        # rad/s
Power_shaft = 22e3              # W

v_b = np.arange(1, 41, 1)
# v_b = np.linspace(1, 40, 40)
B = 2
theta = 0

eta, P, T, CT, M_R, D = np.zeros(len(v_b)), np.zeros(len(v_b)), np.zeros(len(v_b)), np.zeros(len(v_b)), np.zeros(len(v_b)), np.zeros(len(v_b))
diff = 1000
v_end, idx_end = 0, 0
for i in range(len(v_b)):
    C_P, P[i], T[i], M_R[i], pn_lst, pt_lst, omega_, CT[i] = propeller_solver(R, v_b[i], r_lst, twist_lst, c_lst, omega, theta, B)
    eta[i] = T[i]*v_b[i]/Power_shaft
    D[i] = 2.44 * v_b[i] ** 2
    # diff = np.abs(T[i]-D)
    if np.abs(T[i]-D[i]) < diff:
        v_end = v_b[i]
        diff = np.abs(T[i]-D[i])
        idx_end = i

print(v_end)
plt.figure(figsize=(10, 8))
plt.subplot(221)
plt.plot(v_b, P)
plt.axhline(22000, color='r', linestyle='--')
plt.xlabel(r'$v_0$ [m/s]')
plt.ylabel(r'$P$ [W]')
plt.grid()

plt.subplot(222)
plt.plot(v_b, eta)
plt.xlabel(r'$v_0$ [m/s]')
plt.ylabel(r'$\eta$ [-]')
plt.grid()

plt.subplot(223)
plt.plot(v_b, T, label='Thrust')
plt.plot(v_b, D, label='Drag')
plt.plot(v_end, T[idx_end], marker='s', color='r', markersize=5)
plt.xlabel(r'$v_0$ [m/s]')
plt.ylabel(r'Forces [N]')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
