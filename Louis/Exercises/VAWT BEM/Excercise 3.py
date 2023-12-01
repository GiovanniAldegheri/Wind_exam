import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

os.chdir(r'G:\Other computers\Grote Laptop\Desktop\TU Delft\MSc EWEM 1\Q1-2 DTU\45300 Wind turbine technology and aerodynamics\Shared Git\Wind_exam\Louis\Exercises\VAWT BEM')


# Load Input data from blade data file
Airfoil_path = 'airfoil.txt'
Airfoil_data = np.loadtxt(Airfoil_path, skiprows=0)
AoA_lst = Airfoil_data[:, 0]
Cl_lst = Airfoil_data[:, 1]
Cd_lst = Airfoil_data[:, 2]

Exam = True

# %% Function


def VAWT_BEM(V_0, B, S, R, omega, dt, t_end):
    # Initialise t and a
    t = np.arange(0, t_end+dt, dt)
    a = np.zeros(len(t))
    print(t)
    # initialise total rotated angle and local blade angle, relative to azimuth
    theta_gen = np.zeros(len(t))
    theta = np.zeros((len(t), B))

    # Initialise load vectors
    p_x = np.zeros((len(t), B))
    p_y = np.zeros((len(t), B))
    p_y_tot = np.zeros(len(t))
    p_x_tot = np.zeros(len(t))
    p_tot = np.zeros(len(t))

    p_n = np.zeros((len(t), B))
    p_t = np.zeros((len(t), B))

    C_T = np.zeros(len(t))
    C_P = np.zeros(len(t))

    AoA_1 = np.zeros((len(t), B))

    # Solver loop
    for t_i in range(1, len(t)):
        # Update total rotated angle based on constant omega
        theta_gen[t_i] = theta_gen[t_i-1] + omega*dt

        # TODO: Calculate W_0: induced wind at x=0 plane
        W_0 = a[t_i - 1] * V_0

        # Loop over every blade to find individual blade loads
        for b_i in range(B):
            # Calculate local angular position relative to azimuth
            theta[t_i, b_i] = theta_gen[t_i] + 2*np.pi*(b_i-1)/B

            # Calc W_x,i and W_y,i
            x = -R * np.sin(theta[t_i, b_i])
            y = R * np.cos(theta[t_i, b_i])
            W_x = W_0 * (1 - 0.4 * np.sin(theta[t_i, b_i]))
            W_y = 0.4 * W_x * np.cos(theta[t_i, b_i])

            # Calculate all relative wind speeds
            V_rel_x = omega*y + V_0 - W_x
            V_rel_y = -omega*x + W_y
            V_norm = (V_0 - W_x)*np.sin(theta[t_i, b_i]) - W_y*np.cos(theta[t_i, b_i])
            V_tan = (V_0 - W_x)*np.cos(theta[t_i, b_i]) + W_y*np.sin(theta[t_i, b_i]) + omega*R
            V_rel = np.sqrt(V_norm ** 2 + V_tan ** 2)

            # Interpolate for Cd and Cl
            AoA = np.arctan(V_norm/V_tan)
            Cl = np.interp(np.degrees(AoA), AoA_lst, Cl_lst)
            Cd = np.interp(np.degrees(AoA), AoA_lst, Cd_lst)

            # Store AoA for the blade
            AoA_1[t_i, b_i] = AoA

            # Find lift and drag forces
            l = 0.5 * rho * S * R / B * Cl * V_rel**2
            d = 0.5 * rho * S * R / B * Cd * V_rel**2

            # find p_x, p_y --> cb = cos beta, sb = sin beta
            cb = V_rel_y / V_rel
            sb = V_rel_x / V_rel

            p_x[t_i, b_i] = l * cb + d * sb
            p_y[t_i, b_i] = -l * sb + d * cb

            # find p_n, p_t
            p_n[t_i, b_i] = l * np.cos(AoA) + d * np.sin(AoA)
            p_t[t_i, b_i] = l * np.sin(AoA) - d * np.cos(AoA)

        # Get total force in y, x and combined for all blades
        p_y_tot[t_i] = np.sum(p_y[t_i, :])
        p_x_tot[t_i] = np.sum(p_x[t_i, :])
        p_tot[t_i] = np.sqrt(p_y_tot[t_i]**2 + p_x_tot[t_i]**2)

        # Calculate CT and CP
        dTdZ = np.sum(p_x[t_i, :])
        dPdZ = omega*R*np.sum(p_t[t_i, :])

        C_T[t_i] = dTdZ/(rho*R*V_0**2)
        C_P[t_i] = dPdZ/(rho*R*V_0**3)

        if a[t_i-1] <= 1/3:
            f_g = 1
        if a[t_i-1] > 1/3:
            f_g = 1/4 * (5-3*a[t_i-1])

        a_qs = C_T[t_i]/(4*(1-f_g*a[t_i-1]))            # NOQA
        a[t_i] = a_qs + (a[t_i-1] - a_qs)*np.exp(-dt/tau)

    return t[1:], p_x[1:], p_y[1:], p_n[1:], p_t[1:], p_x_tot[1:], p_y_tot[1:], p_tot[1:], C_T[1:], C_P[1:], AoA_1[1:], W_0


def VAWT_EXAM_BEM(V_0, B, S, R, omega, dt):
    # Initialise t and a
    t = np.arange(0, t_end+dt, dt)
    a = np.zeros(len(t))

    # initialise total rotated angle and local blade angle, relative to azimuth
    theta_gen = np.zeros(len(t))
    theta = np.zeros((len(t), B))

    # Initialise load vectors
    p_x = np.zeros((len(t), B))
    p_y = np.zeros((len(t), B))
    p_y_tot = np.zeros(len(t))
    p_x_tot = np.zeros(len(t))
    p_tot = np.zeros(len(t))

    p_n = np.zeros((len(t), B))
    p_t = np.zeros((len(t), B))

    C_T = np.zeros(len(t))
    C_P = np.zeros(len(t))

    AoA_1 = np.zeros((len(t), B))

    # Solver loop
    for t_i in range(1, len(t)):
        # Update total rotated angle based on constant omega
        theta_gen[t_i] = theta_gen[t_i-1] + omega*dt

        # TODO: Calculate W_0: induced wind at x=0 plane
        W_0 = a[t_i - 1] * V_0

        # Loop over every blade to find individual blade loads
        for b_i in range(B):
            # Calculate local angular position relative to azimuth
            theta[t_i, b_i] = theta_gen[t_i] + 2*np.pi*(b_i-1)/B

            # Calc W_x,i and W_y,i
            x = -R * np.sin(theta[t_i, b_i])
            y = R * np.cos(theta[t_i, b_i])
            W_x = W_0 * (1 - 0.4 * np.sin(theta[t_i, b_i]))
            W_y = 0.4 * W_x * np.cos(theta[t_i, b_i])
            # W_x = 3.4
            # W_y = 0

            # Calculate all relative wind speeds
            V_rel_x = omega*y + V_0 - W_x
            V_rel_y = -omega*x + W_y
            V_norm = (V_0 - W_x)*np.sin(theta[t_i, b_i]) - W_y*np.cos(theta[t_i, b_i])
            V_tan = (V_0 - W_x)*np.cos(theta[t_i, b_i]) + W_y*np.sin(theta[t_i, b_i]) + omega*R
            V_rel = np.sqrt(V_norm ** 2 + V_tan ** 2)

            # Interpolate for Cd and Cl
            AoA = np.arctan(V_norm/V_tan)
            # Cl = np.interp(np.degrees(AoA), AoA_lst, Cl_lst)
            # Cd = np.interp(np.degrees(AoA), AoA_lst, Cd_lst)
            Cl = 0.1 * AoA
            Cd = 0.01

            # Store AoA for the blade
            AoA_1[t_i, b_i] = AoA

            # Find lift and drag forces
            l = 0.5 * rho * S * R / B * Cl * V_rel**2
            d = 0.5 * rho * S * R / B * Cd * V_rel**2
            if b_i == 1:
                print(f'{np.degrees(theta[t_i, b_i]):.2f}, {b_i}, lift = {l:.2f}')

            # find p_x, p_y --> cb = cos beta, sb = sin beta
            cb = V_rel_y / V_rel
            sb = V_rel_x / V_rel

            p_x[t_i, b_i] = l * cb + d * sb
            p_y[t_i, b_i] = -l * sb + d * cb

            # find p_n, p_t
            p_n[t_i, b_i] = l * np.cos(AoA) + d * np.sin(AoA)
            p_t[t_i, b_i] = l * np.sin(AoA) - d * np.cos(AoA)

        # Get total force in y, x and combined for all blades
        p_y_tot[t_i] = np.sum(p_y[t_i, :])
        p_x_tot[t_i] = np.sum(p_x[t_i, :])
        p_tot[t_i] = np.sqrt(p_y_tot[t_i]**2 + p_x_tot[t_i]**2)

        # Calculate CT and CP
        dTdZ = np.sum(p_x[t_i, :])
        dPdZ = omega*R*np.sum(p_t[t_i, :])

        C_T[t_i] = dTdZ/(rho*R*V_0**2)
        C_P[t_i] = dPdZ/(rho*R*V_0**3)

        if a[t_i-1] <= 1/3:
            f_g = 1
        if a[t_i-1] > 1/3:
            f_g = 1/4 * (5-3*a[t_i-1])

        a_qs = C_T[t_i]/(4*(1-f_g*a[t_i-1]))            # NOQA
        a[t_i] = a_qs + (a[t_i-1] - a_qs)*np.exp(-dt/tau)

    return t[1:], p_x[1:], p_y[1:], p_n[1:], p_t[1:], p_x_tot[1:], p_y_tot[1:], p_tot[1:], C_T[1:], C_P[1:], AoA_1[1:], W_x


# %% Exam Setup
if Exam:

    # Set up global parameters
    t_end = 7                           # End time [s]
    B = int(2)                          # Number of blades [-]
    R = 3                               # Blade radius [m]
    V_0 = 9                             # Free stream wind speed [m/s]
    omega = 115                         # Rotational speed [RPM]
    omega = np.radians(omega)*6         # Rotation speed [rad/s]
    omega = 14                          # rad/s

    # # Constant shape factor
    # S = 0.2                             # Factor: S = Bc/R
    # c = S * R / B                       # Chord length [m], dependent on #blades, shape factor and total radius

    # # Constant cord
    c = 0.3                             # m
    S = B*c/R                           # [-]
    rho = 1.225                         # Density [kg/m3]
    tau = 2*R/V_0                       # Time factor [-]
    plot_start = 6                      # Indicate start time for plots
    dt = 5e-4

    t, p_x, p_y, p_n, p_t, p_x_tot, p_y_tot, p_tot, C_T, C_P, AoA_1, W_x = VAWT_EXAM_BEM(V_0, B, S, R, omega, dt)

    # # Calc for azimuthal positions
    angles = np.array([0, np.pi/2, np.pi, np.pi*3/2])
    W_x = V_0 * (1 - 0.4 * np.sin(angles))
    for i in range(len(angles)):
        print(f'Angle = {np.degrees(angles[i]):<5} [deg], Wx = {W_x[i]:.5} [m/s]')

    # plt.plot(np.degrees(t*omega), )
    # # Check AoA
    plt.plot(np.degrees(t*omega), AoA_1)
    plt.xlim([0, 360])
    plt.grid()
    plt.show()

# %% Exercise Setup
if not Exam:
    # Set up global parameters
    t_end = 7               # End time [s]
    B = np.array([int(3)])              # Number of blades [-]
    R = 3                   # Blade radius [m]
    V_0 = 8                 # Free stream wind speed [m/s]
    omega = 14              # Rotation speed [rad/s]
    S = 0.2                 # Factor: S = Bc/R
    c = S * R / B           # Chord length [m], dependent on #blades, shape factor and total radius
    rho = 1.225             # Density [kg/m3]
    tau = 2*R/V_0           # Time factor [-]
    plot_start = 0          # Indicate start time for plots
    dt = 1e-3

    B = np.arange(3, 4)
    dt = 1e-3
    fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=2, figsize=(14, 8))
    labels = []
    for b_i in B:
        print(b_i)
        t, p_x, p_y, p_n, p_t, p_x_tot, p_y_tot, p_tot, C_T, C_P, AoA_1, W_0 = VAWT_BEM(V_0, b_i, S, R, omega, dt, t_end)

        # Calculate moving averages
        C_T = pd.Series(C_T)
        C_P = pd.Series(C_P)
        moving_CT = (C_T.ewm(alpha=0.0003, adjust=True).mean()).tolist()
        moving_CP = (C_P.ewm(alpha=0.0003, adjust=True).mean()).tolist()

        # Plot upper left - CT
        l0 = ax0[0].plot(t[int(plot_start/dt):], C_T[int(plot_start/dt):])
        l1 = ax0[0].plot(t[int(plot_start/dt):], moving_CT[int(plot_start/dt):])

        # Plot lower left - CP
        l2 = ax1[0].plot(t[int(plot_start/dt):], C_P[int(plot_start/dt):])
        l2 = ax1[0].plot(t[int(plot_start / dt):], moving_CP[int(plot_start / dt):])

        # Plot upper right - Choose (Currently p_t [for all blades])
        l3 = ax0[1].plot(t[int(plot_start / dt):], p_t[int(plot_start / dt):])
        # l3 = ax1[1].plot(t[int(plot_start/dt):], p_tot[int(plot_start/dt):])
        # l3 = ax0[1].plot(t[int(plot_start / dt):], AoA_1[int(plot_start / dt):])
        # l3 = ax0[1].plot(t[int(plot_start / dt):], p_n[int(plot_start / dt):])

        # Plot lower right - Choose (Currently p_x_tot and p_y_tot)
        l4 = ax1[1].plot(t[int(plot_start / dt):], p_x_tot[int(plot_start / dt):])
        l4 = ax1[1].plot(t[int(plot_start / dt):], p_y_tot[int(plot_start / dt):])
        labels.append(f'{b_i} blade(s)')

    # all x-labels
    ax0[0].set_xlabel('t [s]')
    ax0[1].set_xlabel('t [s]')
    ax1[0].set_xlabel('t [s]')
    ax1[1].set_xlabel('t [s]')

    # all y-labels
    ax0[0].set_ylabel(r'$C_T$ [-]')
    ax0[1].set_ylabel(r'$p_t$ and $p_n$ [-]')
    ax1[0].set_ylabel(r'$C_P$ [-]')
    ax1[1].set_ylabel(r'$p_{tot}$ [N]')

    # random y_limit
    ax0[0].set_ylim([min(C_T), max(C_T)])   # NOQA
    # ax0[1].set_ylim([min(C_T), max(C_T)])
    fig.legend(labels=labels)
    plt.show()
