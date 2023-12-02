#%%
import numpy as np
import os
import matplotlib.pyplot as plt

# os.chdir(r'Louis/Exercises')
os.chdir(r'G:/Other computers/Grote Laptop/Desktop/TU Delft/MSc EWEM 1/Q1-2 DTU/45300 Wind turbine technology and aerodynamics/Shared Git/Wind_exam/Louis/Exercises/Deflection/data')
# path = 'Louis/Exercises/Deflection/data/'

bladestruc = np.loadtxt('bladestruc.txt')

# loads_files = ['loads6.txt', 'loads11.txt', 'loads20.txt']
loads_files = ['loads_custom.txt']
pitch_angles = [np.deg2rad(0.896), 0, np.deg2rad(17.35)]


def deflection(loads, structure, pitch):
    Ty = np.zeros(len(loads[:, 0]))
    Tz = np.zeros(len(loads[:, 0]))
    My = np.zeros(len(loads[:, 0]))
    Mz = np.zeros(len(loads[:, 0]))
    kappay = np.zeros(len(loads[:, 0]))
    kappaz = np.zeros(len(loads[:, 0]))
    angley = np.zeros(len(loads[:, 0]))
    anglez = np.zeros(len(loads[:, 0]))
    deflectiony = np.zeros(len(loads[:, 0]))
    deflectionz = np.zeros(len(loads[:, 0]))

    for i, r in reversed(list(enumerate(loads[:, 0]))):
        if i > 0:
            Ty[i - 1] = Ty[i] + 0.5 * (loads[i, 2] + loads[i - 1, 2]) * (loads[i, 0] - loads[i - 1, 0])
            Tz[i - 1] = Tz[i] + 0.5 * (loads[i, 1] + loads[i - 1, 1]) * (loads[i, 0] - loads[i - 1, 0])
            My[i - 1] = My[i] - Tz[i] * (loads[i, 0] - loads[i - 1, 0]) - (
                        1 / 6 * loads[i - 1, 1] + 1 / 3 * loads[i, 1]) * (loads[i, 0] - loads[i - 1, 0]) ** 2
            Mz[i - 1] = Mz[i] + Ty[i] * (loads[i, 0] - loads[i - 1, 0]) + (
                        1 / 6 * loads[i - 1, 2] + 1 / 3 * loads[i, 2]) * (loads[i, 0] - loads[i - 1, 0]) ** 2

    for i, r in enumerate(loads[:, 0]):
        beta = np.deg2rad(structure[i, -1])
        v = np.deg2rad(structure[i, 1])
        M1 = My[i] * np.cos(beta + v + pitch) - Mz[i] * np.sin(beta + v + pitch)
        M2 = My[i] * np.sin(beta + v + pitch) + Mz[i] * np.cos(beta + v + pitch)
        k1 = M1 / (structure[i, 3])
        k2 = M2 / (structure[i, 4])

        kappay[i] = (k1 * np.cos(beta + v + pitch) + k2 * np.sin(beta + v + pitch))
        kappaz[i] = (-k1 * np.sin(beta + v + pitch) + k2 * np.cos(beta + v + pitch))

    for i, r in enumerate(loads[:-1, 0]):
        angley[i + 1] = angley[i] + 0.5 * (kappay[i + 1] + kappay[i]) * (loads[i + 1, 0] - loads[i, 0])
        anglez[i + 1] = anglez[i] + 0.5 * (kappaz[i + 1] + kappaz[i]) * (loads[i + 1, 0] - loads[i, 0])

        deflectiony[i + 1] = deflectiony[i] + anglez[i] * (loads[i + 1, 0] - loads[i, 0]) + (
                    1 / 6 * kappaz[i + 1] + 1 / 3 * kappaz[i]) * ((loads[i + 1, 0] - loads[i, 0]) ** 2)
        deflectionz[i + 1] = deflectionz[i] - angley[i] * (loads[i + 1, 0] - loads[i, 0]) - (
                    1 / 6 * kappay[i + 1] + 1 / 3 * kappay[i]) * ((loads[i + 1, 0] - loads[i, 0]) ** 2)

    return Ty, Tz, My, Mz, kappay, kappaz, angley, anglez, deflectiony, deflectionz


def nat_freq(r, structure, pitch):
    m = 0.5 * (bladestruc[1:, 2] + bladestruc[:-1, 2])
    M = np.diag(np.concatenate([m, m]))

    F = np.zeros(np.shape(M))
    N = len(r)
    # for y
    for i in range(0, N - 1):
        p = np.zeros((N, 3))
        p[:, 0] = r
        p[i + 1, 1] = 0
        p[i + 1, 2] = 1
        uy, uz = deflection(p, structure, pitch)[-2:]
        F[0:N - 1, i] = uy[1:]
        F[N - 1:, i] = uz[1:]
    # for z:
    for i in range(0, N - 1):
        p = np.zeros((N, 3))
        p[:, 0] = r
        p[i + 1, 1] = 1
        p[i + 1, 2] = 0
        uy, uz = deflection(p, structure, pitch)[-2:]
        F[:N - 1, N + i - 1] = uy[1:]
        F[N - 1:, N + i - 1] = uz[1:]

    D, V = np.linalg.eig(F @ M)
    eig_freq = 1 / np.sqrt(D[:3])
    mode_shapes = np.zeros((N, 6))

    mode_shapes[1:, 0] = V[:N - 1, 0] / np.max(abs(V[:, 0]))
    mode_shapes[1:, 1] = V[:N - 1, 1] / np.max(abs(V[:, 1]))
    mode_shapes[1:, 2] = V[:N - 1, 2] / np.max(abs(V[:, 2]))

    mode_shapes[1:, 3] = V[N - 1:, 0] / np.max(abs(V[:, 1]))
    mode_shapes[1:, 4] = V[N - 1:, 1] / np.max(abs(V[:, 2]))
    mode_shapes[1:, 5] = V[N - 1:, 2] / np.max(abs(V[:, 3]))

    return eig_freq, mode_shapes

Exercise = True
Iterative = False

if Iterative == True:
    #In this case (Exam 2020), iterate for constant pn until z deflection at tip = 5m
    r = bladestruc[:,0]
    pn_const = 0
    deflectionz = np.zeros(len(r))

    pitch_angle = 0

    while deflectionz[-1] <= 5:
        pn_const += 1
        pn = np.full(len(r),pn_const)
        pt = np.zeros(len(r))
        loads = np.transpose([r,pn,pt])

        Ty, Tz, My, Mz, kappay, kappaz, angley, anglez, deflectiony, deflectionz = deflection(loads, bladestruc, pitch_angle)

        print('pn =',pn_const,'[N/m] \t Tip deflection =',deflectionz[-1],'[m]')

    bladestruc[:, 2] = bladestruc [:, 2] * 1.1  #Mass distribution increased by 10%

    eig_freq, mode_shapes = nat_freq(r, bladestruc, pitch_angle)
    print(eig_freq)


if Exercise == True:
    for loads_file, pitch_angle in zip(loads_files, pitch_angles):
        loads = np.loadtxt(loads_file)
        loads[:,1:] = 1000*loads[:,1:]
        Ty, Tz, My, Mz, kappay, kappaz, angley, anglez, deflectiony, deflectionz = deflection(loads, bladestruc,
                                                                                            pitch_angle)
        eig_freq, mode_shapes = nat_freq(loads[:, 0], bladestruc, pitch_angle)
        
        index = np.where(loads[:,0] == 58.5344)
        print(loads_file, 'My =',My[index],'\t Mz =',Mz[index])

        plot = False
        if plot == True:
            plt.rcParams['axes.grid'] = True

            plt.figure(figsize=(16, 5))
            plt.subplot(1, 3, 1)
            plt.plot(loads[:, 0], loads[:, 2]/1000, label='py')
            plt.plot(loads[:, 0], loads[:, 1]/1000, label='pz')
            plt.title(f'Loads for {loads_file} (Pitch Angle: {np.rad2deg(pitch_angle):.2f} degrees)')
            plt.xlabel('Radius (m)')
            plt.ylabel('Load (kN)')
            plt.legend()

            plt.subplot(1, 3, 2)
            plt.plot(loads[:, 0], Ty/1000, label='Ty')
            plt.plot(loads[:, 0], Tz/1000, label='Tz')
            plt.title(f'Shear stress for {loads_file} /n (Pitch Angle: {np.rad2deg(pitch_angle):.2f} degrees)')
            plt.xlabel('Radius (m)')
            plt.ylabel('Shear stress (kN)')
            plt.legend()

            plt.subplot(1, 3, 3)
            plt.plot(loads[:, 0]/1000, My, label='My')
            plt.plot(loads[:, 0]/1000, Mz, label='Mz')
            plt.title(f'Bending Moments for {loads_file} /n (Pitch Angle: {np.rad2deg(pitch_angle):.2f} degrees)')
            plt.xlabel('Radius (m)')
            plt.ylabel('Bending Moment (kNm)')
            plt.legend()

            # Plot deflection and angles
            plt.figure(figsize=(16, 5))
            plt.subplot(1, 3, 1)
            plt.plot(loads[:, 0], np.rad2deg(kappay), label='Kappa Y')
            plt.plot(loads[:, 0], np.rad2deg(kappaz), label='Kappa Z')
            plt.title(f'Bending stiffness for {loads_file} /n (Pitch Angle: {np.rad2deg(pitch_angle):.2f} degrees)')
            plt.xlabel('Radius (m)')
            plt.ylabel('Bending Stiffness (deg/m)')
            plt.legend()

            plt.subplot(1, 3, 2)
            plt.plot(loads[:, 0], np.rad2deg(angley), label='Angle Y')
            plt.plot(loads[:, 0], np.rad2deg(anglez), label='Angle Z')
            plt.title(f'Angles for {loads_file} /n (Pitch Angle: {np.rad2deg(pitch_angle):.2f} degrees)')
            plt.xlabel('Radius (m)')
            plt.ylabel('Angle (deg)')
            plt.legend()

            plt.subplot(1, 3, 3)
            plt.plot(loads[:, 0], deflectiony, label='Deflection Y')
            plt.plot(loads[:, 0], deflectionz, label='Deflection Z')
            plt.title(f'Deflection for {loads_file} /n (Pitch Angle: {np.rad2deg(pitch_angle):.2f} degrees)')
            plt.xlabel('Radius (m)')
            plt.ylabel('Deflection (m)')
            plt.legend()
            plt.tight_layout()

            plt.figure(figsize=(15, 5))

            plt.subplot(1, 3, 1)
            plt.plot(loads[:, 0], mode_shapes[:, 0], label='Tangential Mode Shape')
            plt.plot(loads[:, 0], mode_shapes[:, 3], label='Normal Mode Shape')
            plt.title(f'Mode Shape 1 for {loads_file} /n Eigenfrequency: {eig_freq[0]:.2f} Rad/s')
            plt.xlabel('Radius (m)')
            plt.ylabel('Normalized Displacement (-)')
            plt.legend()

            plt.subplot(1, 3, 2)
            plt.plot(loads[:, 0], mode_shapes[:, 1], label='Tangential Mode Shape')
            plt.plot(loads[:, 0], mode_shapes[:, 4], label='Normal Mode Shape')
            plt.title(f'Mode Shape 2 for {loads_file} /n Eigenfrequency: {eig_freq[1]:.2f} Rad/s')
            plt.xlabel('Radius (m)')
            plt.ylabel('Normalized Displacement (-)')
            plt.legend()

            plt.subplot(1, 3, 3)
            plt.plot(loads[:, 0], mode_shapes[:, 2], label='Tangential Mode Shape')
            plt.plot(loads[:, 0], mode_shapes[:, 5], label='Normal Mode Shape')
            plt.title(f'Mode Shape 3 for {loads_file} /n Eigenfrequency: {eig_freq[2]:.2f} Rad/s')
            plt.xlabel('Radius (m)')
            plt.ylabel('Normalized Displacement (-)')
            plt.legend()

            plt.tight_layout()
        plt.show()