import numpy as np

print(3*6235 -(6235+3*3250))

files = ['Data\FFA-W3-241.txt', 'Data\FFA-W3-301.txt', 'Data\FFA-W3-360.txt', 'Data\FFA-W3-480.txt',
             'Data\FFA-W3-600.txt', 'Data\cylinder.txt']

# TEST OF INTERPOLATION ROUTINE. COMPARE TO INTERP1 IN MATLAB
def force_coeffs_10MW(angle_of_attack, thick, aoa, cl_tab, cd_tab, cm_tab, printing=False):
    # Creating a function which takes the angle of attack and the section thickness:

    thick_prof = np.zeros(6)
    # NOTE THAT IN PYTHON THE INTERPOLATION REQUIRES THAT THE VALUES INCREASE IN THE VECTOR!
    thick_prof[0] = 24.1
    thick_prof[1] = 30.1
    thick_prof[2] = 36.0
    thick_prof[3] = 48.0
    thick_prof[4] = 60.0
    thick_prof[5] = 100.0

    cl_aoa = np.zeros([1, 6])
    cd_aoa = np.zeros([1, 6])
    cm_aoa = np.zeros([1, 6])

    # Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0, i] = np.interp(angle_of_attack, aoa, cl_tab[:, i])
        cd_aoa[0, i] = np.interp(angle_of_attack, aoa, cd_tab[:, i])
        cm_aoa[0, i] = np.interp(angle_of_attack, aoa, cm_tab[:, i])

    # Interpolate to current thickness:
    cl = np.interp(thick, thick_prof, cl_aoa[0, :])
    cd = np.interp(thick, thick_prof, cd_aoa[0, :])
    cm = np.interp(thick, thick_prof, cm_aoa[0, :])

    if printing:
        print('cl:', cl)
        print('cd:', cd)
        print('cm:', cm)

    return cl, cd, cm
