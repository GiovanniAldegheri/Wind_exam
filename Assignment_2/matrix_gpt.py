import numpy as np

ni = 18
F = np.zeros((2 * ni - 2, 2 * ni - 2))  # Initialization

py = np.zeros(ni)
pz = np.zeros(ni)

for j in range(2, ni):
    py[j] = 1.0  # In the first loop, load in the y-direction
    pz[j] = 0.0  # No loads in z

    # Call the deflection function
    uy, _, uz, _, _, _, _, _ = deflec(r, ei1, ei2, twist, pitch, py, pz)
    
    py[j] = 0.0  # Reset the load vector
    pz[j] = 0.0
    
    # Fill the F matrix for y-direction loads
    F[0:ni - 1, j - 1] = uy[1:]  # We use uy[1:] since the deformation in 1=0
    F[ni - 1:, j - 1] = uz[1:]

# Now do the same for z-direction loads
for j in range(2, ni):
    py[j] = 0.0
    pz[j] = 1.0
    
    uy, _, uz, _, _, _, _, _ = deflec(r, ei1, ei2, twist, pitch, py, pz)
    
    py[j] = 0.0
    pz[j] = 0.0
    
    # Fill the F matrix for z-direction loads
    F[0:ni - 1, ni + j - 2] = uy[1:]
    F[ni - 1:, ni + j - 2] = uz[1:]

# Setting up the mass matrix
M = np.zeros((2 * ni - 2, 2 * ni - 2))  # Mass matrix initialization

for i in range(2, ni):
    M[i - 1, i - 1] = mass[i]  # Fill the diagonal for y-direction loads
    M[ni + i - 2, ni + i - 2] = mass[i]  # Fill the diagonal for z-direction loads
