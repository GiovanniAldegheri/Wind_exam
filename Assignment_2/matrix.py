# %setting up the flexibility matrix, ni=18, but we only want from point 2:ni
# %in each direction, so F is square with the size of (ni-1)*2
F=zeros(2*ni-2,2*ni-2); %Initialization


for j=2:ni;
py(j)=1.0; %In the first loop I load in the y-direction
pz(j)=0.0; %No loads in z
# %My deflection script takes in the structural parameters and the loads
# %in y and z. Then it spits out all the inner loads and deformations.
# %You only need the deformations here.
[uy,thetay,uz,thetaz,mbendy,mbendz,tyshear,tzshear]=deflec(r,ei1,ei2,twist,pitch,py,pz);
py(j)=0.0; %Remember to reset the load vector before next step
pz(j)=0.0;
# % Now fill in the F matrix. If we put a load in point j, we fill in
# % column j-1. The first N-1 rows are uy and the rest are uz.
F(1:ni-1,j-1)=uy(2:end); % We only use the 2:end, since the deformation in 1=0.
F(ni:end,j-1)=uz(2:end);
end


# %Now we do exactly the same, but loading in z instead and now filling
# %columns N:(2*N)-2
for j=2:ni;
py(j)=0.0;
pz(j)=1.0;
[uy,thetay,uz,thetaz,mbendy,mbendz,tyshear,tzshear]=deflec(r,ei1,ei2,twist,pitch,py,pz);
py(j)=0.0;
pz(j)=0.0;
F(1:ni-1,ni+j-2)=uy(2:end); % We only use the 2:end, since the deformation in 1=0.
F(ni:end,ni+j-2)=uz(2:end);
end;
# %Setting mass matrix
M=zeros(2*ni-2,2*ni-2); %Mass matix is a diagonal matrix so lets start by filling it with zeros
for i=2:ni; %Filling in the diagonal. Eg. The mass in point 5 goes in M(4,4) and M(ni+3,ni+3)
M(i-1,i-1)=mass(i); %So we take directly the mass from the table and put in the diagonal
M(ni+i-2,ni+i-2)=mass(i);
end