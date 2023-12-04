import numpy as np

a = np.arange(0.641,0.642,0.00001)
CT = np.zeros_like(a)
CT = 4*a*(1-a)
CT = 4*a*(1-(5-3*a)*a/4)

for i in range(len(a)):
    print('a:',a[i],'   CT:',CT[i])