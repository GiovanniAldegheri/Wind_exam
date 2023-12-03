import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

A = 8 #m/s
k = 2
P = np.array([0.7e6,2e6,3e6])       #[W]

V_low = np.array([5,10,13])
V_high = np.array([10,13,25])
V = np.arange(0,25+1,0.1) 

def Weibull():
    h = np.zeros(len(V))
    for i in range(len(h)):
        h[i] = k/A*(V[i]/A)**(k-1)*m.exp(-(V[i]/A)**k)   #Weibull probability density function
    return h

h = Weibull()

f = np.zeros(len(P))
H = np.zeros(len(P))
AEO = np.zeros(len(P))
for i in range(len(f)):
    f[i] = m.exp(-(V_low[i]/A)**k) - m.exp(-(V_high[i]/A)**k)       #[days] integrated Weibull function
    H[i] = 8760 * f[i]                                              #[hours]
    AEO[i] = P[i] * H[i]/1e6                                        #[MWh]

AEO_tot = np.sum(AEO)
print('Total AEO =',AEO_tot, '[kWh]')

plt.figure()
plt.plot(V,h, color = 'grey')
plt.title('Weibull probability density function')
plt.xlabel(r'$V_0$')
plt.ylabel(r'$h(V_0)$')
plt.fill_between(V, h, where=(V >= 0) & (V <= V_low[0]), hatch = '/', color = 'red', alpha=0.2, label = 'Cut-in speed')
for i in range((len(P))):
    plt.fill_between(V, h, where=(V >= V_low[i]) & (V <= V_high[i]), alpha=0.2, label = 'Range '+str(i+1))
plt.fill_between(V, h, where=(V >= V_high[-1]) & (V <= V[-1]), hatch = '/',color = 'red', alpha=0.4, label = 'Cut-out speed')
plt.legend(loc='upper right', fontsize=10)

plt.show()