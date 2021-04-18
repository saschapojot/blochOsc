import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd

#total lattice number(each lattice contains 2 points A and B)
N=1000

#gaussian wavepacket parameters
#center
xc=(2*N-1)/2
# xc=101
#width
sgm=15
#momentum
k0=0

#nonlinearity
g=0
D0=2
d0=0.8
J=-1

# parameters of linear part of Hamiltonian
omega = 0.05
# tilt strength
omegaF = 0.5/3
T = 2 * np.pi / omega
Q = 2 ** 14
tTot = 3*T
dt = tTot / Q

#init gaussian part of the wavefunction
wvFcnt=[np.exp(-(n-xc)**2/(4*sgm**2)) for n in range(0,2*N)]
# plt.figure()
# plt.plot(wvFcnt)
# plt.show()
#normalization const
c02Tmp=0
for elem in wvFcnt:
    c02Tmp+=np.abs(elem)**2
C0=np.sqrt(c02Tmp)


psi0=[]
for n in range(0,2*N):
    psi0.append(wvFcnt[n]*np.exp(1j*n*k0)/C0)

L=len(psi0)

