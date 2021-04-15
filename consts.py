import numpy as np
from datetime import datetime



#total lattice number(each lattice contains 2 points A and B)
N=2**10

#gaussian wavepacket parameters
#center
xc=(2*N-1)/2
#width
sgm=0.05
#momentum
k0=0

#nonlinearity
g=0


# parameters of linear part of Hamiltonian
omega = 0.02
# tilt strength
omegaF = 0
T = 2 * np.pi / omega
Q = 2 ** 11
tTot = T
dt = tTot / Q

#init gaussian part of the wavefunction
wvFcnt=[np.exp(-(n-xc)**2/sgm) for n in range(0,2*N)]
#normalization const
c02Tmp=0
for elem in wvFcnt:
    c02Tmp+=np.abs(elem)**2
C0=np.sqrt(c02Tmp)


psi0=[]
for n in range(0,2*N):
    psi0.append(wvFcnt[n]*np.exp(1j*n*k0)/C0)

L=len(psi0)
