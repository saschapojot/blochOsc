import scipy.linalg as slin
import numpy as np
from scipy import integrate
import  pandas as pd


T=1
Q=10000
dt=T/Q

def u(t):
    return np.sin(t)

def v(t):
    return np.cos(t)**2
x=1
y=2

def H(t):
    rst=np.zeros((2,2),dtype=complex)
    uVal=u(t)
    vVal=v(t)
    rst[0,0]=uVal+x
    rst[0,1]=vVal
    rst[1,0]=vVal
    rst[1,1]=y-uVal
    return rst

psi0=[1/np.sqrt(2),1j/np.sqrt(2)]
psiAll=[]
psiAll.append(psi0)
for q in range(0,Q):
    psiCurr=psiAll[q]
    psiNext=slin.expm(-1j*dt*H((q+1/2)*dt)).dot(psiCurr)
    psiAll.append(psiNext)

dataAll=pd.DataFrame(data=psiAll,columns=["a","b"])
dataAll.to_csv("S2.csv")

aL,bL=psiAll[-1]
print(np.abs(aL)**2+np.abs(bL)**2)
