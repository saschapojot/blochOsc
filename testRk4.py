import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as slin

import pandas as pd


T=1
Q=10000
dt=T/Q

def u(t):
    return np.sin(t)

def v(t):
    return np.cos(t)**2
x=1
y=2
def F(t,a,b):
    return -1j*v(t)*b-1j*u(t)*a-1j*x*a

def G(t,a,b):
    return -1j*v(t)*a+1j*u(t)*b-1j*y*b

def RK4(q,aq,bq):
    K1=F((q+1/2)*dt,aq,bq)
    M1=G((q+1/2)*dt,aq,bq)

    K2=F((q+1/2)*dt,aq+1/2*dt*K1,bq+1/2*dt*M1)
    M2=G((q+1/2)*dt,aq+1/2*dt*K1,bq+1/2*dt*M1)

    K3=F((q+1/2)*dt,aq+1/2*dt*K2,bq+1/2*dt*M2)
    M3=G((q+1/2)*dt,aq+1/2*dt*K2,bq+1/2*dt*M2)

    K4=F((q+1)*dt,aq+dt*K3,bq+dt*M3)
    M4=G((q+1)*dt,aq+dt*K3,bq+dt*M3)

    aNext=aq+1/6*(K1+2*K2+2*K3+K4)*dt
    bNext=bq+1/6*(M1+2*M2+2*M3+M4)*dt
    return aNext,bNext
AAll=[]
BAll=[]
a0=1/np.sqrt(2)
b0=1j/np.sqrt(2)
AAll.append(a0)
BAll.append(b0)
for q in range(0,Q):
    aCurr=AAll[q]
    bCurr=BAll[q]
    aNext,bNext=RK4(q,aCurr,bCurr)
    AAll.append(aNext)
    BAll.append(bNext)

dataAll={"a":AAll,"b":BAll}
df=pd.DataFrame(data=dataAll)
df.to_csv("rk4.csv")

aL=AAll[-1]
bL=BAll[-1]
print(np.abs(aL)**2+np.abs(bL)**2)