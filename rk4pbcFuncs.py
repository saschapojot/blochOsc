from consts import *


# This module contains functions that are used to compute the ode with PBC+linear

def u(t):
    # return np.sin(omega * t) ** 2
    return D0 * np.cos(omega * t)


def v(t):
    return J + d0 * np.sin(omega * t)


def w(t):
    return J - d0 * np.sin(omega * t)


def x(n):
    return 2 * n * omegaF


def y(n):
    return (2 * n + 1) * omegaF


xmValsAll = [x(m) for m in range(0, N)]
ymValsAll = [y(m) for m in range(0, N)]

def F(t,an,bnm1,bn,n):
    return -1j*v(t)*bn-1j*w(t)*bnm1-1j*u(t)*an-1j*xmValsAll[n]*an

def G(t,an,anp1,bn,n):
    return -1j*v(t)*an-1j*w(t)*anp1+1j*u(t)*bn-1j*ymValsAll[n]*bn


def RK4(q,Aq,Bq):
    '''
    one step RK4
    :param q: time step q
    :param Aq: input A vector
    :param Bq: input B vector
    :return:
    '''
    #1st step in RK4
    K1=[]
    M1=[]
    for n in range(0,N):
        nPrev=(n-1)%N
        nNext=(n+1)%N
        K1.append(F(q*dt,Aq[n],Bq[nPrev],Bq[n],n))
        M1.append(G(q*dt,Aq[n],Aq[nNext],Bq[n],n))
    #2nd step in RK4
    K2=[]
    M2=[]
    for n in range(0,N):
        nPrev=(n-1)%N
        nNext=(n+1)%N
        K2.append(F((q+1/2)*dt,Aq[n]+1/2*dt*K1[n],Bq[nPrev]+1/2*dt*M1[nPrev],Bq[n]+1/2*dt*M1[n],n))
        M2.append(G((q+1/2)*dt,Aq[n]+1/2*dt*K1[n],Aq[nNext]+1/2*dt*K1[nNext],Bq[n]+1/2*dt*M1[n],n))
    #3rd step in RK4
    K3=[]
    M3=[]
    for n in range(0,N):
        nPrev=(n-1)%N
        nNext=(n+1)%N
        K3.append(F((q+1/2)*dt,Aq[n]+1/2*dt*K2[n], Bq[nPrev]+1/2*dt*M2[nPrev],Bq[n]+1/2*dt*M2[n],n))
        M3.append(G((q+1/2)*dt,Aq[n]+1/2*dt*K2[n],Aq[nNext]+1/2*dt*K2[nNext],Bq[n]+1/2*dt*M2[n],n))

    #4th step in RK4
    K4=[]
    M4=[]
    for n in range(0,N):
        nPrev=(n-1)%N
        nNext=(n+1)%N
        K4.append(F((q+1)*dt,Aq[n]+dt*K3[n],Bq[nPrev]+dt*M3[nPrev], Bq[n]+dt*M3[n],n))
        M4.append(G((q+1)*dt,Aq[n]+dt*K3[n],Aq[nNext]+dt*K3[nNext],Bq[n]+dt*M3[n],n))
    AqNext=[]
    BqNext=[]
    for n in range(0,N):
        AqNext.append(Aq[n]+1/6*dt*(K1[n]+2*K2[n]+2*K3[n]+K4[n]))
        BqNext.append(Bq[n]+1/6*dt*(M1[n]+2*M2[n]+2*M3[n]+M4[n]))
    return AqNext,BqNext





def meanXAndXWd(psiQ):
    '''

    :param psiQ: wavefunction at time q
    :return: mean position at time q
    '''
    xOut = 0
    for j in range(0, len(psiQ)):
        xOut += j * np.abs(psiQ[j]) ** 2
    return xOut


def reNormalization(vec):
    tmp2 = 0
    for elem in vec:
        tmp2 += np.abs(elem) ** 2
    tmp = np.sqrt(tmp2)
    vec /= tmp
    return vec
