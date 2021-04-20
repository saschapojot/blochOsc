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


def S2(psiQ, q):
    '''
    2nd order Strang splitting
    :param psiQ: input wavefunction at time step q
    :param q: current time step
    :return: wavefunction at next time step, q+1
    '''
    # step 1, tilted potential term
    # mapping, coef a, coef b
    for j in range(0, N):
        psiQ[2 * j] *= np.exp(-1j * dt / 2 * xmValsAll[j])
        psiQ[2 * j + 1] *= np.exp(-1j * dt / 2 * ymValsAll[j])
    # step 2, in frequency space
    aAllTmp = [psiQ[2 * j] for j in range(0, N)]
    bAllTmp = [psiQ[2 * j + 1] for j in range(0, N)]
    # coefs of matrix F0
    vq = v((q + 1 / 2) * dt)
    wq = w((q + 1 / 2) * dt)
    uq = u((q + 1 / 2) * dt)

    AkTilde = np.fft.fft(aAllTmp)
    BkTilde = np.fft.fft(bAllTmp)
    # shift frequency
    AkShiftedAll = np.fft.fftshift(AkTilde)
    BkShiftedAll = np.fft.fftshift(BkTilde)
    # shifted frequency k/N
    freqShifted = np.fft.fftshift(np.fft.fftfreq(n=N))

    AkMappedAll = []
    BkMappedAll = []
    for k in range(0, len(freqShifted)):
        AkTmp = AkShiftedAll[k]
        BkTmp = BkShiftedAll[k]

        frqTmp = freqShifted[k]

        nk = np.sqrt(vq ** 2 + 2 * vq * wq * np.cos(2 * np.pi * frqTmp) + wq ** 2 + uq ** 2)
        z1 = (vq + wq * np.cos(2 * np.pi * frqTmp)) / nk
        z2 = wq * np.sin(2 * np.pi * frqTmp) / nk
        z3 = uq / nk
        lmdk = dt * N * nk

        AkMappedTmp = (np.cos(lmdk) - 1j * z3 * np.sin(lmdk)) * AkTmp + (
                -1j * z1 * np.sin(lmdk) - z2 * np.sin(lmdk)) * BkTmp
        BkMappedTmp = (-1j * z1 * np.sin(lmdk) + z2 * np.sin(lmdk)) * AkTmp + (
                np.cos(lmdk) + 1j * z3 * np.sin(lmdk)) * BkTmp
        AkMappedAll.append(AkMappedTmp)
        BkMappedAll.append(BkMappedTmp)

    # shift back
    ABackAll = np.fft.ifftshift(AkMappedAll)
    BBackAll = np.fft.ifftshift(BkMappedAll)
    # invserse fft
    aMappedAll = np.fft.ifft(ABackAll)
    bMappedAll = np.fft.ifft(BBackAll)

    psiq2 = []
    for j in range(0, len(aMappedAll)):
        psiq2.append(aMappedAll[j])
        psiq2.append(bMappedAll[j])

    # step 3, tilted potential term
    for j in range(0, N):
        psiq2[2 * j] *= np.exp(-1j * dt / 2 * xmValsAll[j])
        psiq2[2 * j + 1] *= np.exp(-1j * dt / 2 * ymValsAll[j])
    return psiq2


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
