from PySHTOOLS.utilities import *
import numpy as np


def SHPowerL(c, l):
    return (c[:,l,:]**2).sum()



def SHPowerSpectrum(c):
    CheckCILM(c.shape)
    spectra = np.zeros(c.shape[1])
    l = 0
    for i in range(len(spectra)):
        spectra[i] = SHPowerL(c, l)
        l += 1
    return spectra


def SHDegreeVariance(c):
    return np.sqrt(SHPowerSpectrum(c))



def SHPowerSpectrumDensity(c):
    spectra = SHPowerSpectrum(c)
    l = 0
    for i in range(len(spectra)):
        spectra[i] = spectra[i]/(2.*l + 1.)
    return spectra