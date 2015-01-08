from PySHTOOLS.utilities import *
import numpy as np


def SHPowerL(c, l):
    """ Calculate the power for degree l. c is the CILM array. 
    RETURNS:
        scalar value for the power
    """
    return (c[:,l,:]**2).sum()



def SHPowerSpectrum(c):
    """ Calculate the power for all degrees in the CILM array c. 
    RETURNS
        1D array
    """
    CheckCILM(c.shape)
    spectra = np.zeros(c.shape[1])
    l = 0
    for i in range(len(spectra)):
        spectra[i] = SHPowerL(c, l)
        l += 1
    return spectra


def SHDegreeVariance(c):
    """ Calculate the degree variance for all degrees in c. 
    RETURNS
        1D array
    """
    return np.sqrt(SHPowerSpectrum(c))



def SHPowerSpectrumDensity(c):
    """ Calculate the degree normalized spectrum density. 
    RETURNS
        1D array
    """
    spectra = SHPowerSpectrum(c)
    l = 0
    for i in range(len(spectra)):
        spectra[i] = spectra[i]/(2.*l + 1.)
    return spectra
