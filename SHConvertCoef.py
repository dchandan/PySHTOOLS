from _shtools import shrtoc, shctor
from utilities import *
import numpy as np


def SHctor(ccilm, degmax=None, convention=1, switchcs=0):
    s = ccilm.shape
    CheckCILM(s)
    if not (degmax == None):
        CheckLmaxNotNegative(degmax)
        if (degmax > s[1]-1):
            print("ERROR !!!")
            print("SHctor: degmax cannot be greater than the degree of ccilm")
            raise PySHTOOLSError()
    else:
        degmax = s[1] - 1


    _ccilm = np.zeros(ccilm.shape, dtype=np.float64)
    # determine the type of ccilm
    # If the user has submitted a complex array, then we need to take its
    # real and imaginary parts separately and compose them into a new array that
    # shctor requires. 
    
    if ccilm.dtype == np.dtype("complex128"):
        _ccilm[0,:,:] = ccilm[0,:,:].real
        _ccilm[1,:,:] = ccilm[0,:,:].imag
    else:
        _ccilm = np.copy(ccilm)
    
    return shctor(_ccilm, degmax, convention=convention, switchcs=switchcs)





def SHrtoc(rcilm, degmax=None, convention=1, switchcs=0):
    s = rcilm.shape
    CheckCILM(s)
    if not (degmax == None):
        CheckLmaxNotNegative(degmax)
        if (degmax > s[1]-1):
            print("ERROR !!!")
            print("SHctor: degmax cannot be greater than the degree of ccilm")
            raise PySHTOOLSError()
    else:
        degmax = s[1] - 1

    ccilm = shrtoc(rcilm, degmax, convention=convention, switchcs=switchcs)
    # The ccilm array is a REAL array containing real and imaginary compoenets
    # of the complex spherical harmonic coefficients. We need to compose this COMPLEX
    # coefficients array from the components. This is done through the relationship
    # between the complex and real components (e.g. equation in http://shtools.ipgp.fr/www/conventions_complex.html)
    _ccilm = np.zeros(ccilm.shape, dtype=np.complex128)
    _ccilm[0,:,:] = ccilm[0,:,:] + 1j*ccilm[1,:,:]
    for m in range(1, degmax+1):
        _ccilm[1,:,m] = np.conjugate(_ccilm[0,:,m]) * ((-1.0)**m)
    _ccilm[1,:,0] = 0.0 + 0.0j
    return _ccilm

