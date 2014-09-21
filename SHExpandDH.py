from _shtools import shexpanddh
from utilities import *


def SHExpandDH(grid, lmax_calc, degmax=None, norm=4, sampling=2, csphase=1):
    s = grid.shape
    if (degmax == None):
        py_lmax = lmax_calc
    else:
        py_lmax = degmax
        if (py_lmax < lmax_calc):
            print("ERROR !!!")
            print(" |- In routine - PySHTOOLS::SHExpandDH")
            print(" |->>> degmax cannot be less than lmax_calc")
            raise PySHTOOLSError()
    
    CheckNorm(norm)
    CheckCSphase(csphase)
    CheckSampling1(sampling)
    CheckSampling2(sampling, s)
    # cilm, lmax = shexpanddh(grid, s[0], lmax_calc, py_lmax, sampling=sampling, csphase=csphase, norm=norm)
    cilm, lmax = shexpanddh(grid, lmax_calc, py_lmax, sampling=sampling, csphase=csphase, norm=norm)
    return cilm