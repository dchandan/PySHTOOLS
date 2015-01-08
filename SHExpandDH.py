from _shtools import shexpanddh
from utilities import *


def SHExpandDH(grid, lmax_calc, degmax=None, norm=4, sampling=2, csphase=1):
    """
    Expand a grid data into its coefficients.
    ARGUMENTS
        grid      - 2D array with data
        lmax_calc - the maximum degree to decompose upto
        degmax    - degmax determines the shape of the returned array (see RETURNS)
        norm      - SH normalization
        sampling  - If 1: grid is shaped (n,n), If 2: grid is shaped (n,2n)
        csphase   - If 1, include the condon-shortley phase factor
    RETURNS
        3D array of shape (2, lmax_calc+1, lmax_calc+1) containing the sine and 
        cosine coefficients. If degmax is specified and is greater than lmax_calc then
        the shape of the returned array is (2, degmax+1, degmax+1). 
    """
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