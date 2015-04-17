from _shtools import makegriddhc
from utilities import *



def MakeGridDHC(cilm, l, norm=4, csphase=1, sampling=2):
    """
    Given the COMPLEX Spherical Harmonic coefficients CILM, this subroutine
    will evalate the function on a grid with an equal number of samples N in 
    both latitude and longitude (or N by 2N by specifying the optional parameter
    SAMPLING = 2). This is the inverse of the routine SHExpandDHC, both of which
    are done quickly using FFTs for each degree of each latitude band. 
    ARGUMENTS
        cilm     - 3D Coefficients array
        l        - the degree upto which to include in the constructed grid
        norm     - the spherical harmonic normalization to use
        csphase  - 
        sampling - if equal to 1, then grid in NxN else Nx2N, where N=2*l+2
    
    """
    CheckNorm(norm)
    CheckCSphase(csphase)
    CheckSampling1(sampling)

    s = cilm.shape
    if (l > s[1]-1):
        print("ERROR !!!")
        print(" |- In routine - PySHTOOLS::MakeGridDH")
        print(" |->>> Degree requested is greater than maximum degree in CILM")
        raise PySHTOOLSError()

    ny = 2*l + 2
    if (sampling == 1):
        nx = ny
    elif (sampling == 2):
        nx = 2*ny
    
    grid = makegriddhc(ny, nx, cilm, csphase=csphase, norm=norm, sampling=sampling)
    return grid