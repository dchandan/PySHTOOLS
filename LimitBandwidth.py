import numpy as np
from PySHTOOLS import SHExpandDH, MakeGrid2D


def LimitBandwidth(grid, lmax, lmin=0, norm=4, csphase=1):
    """
    Restricts the input field to a maximum spherical harmonic degree
    ARGUMENTS
        grid - 2D grid of real data type
        lmax - the maximum degree in the reutrned field
        lmin - (optional) minimum degree. If it is 0 then all degrees
               upto lmax are used. Otherwise degrees below lmin will not be included. 
        norm - SH normalization
        csphase - condon-shortley phase factor
    RETURNS
        A 2D grid. 
    """

    if (lmin < 0):
        print "lmin cannot be less than 0"
        raise ValueError()
    if (lmin > lmax):
        print "lmin cannot be greater than lmax"
        raise ValueError()
    
    cilm = SHExpandDH(grid, lmax)
    if (lmin > 0):
        for i in range(lmin+1):
            cilm[:,i,:] = 0.0
    return MakeGrid2D(cilm, 181, 361)
