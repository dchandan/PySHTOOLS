from _shtools import glqgridcoord

def GLQGridCoord(l):
    """
    Given a maximum spherical harmonic degree, this routine 
    will determine the latitude and longitude coordinates associated with 
    grids that are used in the Gauss-Legendre quadratue spherical harmonic 
    expansion routines. The coordinates are output in DEGREES.
    ARGUMENTS
        l - maximum spherical harmonic degree
    RETURNS
        latitudes, longitudes
        Both are 1D arrays containing 
    """
    return glqgridcoord(l)