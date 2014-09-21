from _shtools import makegrid2d
from utilities import *



def MakeGrid2D(cilm, ny, nx, norm=4, csphase=1, north=90.0, south=-90.0, east=360.0, west=0.0):
    CheckNorm(norm)
    CheckCSphase(csphase)
    
    
    # Calculating the "interval", i.e. the spacing between grid points.
    # We need to (i) know this, (ii) make sure the spacing in both dimensions are equal
    i1 = (north - south)/(ny - 1.)
    i2 = (east - west)/(nx - 1.)
    
    if not (i1 == i2):
        print("ERROR !!!")
        print(" |- In routine - PySHTOOLS::MakeGrid2D")
        print(" |->>> X and Y grid spacing must be equal, but the following spacings were calculated based on function arguments:")
        print(" |->>> {0} and {1}".format(i1, i2))
        raise PySHTOOLSError()
    
    grid, a, b = makegrid2d(cilm, i1, ny, nx, north, south, east, west, norm=norm, csphase=csphase, dealloc=0)
    return grid