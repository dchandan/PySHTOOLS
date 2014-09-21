from _shtools import shexpandlsq
from utilities import *


def SHExpandLSQ(d, lat, lon, lmax, norm=4, csphase=1):
    CheckNorm(norm)
    CheckCSphase(csphase)
    
    return shexpandlsq(d, lat, lon, lmax, csphase, norm=norm)
    