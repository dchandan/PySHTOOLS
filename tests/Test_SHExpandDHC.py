
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.pylab as plt
from mypy.grid_toolkit import shiftgrid

from PySHTOOLS import SHExpandDHC, MakeGrid2D, SHExpandDH

from __init__ import *


def test1():
    mpl.rcParams['contour.negative_linestyle'] = 'solid'

    data = np.loadtxt(join(test_dir, "dyntopoC2E5_l1-22.txt"))
    cilm = SHExpandDHC(data[:-1,:-1], 4)
    print data.dtype
    print "Real"
    print cilm.real
    print "Imag"
    print cilm.imag    


if __name__ =="__main__":
    print("Test 1")
    test1()
