
import numpy as np
from PySHTOOLS import SHExpandDHC, SHExpandDH, SHctor, SHrtoc
from PySHTOOLS import SHCilmToCindex, SHCindexToCilm

from __init__ import *


def test_SHctor():
    deg = 30
    data = np.loadtxt(join(test_dir, "dyntopoC2E5_l1-22.txt"))
    rilm = SHExpandDH(data[:-1,:-1], deg)    # Reference real coefficients
    cilm = SHExpandDHC(data[:-1,:-1], deg)   # Complex coefficients
    rilm2 = SHctor(cilm)                   # Converted real coefficients
    assert ( ((rilm - rilm2) <= 0.00000000000001).all())
    print "PASSED"


def test_SHrtoc():
    deg = 30
    data = np.loadtxt(join(test_dir, "dyntopoC2E5_l1-22.txt"))
    cilm = SHExpandDHC(data[:-1,:-1], deg)   # Reference Complex coefficients
    rilm = SHExpandDH(data[:-1,:-1], deg)    # Real coefficients
    cilm2 = SHrtoc(rilm)                   # Converted Complex coefficients
    assert ( ((cilm - cilm2) <= 0.00000000000001).all())
    print "PASSED"

def test_cilm_cindex():
    deg = 30
    data = np.loadtxt(join(test_dir, "dyntopoC2E5_l1-22.txt"))
    cilm = SHExpandDH(data[:-1,:-1], deg)
    cindex = SHCilmToCindex(cilm)
    cilm2 = SHCindexToCilm(cindex, deg)
    assert( (cilm == cilm2).all())
    print "PASSED"



if __name__ =="__main__":
    print("test_SHctor")
    test_SHctor()
    print("test_SHrtoc")
    test_SHrtoc()
    print("test_cilm_cindex")
    test_cilm_cindex()
