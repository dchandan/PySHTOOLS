
import numpy as np
from PySHTOOLS import SHExpandDHC, SHExpandDH, SHctor, SHrtoc

from __init__ import *


def test_SHctor():
    data = np.loadtxt(join(test_dir, "dyntopoC2E5_l1-22.txt"))
    rilm = SHExpandDH(data[:-1,:-1], 3)    # Reference real coefficients
    cilm = SHExpandDHC(data[:-1,:-1], 3)   # Complex coefficients
    rilm2 = SHctor(cilm)                   # Converted real coefficients
    assert ( ((rilm - rilm2) <= 0.00000000000001).all())
    print "PASSED"


def test_SHrtoc():
    data = np.loadtxt(join(test_dir, "dyntopoC2E5_l1-22.txt"))
    cilm = SHExpandDHC(data[:-1,:-1], 4)   # Reference Complex coefficients
    rilm = SHExpandDH(data[:-1,:-1], 4)    # Real coefficients
    cilm2 = SHrtoc(rilm)                   # Converted Complex coefficients
    assert ( ((cilm - cilm2) <= 0.00000000000001).all())
    print "PASSED"


if __name__ =="__main__":
    print("test_SHctor")
    test_SHctor()
    print("test_SHrtoc")
    test_SHrtoc()
