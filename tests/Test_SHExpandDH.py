import sys
sys.path.append("../")


import numpy as np
from mypy.system import elysium
from os.path import join
from SHExpandDH import SHExpandDH


lmax_calc = 8
degmax = 9

data = np.loadtxt("dyntopoC2E5_l1-22.txt")

cilm = SHExpandDH(data[:-1,:-1], lmax_calc, degmax=degmax)

print cilm.shape

np.set_printoptions(precision=3)
for i in range(degmax+1):
    print(cilm[0,i,:])
        
