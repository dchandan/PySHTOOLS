import sys
sys.path.append("../src")


import numpy as np
import matplotlib.pylab as plt
from mypy.system import elysium
from os.path import join
from SHExpandDH import SHExpandDH


data = np.loadtxt(join(elysium(), "Data1/CommonData/Dynamic_topography/C2E5/dyntopoC2E5_l1-22.txt"))

D = SHExpandDH(8, degmax=10)

cilm, lmax = D.decompose(data[:-1,:-1])
print cilm.shape
print lmax

np.set_printoptions(precision=3)
for i in range(11):
    print(cilm[0,i,:])
        
