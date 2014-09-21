import sys
sys.path.append("../")


import numpy as np
from mypy.system import elysium
from os.path import join
from SHExpandDH import SHExpandDH


data = np.loadtxt("dyntopoC2E5_l1-22.txt")

cilm = SHExpandDH(data[:-1,:-1], 8, degmax=10)

# cilm, lmax = D.decompose(data[:-1,:-1])
print cilm.shape
# print lmax

np.set_printoptions(precision=3)
for i in range(11):
    print(cilm[0,i,:])
        
