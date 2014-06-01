import sys
sys.path.append("../src")

from legendre import *
import numpy as np

print PlmSchmidt(8, np.sin(np.pi/4))