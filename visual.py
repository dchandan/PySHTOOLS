from mayavi import mlab
import numpy as np
from scipy.special import sph_harm
from PySHTOOLS import MakeGrid2D

def SPHarmPlot3D(lm):
    r   = 1
    pi  = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:181j, 0:2 * pi:361j]

    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)

    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(500, 500))
    mlab.clf()

    s = np.zeros(phi.shape)
    for l, m in lm:
        ss = sph_harm(m, l, theta, phi).real
        s += ss
    s[s < 0] *= 0.97
    mlab.mesh(s * x, s * y, s * z, scalars=s, colormap='Spectral')
    mlab.show()


def Plot3D(l, m):
    r   = 1
    pi  = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:181j, 0:2 * pi:361j]

    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)
    
    cilm = np.zeros((2, l+1, l+1))
    if m >= 0:
        cilm[0, l, m] = 1.0
    else:
        cilm[1, l, m] = 1.0
    Composer = MakeGrid2D(1.0, 181, 361)
    grid, a, b = Composer.compose(cilm)
    
    
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(500, 500))
    mlab.clf()

    mlab.mesh(grid * x, grid * y, grid * z, scalars=grid, colormap='Spectral')
    mlab.show()


SPHarmPlot3D([[4,3]])