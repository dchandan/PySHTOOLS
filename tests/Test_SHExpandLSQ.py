import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.pylab as plt
from mypy.grid_toolkit import shiftgrid

from PySHTOOLS import SHExpandLSQ, MakeGrid2D

from __init__ import *


def sparsify(LONS, LATS, data, l=1):
    """ l = sparseness level """
    f = 2*l
    return LONS[0:-1:f,0:-1:f], LATS[0:-1:f,0:-1:f], data[0:-1:f,0:-1:f]
    

latitudes      = np.linspace(90, -90, 181)
centeredLons   = np.linspace(-180, 180, 361)
LONS, LATS     = np.meshgrid(centeredLons, latitudes)
data = np.loadtxt(join(test_dir, "dyntopoC2E5_l1-22.txt"))

_LONS, _LATS, _data = sparsify(LONS, LATS, data, 2)

_LONS = _LONS.ravel()
_LATS = _LATS.ravel()
_data = _data.ravel()

cilm, chi2  = SHExpandLSQ(_data, _LATS, _LONS, 8)
print chi2

# Plotting part --------------------------------

latitudes      = np.linspace(90, -90, 91)
centeredLons   = np.linspace(-180, 180, 181)
LONS, LATS     = np.meshgrid(centeredLons, latitudes)

grid   = MakeGrid2D(cilm, 91, 181)

mpl.rcParams['contour.negative_linestyle'] = 'solid'
Figure = plt.figure(figsize=(22,10))
Map    = Basemap(projection='robin', lon_0=0, resolution='l')
x, y   = Map(LONS, LATS)
fcp    = Map.contourf(x, y, grid, 30, interpolation="bilinear", alpha=1.0, cmap=mpl.cm.Spectral)
levels = fcp.get_array()
cp     = Map.contour(x, y, grid, 30, interpolation="bilinear", linewidth=0.5, colors='k', alpha=0.5)
cb     = Map.colorbar(fcp, "bottom", size="5%", pad='5%', extendrect=False)
cb.ax.tick_params(labelsize=18)
cb.solids.set_edgecolor("face")
cb.set_label("Km",fontsize=18)
cb.ax.set_aspect(0.047)

Map.drawcoastlines(linewidth=1)
Map.drawmapboundary(linewidth=1)
Map.drawmeridians([-150,-100,-50,0,50,100, 150],labels=[1,1,1,0],fontsize=18)
Map.drawparallels([-60,-30,0,30,60],labels=[1,1,1,1],fontsize=18)
plt.show()