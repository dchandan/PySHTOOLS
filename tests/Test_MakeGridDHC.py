
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.pylab as plt
from mypy.grid_toolkit import shiftgrid

from PySHTOOLS import SHExpandDHC, MakeGrid2D, SHExpandDH, MakeGridDHC

from __init__ import *


def test1():
    mpl.rcParams['contour.negative_linestyle'] = 'solid'

    data = np.loadtxt(join(test_dir, "dyntopoC2E5_l1-22.txt"))
    _deg = 8
    cilm = SHExpandDHC(data[:-1,:-1], _deg)
    grid = MakeGridDHC(cilm, _deg)
    
    grid = grid.real
    print grid.max(), grid.min()
    a, b = grid.shape
    latitudes      = np.linspace(90, -90, a)
    centeredLons   = np.linspace(-180, 180, b)
    LONS, LATS     = np.meshgrid(centeredLons, latitudes)
    
    grid = shiftgrid(grid, b/2)
    
    Figure = plt.figure(figsize=(22,10))
    Map = Basemap(projection='robin', lon_0=0, resolution='l')
    x, y = Map(LONS, LATS)
    fcp  = Map.contourf(x, y, grid, 30, interpolation="bicubic", alpha=1.0, cmap=mpl.cm.Spectral)
    levels = fcp.get_array()
    cp  = Map.contour(x, y, grid, 30, interpolation="bicubic", linewidth=0.5, colors='k', alpha=0.5)
    
    cb = Map.colorbar(fcp, "bottom", size="5%", pad='5%', extendrect=False)
    cb.ax.tick_params(labelsize=18)
    cb.solids.set_edgecolor("face")
    cb.set_label("Km",fontsize=18)
    cb.ax.set_aspect(0.047)
    
    Map.drawcoastlines(linewidth=1)
    Map.drawmapboundary(linewidth=1)
    Map.drawmeridians([-150,-100,-50,0,50,100, 150],labels=[1,1,1,0],fontsize=18)
    Map.drawparallels([-60,-30,0,30,60],labels=[1,1,1,1],fontsize=18)


if __name__ =="__main__":
    print("Test 1")
    test1()
    plt.show()