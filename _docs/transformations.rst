.. _transforms:

***************
Spherical Harmonic Transformations
***************

.. default-domain:: py

.. class:: SHEapandDH
.. method:: __init__(lmax_calc, degmax=None, norm=4, sampling=2, csphase=1)

   :param int lmax_calc: maximum spherical harmonic degree of decomposition
   :param degmax: maximum degree represented in the cilm array (if none, maximum degree in the array is lmax_calc)
   :param int norm: spherical harmonic normalization
   :param int sampling: the shape of the data array that will be decomposed in SH coefficients
   :param int csphase: Condon-Shortley phase factor
   :type degmax: integer or None
  
.. method:: decompose(grid)

   :param grid: the 2D grid to decompose into SH coefficients. If self.sampling == 1: (N,N) If self.sampling == 2: (N,2N)
   :type grid: numpy array
   :return: cilm array of shape (2, self.degmax+1, self.degmax+1)
   :rtype: numpy array


.. class:: MakeGrid2D
.. method:: __init__(interval, py_m, py_n, norm=4, csphase=1, north=90.0, south=-90.0, east=360.0, west=0.0)
   
   :param real interval: spacing between latitude and longitude points
   :param int py_m: number of latitude points on the output grid
   :param int py_n: number of longitude points on the output grid
   :param int norm: spherical harmonic normalization
   :param int csphase: Condon-Shortley phase factor
   :param real north: north latitude of domain
   :param real south: south latitude of domain
   :param real east: east latitude of domain
   :param real west: west latitude of domain

.. method:: compose(cilm)

   :param cilm: cilm array of SH coefficients of shape (2, self.degmax+1, self.degmax+1)
   :type cilm: numpy array
   :return: composed grid
   :rtype: numpy array
