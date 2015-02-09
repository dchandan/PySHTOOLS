.. _transforms:

***************
Spherical Harmonic Transformations
***************

.. default-domain:: py

.. _shexpanddh:
.. method:: SHExpandDH(grid, lmax_calc, degmax=None, norm=4, sampling=2, csphase=1)

   This routine will expand a grid containing ``N`` samples in both longitude and latitude (or ``N x 2N``, see sampling) into spherical harmonics. This routine makes use of the sampling theorem presented in Driscoll and Healy (1994) and employs FFTs when calculating the sine and cosine terms. The number of samples, ``N``, must be **EVEN** for this routine to work, and the spherical harmonic expansion is exact if the function is band limited to degree ``N/2-1``. Legendre functions are computed on the fly using the scaling methodology presented in Holmes and Featherston (2002). When ``norm`` is 1,2 or 4, these are accurate to about degree 2800. When ``norm`` is 3, the routine is only stable to about degree 15.

   :param grid: the 2D grid to decompose into SH coefficients. If sampling == 1: (N,N) and if sampling == 2: (N,2N)
   :type grid: numpy array
   :param int lmax_calc: maximum spherical harmonic degree for which the coefficients will be calculated
   :param degmax: maximum degree represented in the cilm array (if ``degmax=None``, maximum degree in the array is ``lmax_calc``). If ``degmax > lmax_calc`` then the coefficients for degrees greater than ``lmax_calc`` are set to zero
   :param int norm: spherical harmonic normalization
   :param int sampling: the shape of the data array that will be decomposed in SH coefficients
   :param int csphase: Condon-Shortley phase factor
   :type degmax: integer or None
   :return: coefficients array of shape (2, l+1, l+1) where ``l=max(degmax, lmax_calc)``
   :rtype: numpy array


.. method:: SHExpandLSQ(d, lat, lon, lmax, norm=4, csphase=1)

   This routine will expand a set of discrete data points into spherical harmonics using a least squares inversion. When there are more data points than spherical harmonic coefficients ``(len(d) > (lmax+1)**2)`` the solution of the overdetermined system will be determined by least squares. If there are more coefficients than data points, then the solution of the underdetermined system will be determined by minimizing the solution norm. (See LAPACK DGELS documentation).

   Note that this routine takes lots of memory (~ ``8*nmax*(lmax+1)**2`` bytes) and is very slow for large ``lmax``.

   :param d: 1D array of raw data points
   :param lat: 1D array of corresponding latitudes
   :param lon: 1D array of corresponding longitudes
   :param int lmax: maximum spherical harmonic degree for which the coefficients will be calculated
   :param int norm: spherical harmonic normalization
   :param int csphase: Condon-Shortley phase factor
   :type d: numpy array
   :type lat: numpy array
   :type lon: numpy array
   :return: coefficients array of shape (2, lmax+1, lmax+1)
   :rtype: numpy array
   :return: residual sum of squares misfit for an overdetermined inversion
   :rtype: real


.. method:: MakeGrid2D(cilm, ny, nx, norm=4, csphase=1, north=90.0, south=-90.0, east=360.0, west=0.0)

   Given the Spherical Harmonic coefficients ``cilm``, this routine will compute a 2D grid with equal latitude and longitude spacings. The spacing is determined based on the bounding coordinate values represented by ``north, south, east, west`` and the number of latitude and longitude points - ``ny``, ``nx`` respectively. If the mix of these values does not produce equal spacing in both latitude and longitude, then the program will raise an error.

   The values in the returned grid are in the order of latitude going from 90 to -90 and longitudes going from 0 to 360 (both these longitudes present in the returned grid). If the optional parameters ``north``, ``south``, ``east`` and ``west`` are specified, the upper-left and lower right coordinates of the output grid are (north, west) and (south, east), respectively.

   Note that since this routined does not use FFTs, this routine is therefore slower than MakeGridDH. However MakeGridDH does not allow for the specification of arbitray dimensions of the returned grid.

   :param cilm: SH coefficients array of shape (2, lmax+1, lmax+1), where ``lmax`` represents the maximum degree for which coefficients are present in the array
   :type cilm: numpy array
   :param int ny: number of latitude points in the output grid
   :param int nx: number of longitude points in the output grid
   :param int norm: spherical harmonic normalization
   :param int csphase: Condon-Shortley phase factor
   :param real north: The north bounding latitude in the returned grid
   :param real south: The south bounding latitude in the returned grid
   :param real east: The east bounding latitude in the returned grid
   :param real west: The west bounding latitude in the returned grid
   :return: 2D grid that has been constructed from the coefficients. If ``sampling==1``, the shape of the array is (N,N) and for ``sampling==2`` the shape is (N,2N), where ``N=2*l+1``.
   :rtype: numpy array


.. method:: MakeGridDH(cilm, l, norm=4, csphase=1, sampling=2)

   :param cilm: SH coefficients array of shape (2, lmax+1, lmax+1), where ``lmax`` represents the maximum degree for which coefficients are present in the array
   :type cilm: numpy array
   :param int norm: spherical harmonic normalization
   :param int csphase: Condon-Shortley phase factor
   :param int sampling: the shape of the data array that will be decomposed in SH coefficients
   :return: 2D grid that has been constructed from the coefficients. If ``sampling==1``, the shape of the array is (N,N) and for ``sampling==2`` the shape is (N,2N), where ``N=2*l+1``.
   :rtype: numpy array
