.. _specials:

***************
Special Functionality Routines
***************

.. default-domain:: py

.. method:: LimitBandwidth(grid, lmax, lmin=0, norm=4, csphase=1)

	This convenience routine returns a grid which contains degrees ``lmin`` to ``lmax`` of the input grid. Basically, it first decoposes the input grid by making a call to SHExpandDH , then sets all coefficients outside the range to zero and then issues a call to MakeGrid2D to recompose the coefficients. 

   :param grid: A 2D grid
   :type grid: numpy array
   :param int lmax: maximum spherical harmonic degree of decomposition
   :param int lmin: minimum spherical harmonic degree of decomposition
   :param int norm: spherical harmonic normalization
   :param int csphase: Condon-Shortley phase factor
   :return: 2D grid which is band-limited between ``lmin`` and ``lmax``
   :rtype: numpy array
