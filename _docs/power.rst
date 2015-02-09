.. _power:


.. default-domain:: py

Calculating power
=======================

.. module:: power
   :synopsis: routines for calculating power contained in the various wavelengths


.. method:: SHPowerL(c, l)

   Calculates the power contained in a particular spherical harmonic degree

   :param c: array of SH coefficients of shape ``(2, lmax+1, lmax+1)``, where ``lmax`` is the maximum degree for which coefficients are stored in the array
   :type c: numpy array
   :param int l: the degree for which to calculate power
   :return: power in degree l


.. method:: SHPowerSpectrum(c)

   Calculates the power for all degrees represented in a CILM array

   :param c: array of SH coefficients of shape ``(2, lmax+1, lmax+1)``, where ``lmax`` is the maximum degree for which coefficients are stored in the array
   :type c: numpy array
   :return: 1D numpy array of length ``lmax``


.. method:: SHDegreeVariance(c)

   Calculate the degree variance for all degrees represented in a CILM array

   :param c: array of SH coefficients of shape ``(2, lmax+1, lmax+1)``, where ``lmax`` is the maximum degree for which coefficients are stored in the array
   :type c: numpy array
   :return: 1D numpy array of length ``lmax``


.. method:: SHPowerSpectrumDensity(c)

   Calculate the degree normalized density for all degrees represented in a CILM array

   :param c: array of SH coefficients of shape ``(2, lmax+1, lmax+1)``, where ``lmax`` is the maximum degree for which coefficients are stored in the array
   :type c: numpy array
   :return: 1D numpy array of length ``lmax``
