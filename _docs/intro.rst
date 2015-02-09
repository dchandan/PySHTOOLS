.. _intro:

***************
Introduction
***************

PySHTOOLS exposes spherical harmonics related functions available in the Fortran SHTOOLS library to python. `SHTOOLS`_ is written by Mark Wieczorek.

PySHTOOLS has been made possible by F2PY which has allowed for quick and painless wrapping of the original Fortran code.

.. _SHTOOLS: http://shtools.ipgp.fr/www/conventions.html

Spherical Harmonics
====================
The definitions for spherical harmonics used in PySHTOOLS are the same as those SHTOOLS and these can found `here`_.

.. _here: http://shtools.ipgp.fr/www/conventions.html


.. glossary::

   cilm
      The ``(2, lmax+1, lmax+1)`` shaped array containing the spherical harmonic coefficients. ``lmax`` is the maximum degree for which the coefficients are represented in the array.
