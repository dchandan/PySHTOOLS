.. _intro:

***************
Introduction
***************

PySHTOOLS is a library that brings the tools and facilities for working with spherical harmonics available in the Fortran SHTOOLS library to python. It consists of f2py'd fortran code and python wrappers, which actually provide the usage API. 

The Fortran `SHTOOLS library`_ is written by Mark Wieczorek.

.. _SHTOOLS library: http://shtools.ipgp.fr/www/conventions.html

Spherical Harmonics
====================
The definitions of spherical harmonics in PySHTOOLS can be found `here`_.

.. _here: http://shtools.ipgp.fr/www/conventions.html


.. glossary::

   cilm
      The (2, lmax+1, lmax+1) shaped array containing the spherical harmonic coefficients
