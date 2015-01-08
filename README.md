# PySHTOOLS
Fast, succinct and intuitive library to perform spherical harmonic analysis with Python. This is implemented as a python
wrapper around the Fortran [SHTOOLS](http://shtools.ipgp.fr) library written by Mark Wieczorek. 

**Who is this for?**
This library is aimed for users who simply want to (i) decompose a 2D field into spherical harmonic coefficients and (ii) compose a 2D field from existing coefficients. Therefore the interface has been kept very simple and inuitive. For example, if you want to restrict a 2D array upto degree and order 10 you can just do this:
```
band_limited_data = LimitBandwidth(original_array, 10)
```
If you want to also specify a lower bound for the degree, modify the above line to:
```
band_limited_data = LimitBandwidth(original_array, 10, lmin=2)
```

The fortran SHTOOLS library is a fairly vast library with lots of specialized routines. These are not included (at least in the present release) in PySHTOOLS. Thefore the "shtools" sub-directory includes only those fortran routines that contain the functionality that are necessary. These fortran files have been modified slightly from the files available in SHTOOLS inorder to assist in wrapping with f2py and in order to make the python interface as simple as possible while still retaining the flexibility afforded by SHTOOLS. The algorithm however is untouched. 

PySHTOOLS also presently only deals with arrays of real quantities; complex numbers are not supported presently. 

## Installation
To install you'll need [f2py](http://docs.scipy.org/doc/numpy-dev/f2py/), a fortran compiler and LAPACK, BLAS and FFTW libraries.

Steps:

1. Download and unpack the source code

2. On a terminal, change to the "shtools" sub-directory and type `make`. 

3. Finally to use this library, either move the complete "PySHTOOLS" folder to the site-packages 
directory of your python installation, or if you dont want to do that, then add the path to the directory
to your PYTHONPATH shell variable. 

## Examples
The tests directory contain simple tests for the library that can be used as examples.
