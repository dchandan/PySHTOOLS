.. _intro:

***************
Installation
***************

You will require a Fortran compiler such as the freely available `gfortran <https://gcc.gnu.org/wiki/GFortran>`_ and `F2Py <http://docs.scipy.org/doc/numpy-dev/f2py/>`_. A copy of the relevant source files from SHTOOLS is included in PySHTOOLS so you *don't need* to install SHTOOLS separately in order to install PySHTOOLS. 

After `downloading <https://github.com/dchandan/PySHTOOLS>`_ and unpacking the source code follow these steps:

1. Open terminal, change to the "shtools" sub-directory within the PySHTOOLS directory.
2. Type ``make``. This will compile the Fortran codes using f2py and gfortran. If you do not have these two executables available in your ``PATH`` then this step will produce an error. Note: It is natural for f2py to produce lots of output on the screen. If the compilation is successful you will have a *_shtools.so* shared library in the top level PySHTOOLS directory. 
3. Finally to use this library, either move the complete "PySHTOOLS" folder to the site-packages directory of your python installation, or if you don't want to do that, then add the path to the directory to your ``PYTHONPATH`` shell variable. 