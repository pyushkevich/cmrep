****************
Compiling CM-Rep
****************

These instructions should work (with minor modifications) for MacOS and Linux. I have never tried building on Windows.

Required Packages
=================

ITK
---
You will need to download and compile ITK v4.13.2. Use the following flags when compiling:

* ``CMAKE_CXX_FLAGS: -Wno-deprecated -std=c++11``
* ``CMAKE_C_FLAGS: -Wno-deprecated``

VTK
---
You will need to download and compile VTK v6.3.0. Use the following flags

* ``CMAKE_CXX_FLAGS: -Wno-deprecated -std=c++11``
* ``CMAKE_C_FLAGS: -Wno-deprecated``
* ``BUILD_SHARED_LIBS: FALSE``
* ``VTK_REQUIRED_OBJCXX_FLAGS:``
 
IPOPT and HSL
-------------
These constrained optimization libraries are required for compiling the newer (2013 IPMI and 2019 IPMI) programs, can be omitted for older ones.

* Download IPOPT and HSL libraries as described in https://www.coin-or.org/Ipopt/documentation/

* When compiling HSL, I used the following command line. Adapt to your platform::

    ./configure F77=/usr/local/Cellar/gcc/4.9.2_1/bin/gfortran-4.9 FFLAGS=-fexceptions -m64 -fbackslash F90=/usr/local/Cellar/gcc/4.9.2_1/bin/gfortran-4.9 FC=/usr/local/Cellar/gcc/4.9.2_1/bin/gfortran-4.9 --prefix=/Users/pauly/tk/ipopt/install_hsl
    make
    make install

* When compiling IPOPT, I used the following command line::

    ./configure --prefix=/Users/pauly/tk/ipopt/install64_so F77=/usr/local/Cellar/gcc/4.9.2_1/bin/gfortran FFLAGS=-fexceptions -m64 -fbackslash CFLAGS=-fno-common -no-cpp-precomp -fexceptions -arch x86_64 -m64 CFLAGS=-fno-common -no-cpp-precomp -fexceptions -arch x86_64 -m64 CXXFLAGS=-fno-common -no-cpp-precomp -fexceptions -arch x86_64 -m64 --with-hsl=/Users/pauly/tk/ipopt/install_hsl/lib/libcoinhsl.a
    make
    make install


TETGEN
------
* Download and build from http://wias-berlin.de/software/tetgen/


NLOPT
-----
This optimization library is needed for the IPMI 2019 method only
* Download from https://nlopt.readthedocs.io/en/latest/
* Compile using CMake


PARDISO
-------
Sparse solver, needed for PDE-based cm-rep programs (2005 IPMI, 2008 NeuroImage)
* Download from https://www.pardiso-project.org/
* Get license and follow instructions for how to use it


Building CM-Rep
===============

* Create a build directory and run ccmake there
* Set USE_XXX flags to enable IPOPT, HSL, NLOPT, TETGEN, PARDISO
* Set the necessary library and include paths for these directories
* Run ``make``


