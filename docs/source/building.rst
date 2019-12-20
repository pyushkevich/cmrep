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

* Recommend downloading the **Full** HSL package (for academic use only, must receive license)

* When compiling HSL on the Mac, I used the following configuration options. Adapt the paths to your system::

    ./configure F77=/usr/local/Cellar/gcc/4.9.2_1/bin/gfortran-4.9 FFLAGS=-fexceptions -m64 -fbackslash F90=/usr/local/Cellar/gcc/4.9.2_1/bin/gfortran-4.9 FC=/usr/local/Cellar/gcc/4.9.2_1/bin/gfortran-4.9 --prefix=/Users/pauly/tk/ipopt/install_hsl
    make
    make install

* When compiling HSL on recent Linux, the following suffice::

    ./configure --prefix=/home/pauly/tk/ipopt/install_hsl
    make
    make install

* Before compiling IPOPT, enter the ``ThirdParty`` directory and run ``get.XXX`` scripts in the directories ``Blas``, ``Lapack``, ``ASL``, ``Mumps`` and ``Metis``

* When compiling IPOPT on the Mac, I used the following command line::

    ./configure --prefix=/Users/pauly/tk/ipopt/install64_so F77=/usr/local/Cellar/gcc/4.9.2_1/bin/gfortran FFLAGS=-fexceptions -m64 -fbackslash CFLAGS=-fno-common -no-cpp-precomp -fexceptions -arch x86_64 -m64 CFLAGS=-fno-common -no-cpp-precomp -fexceptions -arch x86_64 -m64 CXXFLAGS=-fno-common -no-cpp-precomp -fexceptions -arch x86_64 -m64 --with-hsl=/Users/pauly/tk/ipopt/install_hsl/lib/libcoinhsl.a
    make
    make install

* For Linux::

    ./configure --prefix=/home/pauly/tk/ipopt/install --with-hsl=/home/pauly/tk/ipopt/install_hsl
    make
    make install

TETGEN
------
* Download from http://wias-berlin.de/software/tetgen/
* Build using ``CMake``


NLOPT
-----
This optimization library is needed for the IPMI 2019 method only
* Download from https://nlopt.readthedocs.io/en/latest/
* Build using CMake


PARDISO
-------
Sparse solver, needed for PDE-based cm-rep programs (2005 IPMI, 2008 NeuroImage)
* Download from https://www.pardiso-project.org/
* Get license and follow instructions for how to use it
* The easiest is to put the library into ``/usr/local/lib``


Building CM-Rep
===============

* Create a build directory and run ccmake there
* Set the following flags to ``ON``:

  * ``USE_IPOPT`` (recommended, use for 2013 and 2019 boundary-first methods)
  * ``USE_HSL`` (recommended, use for 2013 and 2019 boundary-first methods)
  * ``USE_TETGEN``
  * ``USE_NLOPT`` (recommended, use for 2019 boundary-first method)
  * ``USE_PARDISO`` (optional, use for 2006 and 2008 PDE-based method)
  
* Set the necessary library and include paths:

  * ``ITK_DIR`` (point to ITK build directory)
  * ``VTK_DIR`` (point to VTK build directory)
  * ``IPOPT_LIBRARY`` (e.g., ``/home/pauly/tk/ipopt/install/lib/libipopt.so``)
  * ``IPOPT_INCLUDE_DIR`` (e.g., ``/home/pauly/tk/ipopt/install/include/coin``)
  * ``IPOPT_HSL_LIBRARY`` (e.g., ``/home/pauly/tk/ipopt/install_hsl/libcoinhsl.so``)
  * ``IPOPT_GFORTRAN_LIB`` (point to system ``gfortran``, e.g., ``/usr/lib/gcc/x86_64-linux-gnu/4.9/libgfortran.so``)
  * ``IPOPT_BLAS_LIB`` (point to BLAS, e.g., ``/home/pauly/tk/ipopt/install/lib/libcoinblas.so``)
  * ``TETGEN_LIBRARY:`` (e.g., ``/home/pauly/tk/tetgen/build/libtet.a``)
  * ``TETGEN_INCLUDE_DIR:`` (e.g., ``/home/pauly/tk/tetgen``)
  * ``NLOPT_LIBRARIES:`` (e.g., ``/home/pauly/tk/nlopt/build/lib/libnlopt.so``)
  * ``NLOPT_INCLUDE_DIRS:`` (e.g., ``/home/pauly/tk/nlopt/include``)

* When using PARDISO:

  * ``PARDISO_LIB`` (point to the shared library)
  * ``LAPACK_LIB`` (point to the system Lapack)
  * ``GOMP_LIB`` (point to the system GOMP: GNU OpenMP library)
  
* Set compilation flags:

  * ``CMAKE_CXX_FLAGS``: set to match those for ITK/VTK, i.e., ``-Wno-deprecated -std=c++11``
  * ``CMAKE_C_FLAGS``: set to match those for ITK/VTK, i.e., ``-Wno-deprecated``

* Run ``make``


