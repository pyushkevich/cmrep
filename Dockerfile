# Begin from a pre-built environment with ITK, VTK, IpOpt
FROM pyushkevich/cmrep-env
USER itk

# Copy our code
COPY --chown=itk:itk . cmrep

# Configure
RUN mkdir cmrep_build \
  && cd cmrep_build \
  && cmake \
    -D ITK_DIR:FILEPATH=/home/itk/ITKbuild \
    -D VTK_DIR:FILEPATH=/home/itk/VTKbuild \
    -D USE_IPOPT:BOOL=ON \
    -D USE_TETGEN:BOOL=ON \
    -D CMAKE_CXX_FLAGS:STRING="-std=c++11" \
    -D IPOPT_LAPACK_LIB:FILEPATH=/home/itk/CoinIpopt/install/lib/libcoinlapack.so \
    -D IPOPT_BLAS_LIB=/home/itk/CoinIpopt/install/lib/libcoinblas.so \
    -D IPOPT_GFORTRAN_LIB=/usr/lib/gcc/x86_64-linux-gnu/4.9/libgfortran.so \
    -D IPOPT_INCLUDE_DIR:PATH=/home/itk/CoinIpopt/install/include/coin \
    -D IPOPT_LIBRARY:FILEPATH=/home/itk/CoinIpopt/install/lib/libipopt.so \
    -D TETGEN_INCLUDE_DIR:PATH=/home/itk/tetgen \
    -D TETGEN_LIBRARY=/home/itk/tetgen_build/libtet.a \
    -D USE_NLOPT:BOOL=ON \
    -D NLOPT_LIBRARIES:FILEPATH=/home/itk/nlopt_build/install/lib/libnlopt.so \
    -D NLOPT_INCLUDE_DIRS:PATH=/home/itk/nlopt_build/install/include \
     ../cmrep \
  && make -j $(nproc)
