# Start with Debian
FROM debian:latest

# Install required packages
RUN apt-get update && \
    apt-get install -y build-essential git cmake && \
    apt-get install -y python3 python3-pip && \
    apt-get clean

# Download and build ITK
RUN git clone https://github.com/Kitware/ITK.git /tk/itk/src
RUN cd /tk/itk/src && git checkout v5.2.1
RUN mkdir /tk/itk/build
WORKDIR /tk/itk/build
RUN cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=OFF \
    -DBUILD_EXAMPLES=OFF \
    /tk/itk/src
RUN make -j$(nproc)

RUN apt-get install -y mesa-common-dev libglu1-mesa-dev

# Download and build VTK
RUN git clone https://gitlab.kitware.com/vtk/vtk.git /tk/vtk/src
RUN cd /tk/vtk/src && git checkout v9.1.0
RUN mkdir /tk/vtk/build
WORKDIR /tk/vtk/build
RUN cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DVTK_USE_SYSTEM_LIBRARIES=ON \
    -DVTK_RENDERING_BACKEND=OpenGL2 \
    /tk/vtk/src
RUN make -j$(nproc)

# Install lapack
RUN apt-get install -y libopenblas-dev gfortran pkgconf

# Download VCG
RUN git clone https://github.com/cnr-isti-vclab/vcglib /tk/vcg/src
RUN cd /tk/vcg/src && git checkout 2022.02
RUN mkdir /tk/vcg/build
WORKDIR /tk/vcg/build
RUN cmake \
    -DCMAKE_BUILD_TYPE=Release \
    /tk/vcg/src
RUN make -j$(nproc)

RUN apt-get install -y libmetis-dev autoconf hwloc libhwloc-dev libudev-dev

# Try building SPRAL
RUN git clone https://github.com/ralna/spral /tk/spral/src
RUN cd /tk/spral/src && git checkout v2023.03.29
WORKDIR /tk/spral/src
RUN ./autogen.sh && mkdir /tk/spral/build
WORKDIR /tk/spral/build
RUN CFLAGS="-fPIC" CXXFLAGS="-fPIC" ../src/configure && make && make install

# Download and build IPOPT
RUN git clone https://github.com/coin-or/Ipopt /tk/ipopt/src
RUN cd /tk/ipopt/src && git checkout stable/3.14
RUN mkdir /tk/ipopt/build
WORKDIR /tk/ipopt/build
RUN CFLAGS="-fPIC" CXXFLAGS="-fPIC" ../src/configure \
    --with-spral-cflags="-I/usr/local/include" \
    --with-spral-lflags="-L/usr/local/lib -lspral -lgfortran -lmetis -lgomp -lopenblas -lstdc++ -lhwloc -fopenmp"
RUN make -j$(nproc) && make install

# Install eigen
RUN apt-get install -y libeigen3-dev

# Install qhull
RUN apt-get install -y qhull-bin libqhull-dev

# Download and build tetgen
RUN git clone https://github.com/libigl/tetgen /tk/tetgen/src
RUN mkdir /tk/tetgen/build
WORKDIR /tk/tetgen/build
RUN cmake \
    -DCMAKE_BUILD_TYPE=Release \
    /tk/tetgen/src
RUN make -j$(nproc)

# Build cmrep code
COPY . /tk/cmrep/src
RUN mkdir /tk/cmrep/build
WORKDIR /tk/cmrep/build
RUN cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DCMREP_BUILD_PDE=ON \
    -DCMREP_BUILD_BCMREP=ON \
    -DCMREP_BUILD_GSHOOT=ON \
    -DCMREP_BUILD_TETGEN_UTILS=ON \
    -DCMREP_BUILD_VCG_UTILS=ON \
    -DITK_DIR=/tk/itk/build \
    -DVTK_DIR=/tk/vtk/build \
    -DIPOPT_INCLUDE_DIR=/usr/local/include/coin-or \
    -DIPOPT_LIBRARY=/usr/local/lib/libipopt.so \
    -DVCGLIB_DIR=/tk/vcg/src \
    -DTETGEN_LIBRARY=/tk/tetgen/build/libtetgen.a \
    -DTETGEN_INCLUDE_DIR=/tk/tetgen/src \
    /tk/cmrep/src
RUN make -j$(nproc)

# Set the environment variables needed by SPRAL
ENV OMP_CANCELLATION=TRUE
ENV OMP_NESTED=TRUE
ENV OMP_PROC_BIND=TRUE

