# Requisite CMake commands
cmake_minimum_required(VERSION 3.0)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Project name
PROJECT(CMREP)
SET(CMAKE_CXX_STANDARD 11)

# Testing
ENABLE_TESTING()
INCLUDE(CTest)

# Versioning
SET(CMREP_VERSION_MAJOR 1)
SET(CMREP_VERSION_MINOR 0)
SET(CMREP_VERSION_PATCH 0)
SET(CMREP_VERSION_FULL "${CMREP_VERSION_MAJOR}.${CMREP_VERSION_MINOR}.${CMREP_VERSION_PATCH}")

# Get ITK
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

# Get VTK
FIND_PACKAGE(VTK 9 REQUIRED COMPONENTS
  ChartsCore
  CommonComputationalGeometry
  CommonCore
  CommonDataModel
  CommonExecutionModel
  CommonMath
  CommonTransforms
  FiltersCore
  FiltersGeneral
  FiltersGeometry
  FiltersModeling
  FiltersSources
  IOExport
  IOGeometry
  IOImage
  IOLegacy
  IOPLY
  ImagingCore
  ImagingGeneral
  ViewsContext2D
)

# Where to get additional modules
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")

# Set the include directories
INCLUDE_DIRECTORIES(${CMREP_SOURCE_DIR}/src ${CMREP_SOURCE_DIR}/src/dijkstra)

# Keep track of all binaries that will be compiled
SET(CMREP_BINARY_TOOLS "")

# List of libraries needed to compile different components of CMREP
SET(CMREP_LMSHOOT_LIBS "")
SET(CMREP_FIT_LIBS "")
SET(CMREP_BCMREP_LIBS "")
SET(CMREP_AUGLAG_LIBS "")
SET(CMREP_UTIL_LIBS "")
SET(CMREP_PDE_LIBS "")

# Optionally, use Eigen (provides sparse solvers)
OPTION(USE_EIGEN "Use the EIGEN sparse solvers" OFF)
IF(USE_EIGEN)
  FIND_PACKAGE (Eigen3 3.3 REQUIRED)
  IF(EIGEN3_FOUND)
    ADD_DEFINITIONS(-DHAVE_EIGEN)
    SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} Eigen3::Eigen)
    SET(CMREP_LMSHOOT_LIBS ${CMREP_LMSHOOT_LIBS} Eigen3::Eigen)
    SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} Eigen3::Eigen)
    SET(SOLVER_SRC src/EigenSolverInterface.cxx)
    INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
  ENDIF(EIGEN3_FOUND)
ENDIF(USE_EIGEN)

# Optionally, support PARDISO
OPTION(USE_PARDISO "Use the PARDISO sparse solver library" OFF)
IF(USE_PARDISO)
  IF(WIN32)
    # On Windows, the library can only be loaded dynamically
    ADD_DEFINITIONS(-DHAVE_PARDISO -DPARDISO_DYNLOAD)  
  ELSE(WIN32)

    # On Linux, Lapack is also required
    INCLUDE(${CMREP_SOURCE_DIR}/CMake/FindLapack.cmake)

    # Search for PARDISO
    FIND_LIBRARY(PARDISO_LIB pardiso DOC "The location of the PARDISO shared library")
    IF(PARDISO_LIB)
      ADD_DEFINITIONS(-DHAVE_PARDISO)
      SET(SOLVER_SRC src/PardisoInterface.cxx)

      FIND_LIBRARY(GOMP_LIB gomp DOC "The location of the libgomp library")
      IF(NOT GOMP_LIB)
        MESSAGE(ERROR "GNU OpenMP (gomp) library required for PARDISO")
      ENDIF(NOT GOMP_LIB) 

      IF(NOT LAPACK_LIB)
        MESSAGE(ERROR "Lapack is required for PARDISO")
      ENDIF(NOT LAPACK_LIB)

      # Set the PDE libraries variable
      SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} ${PARDISO_LIB} ${LAPACK_LIB} ${GOMP_LIB})
    ENDIF(PARDISO_LIB)  
  ENDIF(WIN32)
ELSE(USE_PARDISO)
  ADD_DEFINITIONS(-UHAVE_PARDISO)
ENDIF(USE_PARDISO)

# Optionally, use MKL
OPTION(USE_MKL "Use the Intel MKL library" OFF)
IF(USE_MKL)
  FIND_PACKAGE(MKL REQUIRED)
  IF(MKL_LIBRARIES)
    IF(USE_EIGEN)
      ADD_DEFINITIONS(-DHAVE_MKL)
      ADD_DEFINITIONS(-DEIGEN_USE_MKL_ALL)
      INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})

      SET(MKL_LIBRARY_DIR "")
      FOREACH(LIB ${MKL_LIBRARIES})
        GET_FILENAME_COMPONENT(LIBDIR ${LIB} DIRECTORY)
        SET(MKL_LIBRARY_DIR ${MKL_LIBRARY_DIR} ${LIBDIR})
      ENDFOREACH(LIB)

      SET(SOLVER_SRC src/EigenSolverInterface.cxx)
      SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} ${MKL_LIBRARIES})
      SET(CMREP_LMSHOOT_LIBS ${CMREP_LMSHOOT_LIBS} ${MKL_LIBRARIES})
      SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} ${MKL_LIBRARIES})
    ELSE(USE_EIGEN)
      MESSAGE(ERROR "MKL cannot be used without Eigen, set USE_EIGEN to ON")
    ENDIF(USE_EIGEN)
  ENDIF(MKL_LIBRARIES)
ENDIF(USE_MKL)

# Optionally, use BLAS/LAPACK with Eigen
OPTION(USE_LAPACK "Use the Lapack/BLAS libraries" OFF)
IF(USE_LAPACK)
  FIND_PACKAGE(LAPACK REQUIRED)
  IF(LAPACK_FOUND)
    IF(USE_EIGEN)
      ADD_DEFINITIONS(-DEIGEN_USE_BLAS)
      ADD_DEFINITIONS(-DEIGEN_USE_LAPACKE)

      SET(LAPACK_LIBRARY_DIR "")
      FOREACH(LIB ${LAPACK_LIBRARIES})
        GET_FILENAME_COMPONENT(LIBDIR ${LIB} DIRECTORY)
        SET(LAPACK_LIBRARY_DIR ${MKL_LIBRARY_DIR} ${LIBDIR})
      ENDFOREACH(LIB)

      SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} ${LAPACK_LIBRARIES})
      SET(CMREP_LMSHOOT_LIBS ${CMREP_LMSHOOT_LIBS} ${LAPACK_LIBRARIES})
      SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} ${LAPACK_LIBRARIES})
    ELSE(USE_EIGEN)
      MESSAGE(ERROR "Lapack cannot be used without Eigen, set USE_EIGEN to ON")
    ENDIF(USE_EIGEN)
  ENDIF(LAPACK_FOUND)
ENDIF(USE_LAPACK)

INCLUDE(ExternalProject)

OPTION(USE_IPOPT "Use the IpOpt Solver" OFF)
IF(USE_IPOPT)

  # Flags for compiling IPOPT and HSL
  IF(COMMENT_THIS_OUT_FOR_FUTURE)
    IF(APPLE)
      SET(IPOPT_FFLAGS "-fexceptions -m64 -fbackslash")
      SET(IPOPT_CFLAGS "-fno-common -no-cpp-precomp -fexceptions -arch x86_64 -m64 -Wno-sign-conversion")

    ENDIF(APPLE)
   
    ExternalProject_Add(project_ipopt 
      SVN_REPOSITORY "https://projects.coin-or.org/svn/Ipopt/releases/3.10.2"
      PREFIX "${CMAKE_CURRENT_BINARY_DIR}/ipopt-3.10.2"
      CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/ipopt-3.10.2
        FFLAGS=${IPOPT_CFLAGS} CFLAGS=${IPOPT_CFLAGS} CXXFLAGS=${IPOPT_CFLAGS}
      BUILD_COMMAND make
      INSTALL_COMMAND make install)

    ExternalProject_Get_Property(project_ipopt install_dir)
    SET(ipopt_install_dir ${install_dir})

    add_library(ipopt STATIC IMPORTED)
    set_property(TARGET ipopt PROPERTY IMPORTED_LOCATION ${ipopt_install_dir}/lib/libipopt.a)
    add_dependencies(ipopt project_ipopt)

    INCLUDE_DIRECTORIES(${ipopt_install_dir})
  ENDIF(COMMENT_THIS_OUT_FOR_FUTURE)

  FIND_PACKAGE(IpOpt REQUIRED)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
  MESSAGE(STATUS "Using IPOPT libraries ${IPOPT_LIBRARIES}")
  SET(CMREP_BCMREP_LIBS ${CMREP_BCMREP_LIBS} ${IPOPT_LIBRARIES})
ENDIF(USE_IPOPT)

OPTION(USE_TETGEN "Use the TetGen library" OFF)
IF(USE_TETGEN)
  FIND_PACKAGE(TetGen REQUIRED)
  SET(CMREP_BCMREP_LIBS ${CMREP_BCMREP_LIBS} ${TETGEN_LIBRARY})
ENDIF(USE_TETGEN)

# Not sure why we need OpenGL
# FIND_PACKAGE(OpenGL)

# Find NLOPT - a small non-linear optimization package
OPTION(USE_NLOPT "Use the NLOPT optimization library" OFF)
IF(USE_NLOPT)
  FIND_PACKAGE(NLOPT REQUIRED)
  INCLUDE_DIRECTORIES(${NLOPT_INCLUDE_DIRS})
  GET_FILENAME_COMPONENT(NLOPT_LIBRARY_DIR ${NLOPT_LIBRARIES} DIRECTORY)
  SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} ${NLOPT_LIBRARIES})
ENDIF(USE_NLOPT)

# VCG library (for some mesh utilities, like Taubin smooth)
OPTION(USE_VCG "Use the VCG mesh geometry library (for some utilities)")
IF(USE_VCG)
  FIND_PATH(VCGLIB_DIR NAMES vcg/simplex/vertex/base.h)
  IF(VCGLIB_DIR)
    INCLUDE_DIRECTORIES(${VCGLIB_DIR} ${VCGLIB_DIR}/eigenlib)
    ADD_COMPILE_DEFINITIONS(CMREP_USE_VCG)
  ENDIF(VCGLIB_DIR)
ENDIF(USE_VCG)

# Include subdirectories
SUBDIRS(extras/toms611)
SUBDIRS(testing)

# Add triangle to include path
INCLUDE_DIRECTORIES(extras/triangle)

OPTION(BUILD_UTILS "Build some cm-rep utilities" ON)

# Platform-specifics
IF(WIN32)
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_DEPRECATE)
ENDIF(WIN32)

# Sources for the PDE executable
SET(COMMON_SRCS
  src/BasisFunctions2D.cxx
  src/BranchingSubdivisionSurface.cxx
  src/BruteForceSubdivisionMedialModel.cxx
  src/CoefficientMapping.cxx
  src/CartesianMedialModel.cxx
  src/DiffeomorphicEnergyTerm.cxx
  src/GeometryDescriptor.cxx
  src/ITKImageWrapper.cxx
  src/JacobianDistortionPenaltyTerm.cxx
  src/MedialAtom.cxx
  src/MedialAtomGrid.cxx
  src/MeshMedialPDESolver.cxx
  src/MedialModelIO.cxx
  src/MedialPDEMasks.cxx
  src/MedialPDESites.cxx
  src/MeshTraversal.cxx
  src/OptimizationTerms.cxx
  src/OptimizationParameters.cxx
  src/PDESubdivisionMedialModel.cxx
  src/PrincipalComponents.cxx
  src/PrincipalComponentsPenaltyTerm.cxx
  src/Procrustes.cxx
  src/Registry.cxx
  src/ScriptImaging.cxx
  src/ScriptInterface.cxx
  src/SmoothedImageSampler.cxx
  src/SparseMatrix.cxx
  src/SparseSolver.cxx
  src/SubdivisionMedialModel.cxx
  src/SubdivisionSurface.cxx
  src/System.cxx
  src/TestSolver.cxx
  src/VTKMedialMesh.cxx
  src/VTKSubdivision.cxx
  ${SOLVER_SRC})

# Create the CM-REP library
ADD_LIBRARY(cmrep ${COMMON_SRCS})
TARGET_LINK_LIBRARIES(cmrep ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS cmrep MODULES ${VTK_LIBRARIES})

# Source code for the Dijkstra library
SET(DIJKSTRA_SRCS
  src/dijkstra/ShortestPath.cxx
  src/dijkstra/VTKMeshShortestDistance.cxx)

# Create library of shortest path routines
ADD_LIBRARY(cmrep_dijkstra ${DIJKSTRA_SRCS})
TARGET_LINK_LIBRARIES(cmrep_dijkstra ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS cmrep_dijkstra MODULES ${VTK_LIBRARIES})

# Set the list of required libraries
SET(CMREP_LMSHOOT_LIBS ${CMREP_LMSHOOT_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES} toms611)
SET(CMREP_BCMREP_LIBS ${CMREP_BCMREP_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES} toms611)
SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES} toms611)
SET(CMREP_UTIL_LIBS ${CMREP_UTIL_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

# New: Geodesic shooting
IF(EIGEN3_FOUND)

  ADD_EXECUTABLE(shooting_test src/PointSetHamiltonianSystem.cxx src/TestGeodesicShooting.cxx)
  TARGET_LINK_LIBRARIES(shooting_test ${CMREP_LMSHOOT_LIBS})

  ADD_EXECUTABLE(lmshoot 
    src/PointSetHamiltonianSystem.cxx 
    src/PointSetGeodesicShooting.cxx
    src/util/ReadWriteVTK.cxx)
  TARGET_LINK_LIBRARIES(lmshoot ${CMREP_LMSHOOT_LIBS})

  ADD_EXECUTABLE(lmflow
    src/PointSetHamiltonianSystem.cxx
    src/util/ReadWriteVTK.cxx
    src/PointSetGeodesicFlowUtil.cxx)
  TARGET_LINK_LIBRARIES(lmflow ${CMREP_LMSHOOT_LIBS})

  INCLUDE_DIRECTORIES(submodules/greedy/src/ITKFilters/include submodules/greedy/src)
  ADD_EXECUTABLE(lmtowarp
    src/PointSetHamiltonianSystem.cxx
    src/util/ReadWriteVTK.cxx
    src/PointSetGeodesicToWarp.cxx
    submodules/greedy/src/lddmm_data.cxx)
  TARGET_LINK_LIBRARIES(lmtowarp ${CMREP_LMSHOOT_LIBS})

  # Update the list of binaries
  SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} shooting_test lmshoot lmflow lmtowarp)

  IF(NLOPT_FOUND)
    ADD_EXECUTABLE(alcmrep_test
      src/TestAugLagCMRep.cxx
      src/MeshTraversal.cxx
      src/PointSetOptimalControlSystem.cxx
      src/PointSetHamiltonianSystem.cxx
      src/VTKMeshBuilder.cxx)
    TARGET_LINK_LIBRARIES(alcmrep_test cmrep ${CMREP_AUGLAG_LIBS})
    SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} alcmrep_test)
  ENDIF(NLOPT_FOUND)

ENDIF(EIGEN3_FOUND)

# Link the main executables
ADD_EXECUTABLE(cmrep_test src/ScriptTest.cxx)
TARGET_LINK_LIBRARIES(cmrep_test cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_fit src/Fit.cxx)
TARGET_LINK_LIBRARIES(cmrep_fit cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_getmesh src/CMRepToMesh.cxx)
TARGET_LINK_LIBRARIES(cmrep_getmesh cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_subdivide src/CMRepSubdivide.cxx)
TARGET_LINK_LIBRARIES(cmrep_subdivide cmrep  ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_remesh src/CMRepRemesh.cxx)
TARGET_LINK_LIBRARIES(cmrep_remesh cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_meanmodel src/CMRepMeanModel.cxx)
TARGET_LINK_LIBRARIES(cmrep_meanmodel cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_pca_generate src/ComputePCA.cxx)
TARGET_LINK_LIBRARIES(cmrep_pca_generate cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_fillmesh src/CMRepMeshToImage.cxx)
TARGET_LINK_LIBRARIES(cmrep_fillmesh cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_genwarp src/CMRepToWarp.cxx)
TARGET_LINK_LIBRARIES(cmrep_genwarp cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_warpimage src/CMRepWarpImage.cxx)
TARGET_LINK_LIBRARIES(cmrep_warpimage cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_afftran src/CMRepToAffineTransform.cxx)
TARGET_LINK_LIBRARIES(cmrep_afftran cmrep ${CMREP_FIT_LIBS})

ADD_EXECUTABLE(cmrep_vskel src/VoronoiSkeletonTool.cxx src/util/ReadWriteVTK.cxx)
TARGET_LINK_LIBRARIES(cmrep_vskel cmrep cmrep_dijkstra ${CMREP_FIT_LIBS})

# Update the list of binaries
SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} 
  cmrep_test cmrep_fit cmrep_getmesh cmrep_remesh cmrep_meanmodel 
  cmrep_pca_generate cmrep_fillmesh cmrep_genwarp cmrep_warpimage
  cmrep_afftran cmrep_vskel)

IF(USE_IPOPT AND USE_TETGEN)
  ADD_EXECUTABLE(bcmrep_test src/GentleNLP.cxx src/IPOptProblemInterface.cxx src/IPOptTest.cxx)
  ADD_EXECUTABLE(bcmrep_main
                 src/GentleNLP.cxx
                 src/IPOptProblemInterface.cxx
                 src/ConstrainedCMRepObjectives.cxx
                 src/ConstrainedCMRepTest.cxx
                 src/VTKMeshBuilder.cxx
                 src/QCQPProblem.h)
  TARGET_LINK_LIBRARIES(bcmrep_test cmrep ${CMREP_BCMREP_LIBS})
  TARGET_LINK_LIBRARIES(bcmrep_main cmrep ${CMREP_BCMREP_LIBS})

  SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} bcmrep_main bcmrep_test)

ENDIF(USE_IPOPT AND USE_TETGEN)

# Include the utilities
IF(BUILD_UTILS)

  # VCG-based utilities
  IF(USE_VCG)
    SET(VCG_UTIL_SRC src/util/VCGTriMesh.h src/util/VCGTriMesh.cxx ${VCGLIB_DIR}/wrap/ply/plylib.cpp)
  ENDIF(USE_VCG)

  ADD_EXECUTABLE(sphere_splat src/util/SplatOntoSphere.cxx)
  ADD_EXECUTABLE(vtkprocrustes 
    src/util/MeshProcrustes.cxx 
    src/util/ReadWriteVTK.cxx 
    src/Procrustes.cxx)

  ADD_EXECUTABLE(vtklevelset src/util/RealImageToMesh.cxx src/util/ReadWriteVTK.cxx ${VCG_UTIL_SRC})
  TARGET_LINK_LIBRARIES(vtklevelset cmrep ${CMREP_FIT_LIBS})

  ADD_EXECUTABLE(vtkcontour src/util/RealImageToContour.cxx src/util/ReadWriteVTK.cxx)

  ADD_EXECUTABLE(mesh2img src/util/FillSTLMesh.cxx src/util/DrawTriangles.cpp src/util/ReadWriteVTK.cxx)

  ADD_EXECUTABLE(mesh_image_sample src/util/MeshImageSample.cxx)

  ADD_EXECUTABLE(mesh_merge_arrays src/util/MeshMergeArrays.cxx src/util/ReadWriteVTK.cxx)

  ADD_EXECUTABLE(mesh_area_element src/util/MeshAreaElement.cxx)

  ADD_EXECUTABLE(meshglm 
    src/util/SurfaceClusterGLM.cxx 
    src/util/asa243.c
    src/util/ReadWriteVTK.cxx 
    src/Registry.cxx)

  ADD_EXECUTABLE(meshcluster 
    src/util/SurfaceClusterAnalysis.cxx 
    src/util/ReadWriteVTK.cxx 
    src/Registry.cxx)

  ADD_EXECUTABLE(meshdiff src/util/CompareMeshes.cxx src/util/DrawFillInside.cpp src/util/DrawTriangles.cpp)
  #TARGET_LINK_LIBRARIES(meshdiff ${ITK_LIBRARIES} ${VTK_LIBRARIES})

  ADD_EXECUTABLE(warpmesh src/util/WarpMeshBackwards.cxx)

  # This is a one-off tool for sampling from skeleton tetra meshes
  ADD_EXECUTABLE(skel_tetra_max_thick_sphere src/util/SkelTetraMaxThicknessSphere.cxx)

  # This uses triangle library to triangulate a 2d contour
  ADD_EXECUTABLE(contour2surf src/util/ContourMesh2SurfMesh.cxx src/util/ReadWriteVTK.cxx)

  ADD_EXECUTABLE(apply_remeshing_to_another src/util/ApplyRemeshingToAnotherMesh.cxx src/util/ReadWriteVTK.cxx)

  ADD_EXECUTABLE(mesh_convert_format src/util/MeshConvertFormat.cxx src/util/ReadWriteVTK.cxx)

  ADD_EXECUTABLE(ml_affine src/util/MultiLabelAffine.cxx src/GentleNLP.cxx)

  # List of cmrep utilities
  SET(CMREP_UTIL_TOOLS
    sphere_splat vtkprocrustes vtklevelset vtkcontour
    mesh2img mesh_image_sample mesh_merge_arrays mesh_area_element
    meshglm meshcluster meshdiff warpmesh skel_tetra_max_thick_sphere
    contour2surf apply_remeshing_to_another mesh_convert_format ml_affine)

  # Tetrahedron tools
  ADD_EXECUTABLE(tetdumpattr src/util/TetDumpAttr.cxx)
  ADD_EXECUTABLE(tetfill src/util/TetraFill.cxx)
  ADD_EXECUTABLE(tetjac src/util/TetJacobian.cxx)
  ADD_EXECUTABLE(tetmergearr src/util/MergeVTKArrays.cxx)
  ADD_EXECUTABLE(tetsample src/util/TetSampleImage.cxx)
  SET(CMREP_UTIL_TOOLS ${CMREP_UTIL_TOOLS}
    tetdumpattr tetfill tetjac tetmergearr tetsample)

  IF(USE_VCG)
    ADD_EXECUTABLE(mesh_smooth_curv src/util/MeshSmoothAndCurvature.cxx ${VCG_UTIL_SRC} )
    ADD_EXECUTABLE(mesh_poisson_sample src/util/MeshPoissonSample.cxx ${VCG_UTIL_SRC} )
    SET(CMREP_UTIL_TOOLS ${CMREP_UTIL_TOOLS} mesh_smooth_curv mesh_poisson_sample)
  ENDIF(USE_VCG)

  # Link libraries for the tools
  FOREACH(CMREP_TOOL ${CMREP_UTIL_TOOLS})
    vtk_module_autoinit(TARGETS ${CMREP_TOOL} MODULES ${VTK_LIBRARIES})
    TARGET_LINK_LIBRARIES(${CMREP_TOOL} ${CMREP_UTIL_LIBS})
  ENDFOREACH(CMREP_TOOL)

  # Add to the list of binaries
  SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} ${CMREP_UTIL_TOOLS})

  # TetGen tetra tools
  IF(USE_TETGEN)
    ADD_EXECUTABLE(tetgen2vtk 
      src/util/TetGenToVTK.cxx src/util/ReadWriteVTK.cxx src/util/TetraLib.cxx)
    ADD_EXECUTABLE(surf2tet
      src/util/SurfMesh2VolMesh.cxx src/util/ReadWriteVTK.cxx src/util/TetraLib.cxx)
    TARGET_LINK_LIBRARIES(tetgen2vtk ${CMREP_UTIL_LIBS} ${TETGEN_LIBRARY})
    TARGET_LINK_LIBRARIES(surf2tet ${CMREP_UTIL_LIBS} ${TETGEN_LIBRARY})
    SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} tetgen2vtk surf2tet)
  ENDIF(USE_TETGEN)

ENDIF(BUILD_UTILS)

# Install all the targets
FOREACH(CMREP_TOOL ${CMREP_BINARY_TOOLS})

  INSTALL(TARGETS ${CMREP_TOOL} 
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT RUNTIME)

ENDFOREACH(CMREP_TOOL)

INSTALL(DIRECTORY ${CMREP_SOURCE_DIR}/examples DESTINATION
  "share/cmrep/")

# INCLUDE(BundleUtilities)

# Fixup bundle on MACOS
IF(APPLE)
  SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/cmrep_fit")
  INSTALL(CODE "
     function(gp_item_default_embedded_path_override item path) 
       set(path @executable_path PARENT_SCOPE) 
     endfunction() 
     include(BundleUtilities)
     fixup_bundle(\"${APPS}\"   \"\"   \"\")
     " COMPONENT Runtime)
ENDIF(APPLE)

# Fixup bundle on LINUX
IF(UNIX AND NOT APPLE)

  SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/cmrep_fit")
  SET(FIXUP_DIRS "${IPOPT_LIB_DIR}" "${VTK_DIR}/lib" "${NLOPT_LIBRARY_DIR}" "${MKL_LIBRARY_DIR}" "${LAPACK_LIBRARY_DIR}")
  MESSAGE(STATUS "Fixup dirs: ${FIXUP_DIRS}")

  INSTALL(CODE "
     include(BundleUtilities)
     fixup_bundle(\"${APPS}\"   \"\"   \"${FIXUP_DIRS}\")
     " COMPONENT Runtime)

ENDIF(UNIX AND NOT APPLE)


# Add package support
INCLUDE(CMake/Package.cmake)



