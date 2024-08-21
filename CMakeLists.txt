# Requisite CMake commands
cmake_minimum_required(VERSION 3.0)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Project name
PROJECT(CMREP)
SET(CMAKE_CXX_STANDARD 17)

# Testing
ENABLE_TESTING()
INCLUDE(CTest)

# Versioning
SET(CMREP_VERSION_MAJOR 1)
SET(CMREP_VERSION_MINOR 0)
SET(CMREP_VERSION_PATCH 0)
SET(CMREP_VERSION_FULL "${CMREP_VERSION_MAJOR}.${CMREP_VERSION_MINOR}.${CMREP_VERSION_PATCH}")

# Variables that determine the needed packages
SET(NEED_EIGEN OFF)
SET(NEED_IPOPT OFF)
SET(NEED_LAPACK OFF)
SET(NEED_MKL OFF)
SET(NEED_NLOPT OFF)
SET(NEED_PARDISO OFF)
SET(NEED_TETGEN OFF)
SET(NEED_VCG OFF)
SET(NEED_QHULL OFF)

# Set of options that users can select to build different components of cm-rep
# Building PDE requires a sparse solver
OPTION(CMREP_BUILD_PDE "Build cm-rep models formulated as a solution of a PDE (Neuroimage '08)" OFF)
IF(CMREP_BUILD_PDE)
SET(CMREP_PDE_SPARSE_SOLVER "EIGEN" CACHE STRING "Sparse linear solver to use for PDE-based cm-rep model")
SET_PROPERTY(CACHE CMREP_PDE_SPARSE_SOLVER PROPERTY STRINGS EIGEN PARDISO MKL_PARDISO)
IF(CMREP_PDE_SPARSE_SOLVER STREQUAL EIGEN)
SET(NEED_EIGEN ON)
SET(NEED_LAPACK ON)
ELSEIF(CMREP_PDE_SPARSE_SOLVER STREQUAL PARDISO)
SET(NEED_PARDISO ON)
IF(NOT WIN32)
SET(NEED_LAPACK ON)
ENDIF()
ELSEIF(CMREP_PDE_SPARSE_SOLVER STREQUAL MKL_PARDISO)
SET(NEED_EIGEN ON)
SET(NEED_MKL ON)
ENDIF()
ENDIF()

# Building the bcmrep cm-rep model requires IPOpt and a sparse solver
OPTION(CMREP_BUILD_BCMREP "Build boundary-constrained cm-rep model (IPMI '13)" OFF)
IF(CMREP_BUILD_BCMREP)
SET(NEED_IPOPT ON)
SET(NEED_LAPACK ON)
ENDIF()

# Building AUGLAG requires a separate optimization library
OPTION(CMREP_BUILD_AUGLAG "Build augmented lagrangean diffeomorphic cm-rep model (IPMI '19)" OFF)
IF(CMREP_BUILD_AUGLAG)
SET(NEED_EIGEN ON)
SET(NEED_NLOPT ON)
ENDIF()

# Build tet-gen utilities
OPTION(CMREP_BUILD_TETGEN_UTILS "Build utility programs for mesh tetrahedralization" OFF)
IF(CMREP_BUILD_TETGEN_UTILS)
SET(NEED_TETGEN ON)
ENDIF()

# Build tet-gen utilities
OPTION(CMREP_BUILD_VCG_UTILS "Build utility programs for that use the VCG library" OFF)
IF(CMREP_BUILD_VCG_UTILS)
SET(NEED_VCG ON)
ENDIF()

# Build Voronoi skeletonization tool
OPTION(CMREP_BUILD_VSKEL "Build Voronoi skeletonization tool" ON)
IF(CMREP_BUILD_VSKEL)
SET(NEED_QHULL ON)
ENDIF()

# Print the current configuration
MESSAGE(STATUS "Libraries required by the current configuration:")
MESSAGE(STATUS "  EIGEN: ${NEED_EIGEN}")
MESSAGE(STATUS "  IPOPT: ${NEED_IPOPT}")
MESSAGE(STATUS "  LAPACK: ${NEED_LAPACK}")
MESSAGE(STATUS "  MKL: ${NEED_MKL}")
MESSAGE(STATUS "  NLOPT: ${NEED_NLOPT}")
MESSAGE(STATUS "  PARDISO: ${NEED_PARDISO}")
MESSAGE(STATUS "  QHULL: ${NEED_QHULL}")
MESSAGE(STATUS "  TETGEN: ${NEED_TETGEN}")
MESSAGE(STATUS "  VCG: ${NEED_VCG}")

# Get ITK
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

# Get VTK
FIND_PACKAGE(VTK 9 REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  CommonMath
  CommonTransforms
  FiltersCore
  FiltersGeneral
  FiltersGeometry
  FiltersModeling
  FiltersSources
  FiltersVerdict
  IOGeometry
  IOImage
  IOLegacy
  IOPLY)

# Where to get additional modules
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")

# Set the include directories
INCLUDE_DIRECTORIES(${CMREP_SOURCE_DIR}/src ${CMREP_SOURCE_DIR}/src/dijkstra)

# Keep track of all binaries and libraries that will be compiled
SET(CMREP_BINARY_TOOLS "")
SET(CMREP_API_LIBRARIES "")

# List of libraries needed to compile different components of CMREP
SET(CMREP_LMSHOOT_LIBS "")
SET(CMREP_FIT_LIBS "")
SET(CMREP_BCMREP_LIBS "")
SET(CMREP_AUGLAG_LIBS "")
SET(CMREP_UTIL_LIBS "")
SET(CMREP_PDE_LIBS "")

# If we need Eigen, then find Eigen
IF(NEED_EIGEN)
FIND_PACKAGE (Eigen3 3.3 REQUIRED)
IF(EIGEN3_FOUND)
ADD_DEFINITIONS(-DHAVE_EIGEN)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
ENDIF()
ENDIF()


# Optionally, support PARDISO
IF(NEED_PARDISO)
IF(WIN32)
# On Windows, the library can only be loaded dynamically
ADD_DEFINITIONS(-DHAVE_PARDISO -DPARDISO_DYNLOAD)
SET(PARDISO_FOUND ON BOOL)
ELSE()
# Search for PARDISO
FIND_LIBRARY(PARDISO_LIBRARY pardiso DOC "The location of the PARDISO shared library")
IF(PARDISO_LIBRARY)
FIND_LIBRARY(GOMP_LIB gomp DOC "The location of the libgomp library")
IF(GOMP_LIB)
ADD_DEFINITIONS(-DHAVE_PARDISO)
SET(PARDISO_FOUND ON BOOL)
SET(PARDISO_LIBS ${PARDISO_LIBRARY} ${GOMP_LIB})
ELSE()
MESSAGE(ERROR "GNU OpenMP (gomp) library required for PARDISO")
ENDIF()
ENDIF()
ENDIF()
ELSE()
ADD_DEFINITIONS(-UHAVE_PARDISO)
ENDIF()

# Optionally, use MKL
IF(NEED_MKL)
FIND_PACKAGE(MKL REQUIRED)
IF(MKL_LIBRARIES)
ADD_DEFINITIONS(-DHAVE_MKL)
INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})

SET(MKL_LIBRARY_DIR "")
FOREACH(LIB ${MKL_LIBRARIES})
GET_FILENAME_COMPONENT(LIBDIR ${LIB} DIRECTORY)
SET(MKL_LIBRARY_DIR ${MKL_LIBRARY_DIR} ${LIBDIR})
ENDFOREACH(LIB)

# TODO: is this necessary?
IF(EIGEN3_FOUND)
ADD_DEFINITIONS(-DEIGEN_USE_MKL_ALL)

#SET(CMREP_LMSHOOT_LIBS ${CMREP_LMSHOOT_LIBS} ${MKL_LIBRARIES})
#SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} ${MKL_LIBRARIES})
ENDIF()
ENDIF(MKL_LIBRARIES)
ENDIF()

# Optionally, use BLAS/LAPACK with Eigen
IF(NEED_LAPACK)

# Provide an option to use MKL as the source of Lapack functionality
OPTION(CMREP_LAPACK_USE_MKL OFF "Use Intel MKL as source of Lapack functionality")
IF(NEED_MKL)
SET(CMREP_LAPACK_USE_MKL ON)
ENDIF()
IF(CMREP_LAPACK_USE_MKL)
SET(NEED_MKL ON)
SET(LAPACK_LIBRARIES ${MKL_LIBRARIES})
ELSE()
# If Lapack taken outside of MKL, it needs to be found
FIND_PACKAGE(LAPACK REQUIRED)
IF(LAPACK_FOUND)
SET(LAPACK_LIBRARY_DIR "")
FOREACH(LIB ${LAPACK_LIBRARIES})
GET_FILENAME_COMPONENT(LIBDIR ${LIB} DIRECTORY)
SET(LAPACK_LIBRARY_DIR ${MKL_LIBRARY_DIR} ${LIBDIR})
ENDFOREACH(LIB)

# TODO: what about this?
IF(EIGEN3_FOUND)
ADD_DEFINITIONS(-DEIGEN_USE_BLAS)
# ADD_DEFINITIONS(-DEIGEN_USE_LAPACKE)
ENDIF(EIGEN3_FOUND)
ENDIF()
ENDIF()
ENDIF()

# Set a variable that includes libraries that can be linked with for programs that
# need basic LAPACK functionality


INCLUDE(ExternalProject)

IF(NEED_IPOPT)
# Find the IPOpt library itself
FIND_LIBRARY(IPOPT_LIBRARY NAMES ipopt REQUIRED)

# Find the include files
FIND_PATH(IPOPT_INCLUDE_DIR NAMES IpTNLP.hpp PATHS ${IPOPT_ROOT_DIR}/include/coin REQUIRED)
INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})

# Set the IPOpt libraries
SET(IPOPT_LIBS ${IPOPT_LIBRARY})

# FIND_PACKAGE(IpOpt REQUIRED)
# MESSAGE(STATUS "Using IPOPT libraries ${IPOPT_LIBRARIES}")
# SET(CMREP_BCMREP_LIBS ${CMREP_BCMREP_LIBS} ${IPOPT_LIBRARIES})
ENDIF()

IF(NEED_TETGEN)
FIND_PACKAGE(TetGen REQUIRED)
ENDIF()

# Not sure why we need OpenGL
# FIND_PACKAGE(OpenGL)

# Find NLOPT - a small non-linear optimization package
IF(NEED_NLOPT)
FIND_PACKAGE(NLOPT REQUIRED)
INCLUDE_DIRECTORIES(${NLOPT_INCLUDE_DIRS})
GET_FILENAME_COMPONENT(NLOPT_LIBRARY_DIR ${NLOPT_LIBRARIES} DIRECTORY)
# SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} ${NLOPT_LIBRARIES})
ENDIF()

# VCG library (for some mesh utilities, like Taubin smooth)
IF(NEED_VCG)
FIND_PATH(VCGLIB_DIR NAMES vcg/simplex/vertex/base.h)
IF(VCGLIB_DIR)
INCLUDE_DIRECTORIES(${VCGLIB_DIR} ${VCGLIB_DIR}/eigenlib)
ADD_COMPILE_DEFINITIONS(CMREP_USE_VCG)
ENDIF(VCGLIB_DIR)
ENDIF()

# Qhull
IF(NEED_QHULL)
  INCLUDE(CMake/FindQhull.cmake)
  SET(QHULL_LIBRARY Qhull::qhullstatic_r)
ENDIF()

# Include subdirectories
SUBDIRS(testing)

# Add triangle to include path
INCLUDE_DIRECTORIES(extras/triangle)

# Platform-specifics
IF(WIN32)
ADD_DEFINITIONS(-D_SCL_SECURE_NO_DEPRECATE)
ENDIF(WIN32)

# Create a library for toms611
ADD_LIBRARY(cmrep_toms611 extras/toms611/toms611.c)

# Determine which solver wrapper to build for PDE
SET(CMREP_FIT_LIBS ${ITK_LIBRARIES} ${VTK_LIBRARIES} cmrep_toms611)
IF(CMREP_BUILD_PDE)
IF(CMREP_PDE_SPARSE_SOLVER STREQUAL "EIGEN")
SET(SOLVER_SRC ${SOLVER_SRC} src/EigenSolverInterface.cxx)
SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} Eigen3::Eigen ${})
ELSEIF(CMREP_PDE_SPARSE_SOLVER STREQUAL "PARDISO")
SET(SOLVER_SRC ${SOLVER_SRC} src/PardisoInterface.cxx)
SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} ${PARDISO_LIBS} ${LAPACK_LIBRARIES})
ELSEIF(CMREP_PDE_SPARSE_SOLVER STREQUAL "MKL_PARDISO")
SET(SOLVER_SRC ${SOLVER_SRC} src/EigenSolverInterface.cxx)
SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} Eigen3::Eigen ${MKL_LIBRARIES})
ENDIF()
ENDIF()

#SET(CMREP_FIT_LIBS ${CMREP_FIT_LIBS} ${LAPACK_LIBRARIES})
#SET(CMREP_LMSHOOT_LIBS ${CMREP_LMSHOOT_LIBS} ${LAPACK_LIBRARIES})
#SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} ${LAPACK_LIBRARIES})


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
  src/SubdivisionMedialModel.cxx
  src/SubdivisionSurface.cxx
  src/System.cxx
  src/VTKMedialMesh.cxx
  src/VTKSubdivision.cxx
  src/SparseSolver.cxx
  src/TestSolver.cxx
  src/MeshMedialPDESolver.cxx
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

# Add these to the list of all libraries we will export
SET(CMREP_API_LIBRARIES ${CMREP_API_LIBRARIES} cmrep cmrep_dijkstra cmrep_toms611)

# Set the list of required libraries
SET(CMREP_LMSHOOT_LIBS ${CMREP_LMSHOOT_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
SET(CMREP_BCMREP_LIBS ${CMREP_BCMREP_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES} cmrep_toms611)
SET(CMREP_AUGLAG_LIBS ${CMREP_AUGLAG_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES} cmrep_toms611)
SET(CMREP_UTIL_LIBS ${CMREP_UTIL_LIBS} ${ITK_LIBRARIES} ${VTK_LIBRARIES})

# Augmented Lagrangean cm-rep
IF(CMREP_BUILD_AUGLAG)
SET(CMREP_LMSHOOT_LIBS Eigen3::Eigen ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${LAPACK_LIBRARIES} ${NLOPT_LIBRARIES})

ADD_EXECUTABLE(alcmrep_test
  src/TestAugLagCMRep.cxx
  src/MeshTraversal.cxx
  src/PointSetOptimalControlSystem.cxx
  submodules/greedy/src/lmshoot/PointSetHamiltonianSystem.cxx
  src/VTKMeshBuilder.cxx)
TARGET_LINK_LIBRARIES(alcmrep_test cmrep ${CMREP_AUGLAG_LIBS})
SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} alcmrep_test)
ENDIF()

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

# Update the list of binaries
SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS}
  cmrep_test cmrep_fit cmrep_getmesh cmrep_remesh cmrep_meanmodel
  cmrep_pca_generate cmrep_fillmesh cmrep_genwarp cmrep_warpimage
  cmrep_afftran)

IF(CMREP_BUILD_VSKEL)
ADD_LIBRARY(cmrep_vskel_api src/VoronoiSkeletonTool.cxx src/util/ReadWriteVTK.cxx)
TARGET_LINK_LIBRARIES(cmrep_vskel_api cmrep cmrep_dijkstra ${CMREP_FIT_LIBS} ${QHULL_LIBRARY})
SET(CMREP_API_LIBRARIES ${CMREP_API_LIBRARIES} cmrep_vskel_api)

ADD_EXECUTABLE(cmrep_vskel src/thinexe/cmrep_vskel.cxx)
TARGET_LINK_LIBRARIES(cmrep_vskel cmrep_vskel_api)
SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} cmrep_vskel)
ENDIF()


IF(CMREP_BUILD_BCMREP)
SET(CMREP_BCMREP_LIBS ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${IPOPT_LIBS} ${LAPACK_LIBRARIES} cmrep_toms611)
ADD_EXECUTABLE(bcmrep_test
  src/GentleNLP.cxx
  src/IPOptProblemInterface.cxx
  src/IPOptTest.cxx)
TARGET_LINK_LIBRARIES(bcmrep_test cmrep ${CMREP_BCMREP_LIBS})

ADD_EXECUTABLE(bcmrep_main
  src/GentleNLP.cxx
  src/IPOptProblemInterface.cxx
  src/IPOptQCQPProblemInterface.cxx
  src/BoundaryConstrainedCMRepObjectives.cxx
  src/BoundaryConstrainedCMRep.cxx
  src/VTKMeshBuilder.cxx
  src/QCQPProblem.h)
TARGET_LINK_LIBRARIES(bcmrep_main cmrep ${CMREP_BCMREP_LIBS})

SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} bcmrep_main bcmrep_test)
ENDIF()

# Include the utilities
ADD_EXECUTABLE(sphere_splat
  src/util/SplatOntoSphere.cxx)

ADD_EXECUTABLE(vtkprocrustes
  src/util/MeshProcrustes.cxx
  src/util/ReadWriteVTK.cxx
  src/Procrustes.cxx)

ADD_EXECUTABLE(vtkcontour
  src/util/RealImageToContour.cxx
  src/util/ReadWriteVTK.cxx)

ADD_EXECUTABLE(mesh2img
  src/util/FillSTLMesh.cxx
  src/util/DrawTriangles.cpp
  src/util/ReadWriteVTK.cxx)

ADD_EXECUTABLE(mesh_image_sample
  src/util/MeshImageSample.cxx)

ADD_EXECUTABLE(mesh_tetra_sample
  src/util/MeshTetraSample.cxx)

ADD_EXECUTABLE(mesh_merge_arrays
  src/util/MeshMergeArrays.cxx
  src/util/ReadWriteVTK.cxx)

ADD_EXECUTABLE(mesh_area_element
  src/util/MeshAreaElement.cxx)

ADD_EXECUTABLE(meshglm
  src/util/SurfaceClusterGLM.cxx
  src/util/asa243.c
  src/util/ReadWriteVTK.cxx
  src/Registry.cxx)

ADD_EXECUTABLE(meshdiff
  src/util/CompareMeshes.cxx
  src/util/DrawFillInside.cpp
  src/util/DrawTriangles.cpp)

ADD_EXECUTABLE(warpmesh
  src/util/WarpMeshBackwards.cxx)

# This is a one-off tool for sampling from skeleton tetra meshes
ADD_EXECUTABLE(skel_tetra_max_thick_sphere
  src/util/SkelTetraMaxThicknessSphere.cxx)

# This uses triangle library to triangulate a 2d contour
ADD_EXECUTABLE(contour2surf
  src/util/ContourMesh2SurfMesh.cxx
  src/util/ReadWriteVTK.cxx)

ADD_EXECUTABLE(apply_remeshing_to_another
  src/util/ApplyRemeshingToAnotherMesh.cxx
  src/util/ReadWriteVTK.cxx)

ADD_EXECUTABLE(mesh_convert_format
  src/util/MeshConvertFormat.cxx
  src/util/ReadWriteVTK.cxx)

ADD_EXECUTABLE(ml_affine
  src/util/MultiLabelAffine.cxx
  src/GentleNLP.cxx)

# List of cmrep utilities
SET(CMREP_UTIL_TOOLS
  sphere_splat vtkprocrustes vtkcontour
  mesh2img mesh_image_sample mesh_tetra_sample mesh_merge_arrays mesh_area_element
  meshglm meshdiff warpmesh skel_tetra_max_thick_sphere
  contour2surf apply_remeshing_to_another mesh_convert_format ml_affine)

# Tetrahedron tools
ADD_EXECUTABLE(tetdumpattr src/util/TetDumpAttr.cxx)
ADD_EXECUTABLE(tetfill src/util/TetraFill.cxx)
ADD_EXECUTABLE(tetjac src/util/TetJacobian.cxx)
ADD_EXECUTABLE(tetmergearr src/util/MergeVTKArrays.cxx)
ADD_EXECUTABLE(tetsample src/util/TetSampleImage.cxx)
SET(CMREP_UTIL_TOOLS ${CMREP_UTIL_TOOLS} tetdumpattr tetfill tetjac tetmergearr tetsample)

IF(CMREP_BUILD_VCG_UTILS)
SET(VCG_UTIL_SRC
  src/util/VCGTriMesh.h
  src/util/VCGTriMesh.cxx
  ${VCGLIB_DIR}/wrap/ply/plylib.cpp)

ADD_EXECUTABLE(mesh_smooth_curv
  src/util/MeshSmoothAndCurvature.cxx
  ${VCG_UTIL_SRC} )
TARGET_LINK_LIBRARIES(mesh_smooth_curv ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(mesh_poisson_sample
  src/util/MeshPoissonSample.cxx
  ${VCG_UTIL_SRC})
TARGET_LINK_LIBRARIES(mesh_poisson_sample ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(mesh_decimate_vcg
  src/util/MeshDecimateVCG.cxx
  src/util/ReadWriteVTK.cxx
  ${VCG_UTIL_SRC})
TARGET_LINK_LIBRARIES(mesh_decimate_vcg ${LAPACK_LIBRARIES})

ADD_EXECUTABLE(vtklevelset
  src/util/RealImageToMesh.cxx
  src/util/ReadWriteVTK.cxx
  ${VCG_UTIL_SRC})

SET(CMREP_UTIL_TOOLS ${CMREP_UTIL_TOOLS} mesh_smooth_curv mesh_poisson_sample mesh_decimate_vcg vtklevelset)
TARGET_LINK_LIBRARIES(vtklevelset cmrep ${CMREP_FIT_LIBS} ${LAPACK_LIBRARIES})
ELSE()
# vtklevelset can be built without VCG support
ADD_EXECUTABLE(vtklevelset
  src/util/RealImageToMesh.cxx
  src/util/ReadWriteVTK.cxx)
TARGET_LINK_LIBRARIES(vtklevelset cmrep ${CMREP_FIT_LIBS})
SET(CMREP_UTIL_TOOLS ${CMREP_UTIL_TOOLS} vtklevelset)
ENDIF()


IF(CMREP_BUILD_TETGEN_UTILS)
ADD_EXECUTABLE(tetgen2vtk
  src/util/TetGenToVTK.cxx
  src/util/ReadWriteVTK.cxx
  src/util/TetraLib.cxx)
ADD_EXECUTABLE(surf2tet
  src/util/SurfMesh2VolMesh.cxx
  src/util/ReadWriteVTK.cxx
  src/util/TetraLib.cxx)
TARGET_LINK_LIBRARIES(tetgen2vtk ${CMREP_UTIL_LIBS} ${TETGEN_LIBRARY})
TARGET_LINK_LIBRARIES(surf2tet ${CMREP_UTIL_LIBS} ${TETGEN_LIBRARY})
SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} tetgen2vtk surf2tet)
ENDIF()

# Link libraries for the tools
FOREACH(CMREP_TOOL ${CMREP_UTIL_TOOLS})
vtk_module_autoinit(TARGETS ${CMREP_TOOL} MODULES ${VTK_LIBRARIES})
TARGET_LINK_LIBRARIES(${CMREP_TOOL} ${CMREP_UTIL_LIBS})
ENDFOREACH(CMREP_TOOL)

# Add to the list of binaries
SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} ${CMREP_UTIL_TOOLS})

# Install all the targets
FOREACH(CMREP_TOOL ${CMREP_BINARY_TOOLS})
INSTALL(TARGETS ${CMREP_TOOL}
  BUNDLE DESTINATION . COMPONENT Runtime
  RUNTIME DESTINATION bin COMPONENT RUNTIME)
ENDFOREACH()

FOREACH(CMREP_LIB ${CMREP_API_LIBRARIES})
INSTALL(TARGETS ${CMREP_LIB} DESTINATION lib EXPORT CMRepTargets)
ENDFOREACH()

# Create files that enable FindPackage to work
EXPORT(EXPORT CMRepTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/CMRepTargets.cmake"
  NAMESPACE CMRep::)

INSTALL(EXPORT CMRepTargets
  FILE CMRepTargets.cmake
  NAMESPACE CMRep::
  DESTINATION ${CMREP_BINARY_DIR})

include(CMakePackageConfigHelpers)

install(EXPORT CMRepTargets
  FILE CMRepTargets.cmake
  DESTINATION lib/cmake/CMRep)

# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/CMRepConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/CMRep"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/CMRepConfigVersion.cmake"
  VERSION "${CMREP_VERSION_FULL}"
  COMPATIBILITY AnyNewerVersion)

# install the configuration file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CMRepConfig.cmake
  DESTINATION lib/cmake/CMRep)

INSTALL(DIRECTORY ${CMREP_SOURCE_DIR}/examples DESTINATION "share/cmrep/")


OPTION(CMREP_FIXUP_BUNDLE "Whether to run FixupBundle code at install/package time" ON)

# INCLUDE(BundleUtilities)

# Fixup bundle on MACOS
IF(CMREP_FIXUP_BUNDLE)
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
ENDIF()


# Add package support
INCLUDE(CMake/Package.cmake)



