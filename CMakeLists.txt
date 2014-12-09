# Requisite CMake commands
cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Project name
PROJECT(CMREP)

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
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# Where to get additional modules
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")

# Optionally, support PARDISO
OPTION(USE_PARDISO "Use the PARDISO sparse solver library" OFF)
IF(USE_PARDISO)
  IF(WIN32)
    # On Windows, the library can only be loaded dynamically
    ADD_DEFINITIONS(-DHAVE_PARDISO -DPARDISO_DYNLOAD)  
  ELSE(WIN32)

    # Search for PARDISO
    FIND_LIBRARY(PARDISO_LIB pardiso DOC "The location of the PARDISO shared library")
    IF(PARDISO_LIB)
      ADD_DEFINITIONS(-DHAVE_PARDISO)
      SET(SOLVER_SRC src/PardisoInterface.cxx)
    ENDIF(PARDISO_LIB)  

    # On Linux, Lapack is also required
    INCLUDE(${CMREP_SOURCE_DIR}/CMake/FindLapack.cmake)

    # Set the PDE libraries variable
    IF(PARDISO_LIB AND LAPACK_FOUND)
      SET(PDE_LIBS ${PARDISO_LIB} ${LAPACK_LIBS})
    ENDIF(PARDISO_LIB AND LAPACK_FOUND)

  ENDIF(WIN32)
ELSE(USE_PARDISO)
  ADD_DEFINITIONS(-UHAVE_PARDISO)
ENDIF(USE_PARDISO)

# Optionally, support TAUCS
OPTION(USE_TAUCS "Use the TAUCS sparse solver library" OFF)
IF(USE_TAUCS)
  FIND_PATH(TAUCS_BUILD_DIR "taucs_config_build.h")
  IF(TAUCS_BUILD_DIR)
    GET_FILENAME_COMPONENT(TAUCS_BUILD_DIR_PARENT ${TAUCS_BUILD_DIR} PATH)
    GET_FILENAME_COMPONENT(TAUCS_ROOT_DIR ${TAUCS_BUILD_DIR_PARENT} PATH)
    GET_FILENAME_COMPONENT(TAUCS_SYS ${TAUCS_BUILD_DIR} NAME_WE)
    INCLUDE_DIRECTORIES(${TAUCS_BUILD_DIR} ${TAUCS_ROOT_DIR}/src)
    FIND_LIBRARY(TAUCS_LIB taucs PATHS ${TAUCS_ROOT_DIR}/lib/${TAUCS_SYS} DOC "The location of the TAUCS library")
    IF(TAUCS_LIB)
      ADD_DEFINITIONS(-DHAVE_TAUCS)
      SET(SOLVER_SRC src/TaucsInterface.cxx)
    ENDIF(TAUCS_LIB)  

    # On Linux, Lapack is also required
    INCLUDE(${CMREP_SOURCE_DIR}/CMake/FindLapack.cmake)
    FIND_LIBRARY(METIS_LIB metis)

    # Set the PDE libraries variable
    IF(TAUCS_LIB AND LAPACK_FOUND AND METIS_LIB)
      SET(PDE_LIBS ${TAUCS_LIB} ${LAPACK_LIBS} ${METIS_LIB})
    ENDIF(TAUCS_LIB AND LAPACK_FOUND AND METIS_LIB)

  ENDIF(TAUCS_BUILD_DIR)
ENDIF(USE_TAUCS)

# Optionally, use MKL
OPTION(USE_MKL "Use the Intel MKL library" OFF)
IF(USE_MKL)
  FIND_PACKAGE(MKL REQUIRED)
  IF(MKL_LIBRARIES)

    ADD_DEFINITIONS(-DHAVE_MKL)
    INCLUDE_DIRECTORIES(${MKL_INCLUDES})

    SET(SOLVER_SRC src/MKLSolverInterface.cxx)
    SET(PDE_LIBS ${MKL_LIBRARIES})

  ENDIF(MKL_LIBRARIES)
ENDIF(USE_MKL)

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
  LINK_LIBRARIES(${IPOPT_LIBRARIES})
ENDIF(USE_IPOPT)

OPTION(USE_TETGEN "Use the TetGen library" OFF)
IF(USE_TETGEN)
  FIND_PACKAGE(TetGen REQUIRED)
ENDIF(USE_TETGEN)

# Set the include directories
INCLUDE_DIRECTORIES(${CMREP_SOURCE_DIR}/src ${CMREP_SOURCE_DIR}/src/dijkstra)

# Find BOOST
IF(USE_CGAL)
  FIND_PACKAGE(BOOST)
  IF(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
  ENDIF(Boost_FOUND)

  # Find CGAL
  FIND_PATH(CGAL_ROOT_DIR NAMES include/CGAL/Surface_mesh_default_triangulation_3.h DOC CGAL main directory)
  IF(CGAL_ROOT_DIR)

    FIND_PATH(CGAL_INCLUDE_DIR NAMES CGAL/Surface_mesh_default_triangulation_3.h PATHS ${CGAL_ROOT_DIR}/include)
    FIND_PATH(CGAL_GMP_INCLUDE_DIR NAMES gmp.h PATHS ${CGAL_ROOT_DIR}/auxiliary/gmp/include)

    INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR})
    INCLUDE_DIRECTORIES(${CGAL_GMP_INCLUDE_DIR})
    LINK_DIRECTORIES(${CGAL_ROOT_DIR}/lib)
    LINK_DIRECTORIES(${CGAL_ROOT_DIR}/auxiliary/gmp/lib)
    
    IF(MSVC)
      FIND_PATH(CGAL_SYSTEM_INCLUDE_DIR CGAL/compiler_config.h PATHS ${CGAL_INCLUDE_DIR}/CGAL/config/msvc)
      INCLUDE_DIRECTORIES(${CGAL_SYSTEM_INCLUDE_DIR})
    ENDIF(MSVC)
    
  ENDIF(CGAL_ROOT_DIR)
ENDIF(USE_CGAL)

# Get GLUT
FIND_PACKAGE(GLUT)
IF(GLUT_FOUND)
  INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
  LINK_DIRECTORIES(${GLUT_INCLUDE_DIR})
ENDIF(GLUT_FOUND)

#  ${CMREP_SOURCE_DIR}/extras/afront/tmp/

# Include subdirectories
SUBDIRS(extras/toms611)
SUBDIRS(testing)


OPTION(BUILD_UTILS "Build some cm-rep utilities" ON)

IF(BUILD_UTILS)
  SUBDIRS(src/util)
ENDIF(BUILD_UTILS)

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

# Source code for the Dijkstra library
SET(DIJKSTRA_SRCS
  src/dijkstra/ShortestPath.cxx
  src/dijkstra/VTKMeshShortestDistance.cxx)

# Create library of shortest path routines
ADD_LIBRARY(cmrep_dijkstra ${DIJKSTRA_SRCS})

# Define the VTK libraries
SET(VTK_LIBRARIES
  vtkIO vtkHybrid vtkRendering vtkImaging vtkGraphics vtkFiltering vtkCommon)


# Set the list of required libraries
LINK_LIBRARIES(
  cmrep
  toms611
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES})

# Link the main executables
ADD_EXECUTABLE(cmrep_test src/ScriptTest.cxx)
TARGET_LINK_LIBRARIES(cmrep_test ${PDE_LIBS})

ADD_EXECUTABLE(cmrep_fit src/Fit.cxx)
TARGET_LINK_LIBRARIES(cmrep_fit ${PDE_LIBS})

ADD_EXECUTABLE(cmrep_getmesh src/CMRepToMesh.cxx)
TARGET_LINK_LIBRARIES(cmrep_getmesh ${PDE_LIBS})

ADD_EXECUTABLE(cmrep_subdivide src/CMRepSubdivide.cxx)
TARGET_LINK_LIBRARIES(cmrep_subdivide ${PDE_LIBS})

ADD_EXECUTABLE(cmrep_remesh src/CMRepRemesh.cxx)
TARGET_LINK_LIBRARIES(cmrep_remesh ${PDE_LIBS})

ADD_EXECUTABLE(cmrep_meanmodel src/CMRepMeanModel.cxx)
TARGET_LINK_LIBRARIES(cmrep_meanmodel ${PDE_LIBS})

ADD_EXECUTABLE(cmrep_pca_generate src/ComputePCA.cxx)
TARGET_LINK_LIBRARIES(cmrep_pca_generate ${PDE_LIBS})

ADD_EXECUTABLE(cmrep_fillmesh src/CMRepMeshToImage.cxx)
ADD_EXECUTABLE(cmrep_genwarp src/CMRepToWarp.cxx)
ADD_EXECUTABLE(cmrep_warpimage src/CMRepWarpImage.cxx)
ADD_EXECUTABLE(cmrep_afftran src/CMRepToAffineTransform.cxx)

ADD_EXECUTABLE(cmrep_vskel src/VoronoiSkeletonTool.cxx)
TARGET_LINK_LIBRARIES(cmrep_vskel cmrep_dijkstra)

IF(USE_IPOPT)
  ADD_EXECUTABLE(opttest src/GentleNLP.cxx src/IPOptProblemInterface.cxx src/IPOptTest.cxx)
  ADD_EXECUTABLE(contest
                 src/GentleNLP.cxx
                 src/IPOptProblemInterface.cxx
                 src/ConstrainedCMRepObjectives.cxx
                 src/ConstrainedCMRepTest.cxx)
  TARGET_LINK_LIBRARIES(opttest ${PDE_LIBS})
  TARGET_LINK_LIBRARIES(contest ${PDE_LIBS})

  INSTALL(TARGETS contest 
     BUNDLE DESTINATION . COMPONENT Runtime
     RUNTIME DESTINATION bin COMPONENT RUNTIME)

ENDIF(USE_IPOPT)

# Add installation targets
INSTALL(TARGETS cmrep_fit cmrep_test cmrep_afftran cmrep_getmesh
  cmrep_meanmodel cmrep_subdivide cmrep_remesh cmrep_fillmesh cmrep_genwarp
  cmrep_vskel cmrep_warpimage 
     BUNDLE DESTINATION . COMPONENT Runtime
     RUNTIME DESTINATION bin COMPONENT RUNTIME)

INSTALL(DIRECTORY ${CMREP_SOURCE_DIR}/examples DESTINATION
  "share/cmrep/")

INCLUDE(BundleUtilities)

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

# Add package support
INCLUDE(CMake/Package.cmake)
