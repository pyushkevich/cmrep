PROJECT(CMREP)

# New CMake compatibility
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


# Get ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


# Get VTK
FIND_PACKAGE(VTK)
IF(VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
  SET(VTK_LIBRARIES
    vtkCommon vtkIO vtkGraphics vtkRendering vtkFiltering vtkImaging vtkHybrid)
  IF(VTK_MAJOR_VERSION EQUAL 4)
    SET(VTK_LIBRARIES ${VTK_LIBRARIES} vtkPatented)
  ENDIF(VTK_MAJOR_VERSION EQUAL 4)
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without VTK.  Please set VTK_DIR.")
ENDIF(VTK_FOUND)

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
      LINK_LIBRARIES(${PARDISO_LIB})
    ENDIF(PARDISO_LIB)  

    # On Linux, Lapack is also required
    INCLUDE(${CMREP_SOURCE_DIR}/CMake/FindLapack.cmake)
    IF(LAPACK_FOUND)
      LINK_LIBRARIES(${LAPACK_LIBS})
    ENDIF(LAPACK_FOUND)

  ENDIF(WIN32)
ELSE(USE_PARDISO)
  ADD_DEFINITIONS(-UHAVE_PARDISO)
ENDIF(USE_PARDISO)


# Set the include directories
INCLUDE_DIRECTORIES(${CMREP_SOURCE_DIR}/src)

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
### SUBDIRS(extras/afront/tmp/afront)
SUBDIRS(testing)

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
  src/PardisoInterface.cxx
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
  src/TestSolver.cxx
  src/VTKMedialMesh.cxx
  src/VTKSubdivision.cxx)

# Create the CM-REP library
ADD_LIBRARY(cmrep ${COMMON_SRCS})

# Define the ITK libraries
SET(ITK_LIBRARIES
  ITKNumerics ITKIO ITKCommon ITKBasicFilters)


# Set the list of required libraries
LINK_LIBRARIES(
  cmrep
  toms611
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES})

# Link the main executables
ADD_EXECUTABLE(cmrep_test src/ScriptTest.cxx)
ADD_EXECUTABLE(cmrep_fit src/Fit.cxx)
ADD_EXECUTABLE(cmrep_getmesh src/CMRepToMesh.cxx)
ADD_EXECUTABLE(cmrep_subdivide src/CMRepSubdivide.cxx)
ADD_EXECUTABLE(cmrep_remesh src/CMRepRemesh.cxx)
ADD_EXECUTABLE(cmrep_fillmesh src/CMRepMeshToImage.cxx)
ADD_EXECUTABLE(cmrep_genwarp src/CMRepToWarp.cxx)
ADD_EXECUTABLE(cmrep_warpimage src/CMRepWarpImage.cxx)
ADD_EXECUTABLE(cmrep_afftran src/CMRepToAffineTransform.cxx)

# Add installation targets
INSTALL_TARGETS("/bin" cmrep_fit cmrep_test cmrep_afftran cmrep_getmesh cmrep_subdivide 
  cmrep_remesh cmrep_fillmesh cmrep_genwarp cmrep_warpimage)
