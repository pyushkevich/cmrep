PROJECT(CMREP)

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
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without VTK.  Please set VTK_DIR.")
ENDIF(VTK_FOUND)

# Find BOOST
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


# Get GLUT
FIND_PACKAGE(GLUT)
IF(GLUT_FOUND)
  INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
  LINK_DIRECTORIES(${GLUT_INCLUDE_DIR})
ENDIF(GLUT_FOUND)

# Set the include directories
INCLUDE_DIRECTORIES(
  ${CMREP_SOURCE_DIR}/src
  )

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

# Define the VTK libraries
SET(VTK_LIBRARIES
  vtkCommon vtkIO vtkGraphics vtkRendering vtkFiltering vtkImaging vtkHybrid)

# Set the list of required libraries
LINK_LIBRARIES(
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  toms611
  cmrep)

# Link the main executables
ADD_EXECUTABLE(testpde src/ScriptTest.cxx)
ADD_EXECUTABLE(cmr_fit src/Fit.cxx)
ADD_EXECUTABLE(cmrep2mesh src/CMRepToMesh.cxx)
ADD_EXECUTABLE(cmrep_subdivide src/CMRepSubdivide.cxx)
ADD_EXECUTABLE(cmrep_remesh src/CMRepRemesh.cxx)

# Add installation targets
INSTALL_TARGETS("/bin" cmr_fit testpde cmrep2mesh cmrep_subdivide cmrep_remesh)