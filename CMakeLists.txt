PROJECT(CMREP)

# Get ITK
FIND_PACKAGE(ITK)
  IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)
ENDIF(USE_ITK)

# Get VTK
FIND_PACKAGE(VTK)
  IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
  ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR
      "Cannot build without VTK.  Please set VTK_DIR.")
  ENDIF(VTK_FOUND)
ENDIF(USE_VTK)

# Set the include directories
INCLUDE_DIRECTORIES(
  ${CMREP_SOURCE_DIR}/src)

# Include subdirectories
SUBDIRS(extras/toms611)

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
  src/MedialAtom.cxx
  src/MedialAtomGrid.cxx
  src/MeshMedialPDESolver.cxx
  src/MedialModelIO.cxx
  src/MedialPDEMasks.cxx
  src/MedialPDESites.cxx
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
