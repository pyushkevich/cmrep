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

# Keep track of all binaries that will be compiled
SET(CMREP_BINARY_TOOLS "")

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

    FIND_LIBRARY(GOMP_LIB gomp DOC "The location of the libgomp library")

    # Set the PDE libraries variable
    IF(PARDISO_LIB AND LAPACK_FOUND)
      SET(PDE_LIBS ${PARDISO_LIB} ${LAPACK_LIBS} ${GOMP_LIB})
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
  MESSAGE(STATUS "Using IPOPT libraries ${IPOPT_LIBRARIES}")
  SET(PDE_LIBS ${PDE_LIBS} ${IPOPT_LIBRARIES})
ENDIF(USE_IPOPT)

OPTION(USE_TETGEN "Use the TetGen library" OFF)
IF(USE_TETGEN)
  FIND_PACKAGE(TetGen REQUIRED)
ENDIF(USE_TETGEN)

# Set the include directories
INCLUDE_DIRECTORIES(${CMREP_SOURCE_DIR}/src ${CMREP_SOURCE_DIR}/src/dijkstra)

# Not sure why we need OpenGL
FIND_PACKAGE(OpenGL)

# Find NLOPT - a small non-linear optimization package
OPTION(USE_NLOPT "Use the NLOPT optimization library" OFF)
IF(USE_NLOPT)
  FIND_PACKAGE(NLOPT REQUIRED)
  INCLUDE_DIRECTORIES(${NLOPT_INCLUDE_DIRS})
  GET_FILENAME_COMPONENT(NLOPT_LIBRARY_DIR ${NLOPT_LIBRARIES} DIRECTORY)
ENDIF(USE_NLOPT)

# Include subdirectories
SUBDIRS(extras/toms611)
SUBDIRS(testing)

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

# Source code for the Dijkstra library
SET(DIJKSTRA_SRCS
  src/dijkstra/ShortestPath.cxx
  src/dijkstra/VTKMeshShortestDistance.cxx)

# Create library of shortest path routines
ADD_LIBRARY(cmrep_dijkstra ${DIJKSTRA_SRCS})

# Set the list of required libraries
LINK_LIBRARIES(
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${SYSTEM_LIBS})

# New: Geodesic shooting
IF(LAPACK_FOUND OR IPOPT_LAPACK_FOUND)

  ADD_EXECUTABLE(shooting_test src/PointSetHamiltonianSystem.cxx src/TestGeodesicShooting.cxx)
  TARGET_LINK_LIBRARIES(shooting_test cmrep ${PDE_LIBS})

  ADD_EXECUTABLE(lmshoot 
    src/PointSetHamiltonianSystem.cxx 
    src/PointSetGeodesicShooting.cxx
    src/util/ReadWriteVTK.cxx)
  TARGET_LINK_LIBRARIES(lmshoot cmrep ${PDE_LIBS})

  ADD_EXECUTABLE(lmflow
    src/PointSetHamiltonianSystem.cxx
    src/util/ReadWriteVTK.cxx
    src/PointSetGeodesicFlowUtil.cxx)
  TARGET_LINK_LIBRARIES(lmflow cmrep ${PDE_LIBS})

  INCLUDE_DIRECTORIES(submodules/greedy/src/ITKFilters/include submodules/greedy/src)
  ADD_EXECUTABLE(lmtowarp
    src/PointSetHamiltonianSystem.cxx
    src/util/ReadWriteVTK.cxx
    src/PointSetGeodesicToWarp.cxx
    submodules/greedy/src/lddmm_data.cxx)
  TARGET_LINK_LIBRARIES(lmtowarp cmrep ${PDE_LIBS})

  # Update the list of binaries
  SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} shooting_test lmshoot lmtowarp)

  IF(NLOPT_FOUND)
    ADD_EXECUTABLE(alcmrep_test
      src/TestAugLagCMRep.cxx
      src/MeshTraversal.cxx
      src/PointSetOptimalControlSystem.cxx
      src/PointSetHamiltonianSystem.cxx
      src/VTKMeshBuilder.cxx)
    TARGET_LINK_LIBRARIES(alcmrep_test cmrep ${PDE_LIBS} ${NLOPT_LIBRARIES})
    SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} alcmrep_test)
  ENDIF(NLOPT_FOUND)

ENDIF(LAPACK_FOUND OR IPOPT_LAPACK_FOUND)

# Link the main executables
ADD_EXECUTABLE(cmrep_test src/ScriptTest.cxx)
TARGET_LINK_LIBRARIES(cmrep_test cmrep ${PDE_LIBS} toms611)

ADD_EXECUTABLE(cmrep_fit src/Fit.cxx)
TARGET_LINK_LIBRARIES(cmrep_fit cmrep cmrep ${PDE_LIBS} toms611)

ADD_EXECUTABLE(cmrep_getmesh src/CMRepToMesh.cxx)
TARGET_LINK_LIBRARIES(cmrep_getmesh cmrep ${PDE_LIBS} toms611)

ADD_EXECUTABLE(cmrep_subdivide src/CMRepSubdivide.cxx)
TARGET_LINK_LIBRARIES(cmrep_subdivide cmrep ${PDE_LIBS} toms611)

ADD_EXECUTABLE(cmrep_remesh src/CMRepRemesh.cxx)
TARGET_LINK_LIBRARIES(cmrep_remesh cmrep ${PDE_LIBS} toms611)

ADD_EXECUTABLE(cmrep_meanmodel src/CMRepMeanModel.cxx)
TARGET_LINK_LIBRARIES(cmrep_meanmodel cmrep ${PDE_LIBS} toms611)

ADD_EXECUTABLE(cmrep_pca_generate src/ComputePCA.cxx)
TARGET_LINK_LIBRARIES(cmrep_pca_generate cmrep ${PDE_LIBS} toms611)

ADD_EXECUTABLE(cmrep_fillmesh src/CMRepMeshToImage.cxx)
TARGET_LINK_LIBRARIES(cmrep_fillmesh cmrep)

ADD_EXECUTABLE(cmrep_genwarp src/CMRepToWarp.cxx)
TARGET_LINK_LIBRARIES(cmrep_genwarp cmrep)

ADD_EXECUTABLE(cmrep_warpimage src/CMRepWarpImage.cxx)
TARGET_LINK_LIBRARIES(cmrep_warpimage cmrep)

ADD_EXECUTABLE(cmrep_afftran src/CMRepToAffineTransform.cxx)
TARGET_LINK_LIBRARIES(cmrep_afftran cmrep)

ADD_EXECUTABLE(cmrep_vskel src/VoronoiSkeletonTool.cxx)
TARGET_LINK_LIBRARIES(cmrep_vskel cmrep cmrep_dijkstra)

# Update the list of binaries
SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} 
  cmrep_test cmrep_fit cmrep_getmesh cmrep_remesh cmrep_meanmodel 
  cmrep_pca_generate cmrep_fillmesh cmrep_genwarp cmrep_warpimage
  cmrep_afftran cmrep_vskel)

IF(USE_IPOPT AND USE_TETGEN)
  ADD_EXECUTABLE(opttest src/GentleNLP.cxx src/IPOptProblemInterface.cxx src/IPOptTest.cxx)
  ADD_EXECUTABLE(contest
                 src/GentleNLP.cxx
                 src/IPOptProblemInterface.cxx
                 src/ConstrainedCMRepObjectives.cxx
                 src/ConstrainedCMRepTest.cxx
                 src/VTKMeshBuilder.cxx)
  TARGET_LINK_LIBRARIES(opttest cmrep ${PDE_LIBS})
  TARGET_LINK_LIBRARIES(contest cmrep ${PDE_LIBS} toms611)

  SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} opttest contest)

ENDIF(USE_IPOPT AND USE_TETGEN)

# Include the utilities
IF(BUILD_UTILS)

  ADD_EXECUTABLE(sphere_splat src/util/SplatOntoSphere.cxx)
  ADD_EXECUTABLE(vtkprocrustes 
    src/util/MeshProcrustes.cxx 
    src/util/ReadWriteVTK.cxx 
    src/Procrustes.cxx)

  ADD_EXECUTABLE(vtklevelset src/util/RealImageToMesh.cxx src/util/ReadWriteVTK.cxx)
  TARGET_LINK_LIBRARIES(vtklevelset cmrep ${PDE_LIBS} toms611)

  ADD_EXECUTABLE(vtkcontour src/util/RealImageToContour.cxx src/util/ReadWriteVTK.cxx)

  ADD_EXECUTABLE(mesh2img src/util/FillSTLMesh.cxx src/util/DrawTriangles.cpp)

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

  # Add to the list of binaries
  SET(CMREP_BINARY_TOOLS ${CMREP_BINARY_TOOLS} 
    sphere_splat vtkprocrustes vtklevelset vtkcontour
    mesh2img mesh_image_sample mesh_merge_arrays mesh_area_element
    meshglm meshcluster meshdiff warpmesh)

ENDIF(BUILD_UTILS)

# Install all the targets
FOREACH(CMREP_TOOL ${CMREP_BINARY_TOOLS})

  INSTALL(TARGETS ${CMREP_TOOL} 
    BUNDLE DESTINATION . COMPONENT Runtime
    RUNTIME DESTINATION bin COMPONENT RUNTIME)

ENDFOREACH(CMREP_TOOL)

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

# Fixup bundle on LINUX
IF(UNIX AND NOT APPLE)

  SET(APPS "\${CMAKE_INSTALL_PREFIX}/bin/cmrep_fit")
  SET(FIXUP_DIRS "${IPOPT_LIB_DIR}" "${VTK_DIR}/lib" "${NLOPT_LIBRARY_DIR}")

  INSTALL(CODE "
     include(BundleUtilities)
     fixup_bundle(\"${APPS}\"   \"\"   \"${FIXUP_DIRS}\")
     " COMPONENT Runtime)

ENDIF(UNIX AND NOT APPLE)


# Add package support
INCLUDE(CMake/Package.cmake)


#--------------------------------------------------------------------------------
# Uploading code to SourceForge
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Configure SCP

FIND_PROGRAM(SCP_PROGRAM NAMES scp DOC "Location of the scp program (optional)")
MARK_AS_ADVANCED(SCP_PROGRAM)

SET(SCP_ARGUMENTS "-v" CACHE STRING "Optional arguments to the scp command for uploads to SourceForge")
MARK_AS_ADVANCED(SCP_ARGUMENTS)

SET(SCP_USERNAME "" CACHE STRING "SourceForge.net account id for uploads")
MARK_AS_ADVANCED(SCP_USERNAME)

SET(NIGHTLY_TARGET "cmrep-nightly-${CPACK_SYSTEM_NAME}.${CPACK_EXTENSION}")
SET(EXPERIMENTAL_TARGET "cmrep-experimental-${CPACK_SYSTEM_NAME}.${CPACK_EXTENSION}")

SET(SCP_ROOT "frs.sourceforge.net:/home/frs/project/c/cm/cmrep/cmrep")

#--------------------------------------------------------------------------------
# Create targets

SET(CPACK_PACKAGE_FILE_NAME_WEXT "${CPACK_PACKAGE_FILE_NAME}.${CPACK_EXTENSION}")

ADD_CUSTOM_TARGET(cmrep_upload_nightly
  VERBATIM COMMAND "${SCP_PROGRAM}" ${SCP_ARGUMENTS}
  ${CPACK_PACKAGE_FILE_NAME_WEXT} ${SCP_USERNAME},itk-snap@${SCP_ROOT}/Nightly/${NIGHTLY_TARGET}
  DEPENDS ${CPACK_TARGET}
  WORKING_DIRECTORY ${CMREP_BINARY_DIR}
  COMMENT "Uploading package ${CPACK_PACKAGE_FILE_NAME_WEXT} to SourceForge.net as ${NIGHTLY_TARGET}")

ADD_CUSTOM_TARGET(cmrep_upload_experimental
  VERBATIM COMMAND "${SCP_PROGRAM}" ${SCP_ARGUMENTS} 
    ${CPACK_PACKAGE_FILE_NAME_WEXT} ${SCP_USERNAME},itk-snap@${SCP_ROOT}/Experimental/${EXPERIMENTAL_TARGET}
  DEPENDS ${CPACK_TARGET}
  WORKING_DIRECTORY ${CMREP_BINARY_DIR}
  COMMENT "Uploading package ${CPACK_PACKAGE_FILE_NAME_WEXT} to SourceForge.net as ${EXPERIMENTAL_TARGET}")



