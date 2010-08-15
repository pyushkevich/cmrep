ADD_EXECUTABLE(vtklevelset RealImageToMesh.cxx ReadWriteVTK.cxx)
TARGET_LINK_LIBRARIES(vtklevelset ${ITK_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(mesh2img FillSTLMesh.cxx DrawTriangles.cpp)
TARGET_LINK_LIBRARIES(mesh2img ${ITK_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(meshcluster 
  SurfaceClusterAnalysis.cxx 
  ReadWriteVTK.cxx 
  ${CMREP_SOURCE_DIR}/src/Registry.cxx)
TARGET_LINK_LIBRARIES(meshcluster ${ITK_LIBRARIES} ${VTK_LIBRARIES})


# Add installation targets
INSTALL_TARGETS("/bin" vtklevelset mesh2img meshcluster)

