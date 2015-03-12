ADD_EXECUTABLE(vtklevelset RealImageToMesh.cxx ReadWriteVTK.cxx)

ADD_EXECUTABLE(vtkcontour RealImageToContour.cxx ReadWriteVTK.cxx)

ADD_EXECUTABLE(mesh2img FillSTLMesh.cxx DrawTriangles.cpp)

ADD_EXECUTABLE(mesh_image_sample MeshImageSample.cxx)

ADD_EXECUTABLE(mesh_merge_arrays MeshMergeArrays.cxx ReadWriteVTK.cxx)

ADD_EXECUTABLE(meshglm 
  SurfaceClusterGLM.cxx 
  asa243.c
  ReadWriteVTK.cxx 
  ${CMREP_SOURCE_DIR}/src/Registry.cxx)

ADD_EXECUTABLE(meshcluster 
  SurfaceClusterAnalysis.cxx 
  ReadWriteVTK.cxx 
  ${CMREP_SOURCE_DIR}/src/Registry.cxx)

# ADD_EXECUTABLE(meshdiff CompareMeshes.cxx DrawFillInside.cpp DrawTriangles.cpp)
# TARGET_LINK_LIBRARIES(meshdiff ${ITK_LIBRARIES} ${VTK_LIBRARIES})


# Add installation targets
INSTALL_TARGETS("/bin" vtklevelset mesh2img meshcluster meshglm mesh_image_sample vtkcontour mesh_merge_arrays)

