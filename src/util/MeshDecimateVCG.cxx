#include "VCGTriMesh.h"
#include "ReadWriteVTK.h"

int usage()
{
  printf("mesh_decimate_vcg: reduce the size of a surface mesh using VCG library\n");
  printf("Usage: \n");
  printf("  mesh_decimate_vcg <input_mesh> <output_mesh> <reduction_factor>\n");
  printf("Required parameters:\n");
  printf("  reduction_factor    : If >1, desired number of faces; if <1, reduction factor\n");
  return -1;
}

int main(int argc,char ** argv)
{
  if(argc < 4)
  {
    return usage();
  }

  vtkSmartPointer<vtkPolyData> mesh = ReadVTKData(argv[1]);
  vtkSmartPointer<vtkPolyData> mesh_reduced = vtkPolyData::New();
  double reduction_factor = atof(argv[3]);

  VCGTriMesh tri_mesh;
  tri_mesh.ImportFromVTK(mesh);
  tri_mesh.CleanMesh();

  vcg::tri::TriEdgeCollapseQuadricParameter qparams;
  qparams.QualityThr = 0.3;
  qparams.PreserveBoundary=true;
  qparams.BoundaryQuadricWeight=1.0;
  qparams.PreserveTopology=true;
  qparams.QualityWeight=false;
  qparams.NormalCheck=true;
  qparams.OptimalPlacement=true;
  qparams.QualityQuadric=true;
  qparams.QualityQuadricWeight=0.001;

  // Perform decimation
  tri_mesh.QuadricEdgeCollapseRemeshing(reduction_factor, qparams);

  // Postprocess by cleaning and fixing normals
  tri_mesh.CleanMesh();
  tri_mesh.RecomputeNormals();
  tri_mesh.ExportToVTK(mesh_reduced);

  // Write the mesh
  WriteVTKData(mesh_reduced, argv[2]);
  return 0;
}
