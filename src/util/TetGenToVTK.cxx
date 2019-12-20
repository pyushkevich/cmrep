#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>
#include "ReadWriteVTK.h"
#include "TetraLib.h"

using namespace std;

typedef vnl_matrix_fixed<double, 4, 4> MatrixType;

int usage()
{
  cout << "tetgen2vtk - Converts a TETGEN tetra mesh to a VTK unstructured grid" << endl;
  cout << "usage: " << endl;
  cout << "   tetgen2vtk tetinput output.vtk" << endl;
  cout << "params: " << endl;
  cout << "   tetinput      The tetgen output base name without extension" << endl;
  return -1;
}

int main(int argc, char **argv)
{
  tetgenio in;

  // Check the parameters
  if(argc < 3) return usage();

  // Read the mesh
  in.load_tetmesh(argv[1]);

  // Write to a mesh
  WriteTetgenOutputAsUnstructuredMesh(in, argv[2]);

  return 0;
}
