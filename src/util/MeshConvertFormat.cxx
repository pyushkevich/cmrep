#include <cstdio>
#include <iostream>
#include <vtkSmartPointer.h>
#include "ReadWriteVTK.h"

using namespace std;

int usage()
{
  cout << "mesh_convert_format: format conversion for polydata formats\n"
       << "usage:\n"
       << "  mesh_convert_format [options] <input> <output>\n"
       << "options:\n";
  return -1;
}

int main(int argc, char *argv[])
{
  typedef vtkSmartPointer<vtkPolyData> PolyDataPtr;

  if(argc < 2)
    return usage();

  // Read the three meshes
  PolyDataPtr pd_src = ReadVTKData(argv[argc-2]);

  // Save output mesh
  WriteVTKData(pd_src, argv[argc-1]);

  return 0;
}
