#include <iostream>
#include <string>
#include <vector>
#include <vtksys/RegularExpression.hxx>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include "ReadWriteVTK.h"

using namespace std;
using namespace vtksys;

int usage()
{
  cout <<
    "vtkmergearr: merge arrays across many VTK meshes\n"
    "usage:\n"
    "  vtkmergearr ref.vtk out.vtk ArrayName [cell|point] mesh1.vtk mesh2.vtk ... meshN.vtk\n"
    "parameters:\n"
    "  ref.vtk  Mesh that will be used as reference. Must match all inputs\n"
    "  out.vtk  Output filename\n";
  return -1;
}

template <class TMeshType>
int MyMain(int argc, char *argv[])
{
  if(argc < 6)
    return usage();

  const char *fn_ref = argv[1];
  const char *fn_out = argv[2];
  const char *nm_array = argv[3];
  const char *type = argv[4];

  size_t ifarg = 5;

  // Array of vtk objects
  std::vector<TMeshType *> mesh;

  // Read the reference object
  TMeshType *mref = ReadMesh<TMeshType>(fn_ref);

  // Read each of the meshes
  for(size_t i = ifarg; i < argc; i++)
    mesh.push_back(ReadMesh<TMeshType>(argv[i]));

  // Find the reference array
  vtkDataArray *arr0 = NULL;
  if(strcmp(type,"cell")==0)
    arr0 = mesh.front()->GetCellData()->GetArray(nm_array);
  else if(strcmp(type,"point")==0)
    arr0 = mesh.front()->GetPointData()->GetArray(nm_array);

  std::cout << "Array " << arr0->GetName() << " with " << arr0->GetNumberOfTuples() << " tuples" << std::endl;

  if(!arr0)
    return usage();

  // Create an output array in the object
  vtkFloatArray *arr = vtkFloatArray::New();
  arr->SetNumberOfComponents(mesh.size());
  arr->SetNumberOfTuples(arr0->GetNumberOfTuples());
  arr->SetName(arr0->GetName());

  // Search all filenames
  for(int i = 0; i < mesh.size(); i++)
  {
    vtkDataArray *arri = NULL;
    if(strcmp(type,"cell")==0)
      arri = mesh[i]->GetCellData()->GetArray(nm_array);
    else if(strcmp(type,"point")==0)
      arri = mesh[i]->GetPointData()->GetArray(nm_array);

    for(int k = 0; k < arr0->GetNumberOfTuples(); k++)
    {
    arr->SetComponent(k, i, arri->GetComponent(k, 0));
    }
  }

  // Attach the new array to the reference dataset
  if(strcmp(type,"cell")==0)
    mref->GetCellData()->AddArray(arr);
  else if(strcmp(type,"point")==0)
    mref->GetPointData()->AddArray(arr);

  // Save the reference mesh
  WriteMesh<TMeshType>(mref, fn_out);

  return 0;
}

int main(int argc, char *argv[])
{
  return MyMain<vtkUnstructuredGrid>(argc, argv);
}

