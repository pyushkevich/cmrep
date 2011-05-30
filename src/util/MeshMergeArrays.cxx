#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include "ReadWriteVTK.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "Usage: mesh_merge_arrays [options] output_mesh array_name input_meshes" << endl;
  cout << "Options: " << endl;
  cout << "  -n name        Name for output array (default: same as input)" << endl;
  cout << "  -r file.vtk    Use file as reference geometry (default: first mesh)" << endl;
  cout << "  -c             Apply to cell arrays (default: point arrays)" << endl;
  return -1;
}

template <class TMeshType>
TMeshType * ReadMesh(const char *fname)
{ return NULL; }

template <>
vtkUnstructuredGrid *ReadMesh<>(const char *fname)
{
  vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}

template <>
vtkPolyData *ReadMesh<>(const char *fname)
{
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}


template <class TMeshType>
void WriteMesh(TMeshType *mesh, const char *fname)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname)
{
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fname);
  writer->SetInput(mesh);
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fname);
  writer->SetInput(mesh);
  writer->Update();
}

template <class TMeshType>
int MeshMergeArrays(int argc, char *argv[])
{
  int i;
  std::string arr_name, out_arr_name, fnout, fnref;
  bool flag_cell = false;
  
  // Read the optional parameters
  for (i = 1; i < argc; i++)
    {
    if (0 == strcmp(argv[i],"-o"))
      out_arr_name = argv[++i];
    else if (0 == strcmp(argv[i],"-r"))
      fnref = argv[++i];
    else if (0 == strcmp(argv[i],"-c"))
      flag_cell = true;
    else
      break;
    }

  // Check that there is enough input left
  if(i > argc-3)
    return usage();

  // Get the output filename
  fnout = argv[i++];
  
  // Get the input array
  arr_name = argv[i++];

  // Start of the input filenames
  char **fnin = argv + i;
  int nin = argc - i;

  // Read the reference mesh
  TMeshType *ref = fnref.length() ? ReadMesh<TMeshType>(fnref.c_str()) : ReadMesh<TMeshType>(fnin[0]);

  // Create output array
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetNumberOfComponents(nin);
  array->SetNumberOfTuples(flag_cell ? ref->GetNumberOfCells() : ref->GetNumberOfPoints());

  // Read each of the input meshes
  for(int j = 0; j < argc - i; j++)
    {
    // Read one of the meshes
    TMeshType *src = ReadMesh<TMeshType>(fnin[j]);

    // Get the corresponding array
    int idx = 0;
    vtkDataArray *arr = (flag_cell) 
      ? src->GetCellData()->GetArray(arr_name.c_str(), idx)
      : src->GetPointData()->GetArray(arr_name.c_str(), idx);

    // If no array, crap out
    if (arr == NULL)
      {
      cerr << "Missing array " << arr_name << " in mesh " << fnin << endl;
      return -1;
      }

    // Add array to main accumulator
    for(int k = 0; k < arr->GetNumberOfTuples(); k++)
      array->SetComponent(k, j, arr->GetComponent(k, 0));

    // Done with it!
    src->Delete();
    }

  // Stick the array in the output
  array->SetName(out_arr_name.length() > 0 ? out_arr_name.c_str() : arr_name.c_str());
  if(flag_cell)
    ref->GetCellData()->AddArray(array);
  else
    ref->GetPointData()->AddArray(array);

  // Write the output
  WriteMesh<TMeshType>(ref, fnout.c_str());
  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  if(argc < 4)
    return usage();

  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(argv[argc-1]);
  reader->OpenVTKFile();
  reader->ReadHeader();

  bool isPolyData = true;

  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    isPolyData = false;
    return MeshMergeArrays<vtkUnstructuredGrid>( argc, argv );
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return MeshMergeArrays<vtkPolyData>( argc, argv );
    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}

