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
  cout << "  -B             Save mesh in binary format (default: ASCII format)" << endl;
  cout << "  -s <file>      Only add selected arrays based on file. The file must" << endl;
  cout << "                 have a 0 (drop) or 1 (select) on a separate line for each" << endl;
  cout << "                 input array." << endl;
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
void WriteMesh(TMeshType *mesh, const char *fname, bool flag_write_binary = false)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname, bool flag_write_binary)
{
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(flag_write_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname, bool flag_write_binary)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(flag_write_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

template <class TMeshType>
int MeshMergeArrays(int argc, char *argv[])
{
  int i;
  std::string arr_name, out_arr_name, fnout, fnref, fnsel;
  bool flag_cell = false;
  bool flag_binary = false;
  
  // Read the optional parameters
  for (i = 1; i < argc; i++)
    {
    if (0 == strcmp(argv[i],"-o"))
      out_arr_name = argv[++i];
    else if (0 == strcmp(argv[i],"-r"))
      fnref = argv[++i];
    else if (0 == strcmp(argv[i],"-s"))
      fnsel = argv[++i];
    else if (0 == strcmp(argv[i],"-c"))
      flag_cell = true;
    else if (0 == strcmp(argv[i],"-B"))
      flag_binary = true;
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

  // Read each of the input meshes
  std::vector<vtkDataArray *> da;
  unsigned int comp_total = 0;
  for(int j = 0; j < argc - i; j++)
    {
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

    // Get number of components
    unsigned int nc = arr->GetNumberOfComponents();
    cout << "Read " << fnin[j] << " with " << nc << " components." << std::endl;

    // Push back the array
    da.push_back(arr);
    comp_total += nc;
    }

  // Read the selection file
  vector<bool> selection;
  unsigned int comp_selected = 0;
  if(fnsel.size())
    {
    ifstream sfile(fnsel);
    bool sel;
    while(sfile >> sel)
      {
      selection.push_back(sel);
      if(sel)
        comp_selected++;
      }

    if(selection.size() != comp_total)
      {
      cerr << "Selection file size " << selection.size() 
        << " does not match total components " << comp_total << endl;
      return -1;
      }
    }
  else
    comp_selected = comp_total;


  // Create output array
  cout << "Output array will contain " << comp_selected << " components" << endl;
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetNumberOfComponents(comp_selected);
  array->SetNumberOfTuples(flag_cell ? ref->GetNumberOfCells() : ref->GetNumberOfPoints());

  // Read each of the input meshes and merge their arrays
  int target_index = 0, source_index = 0;
  for(int j = 0; j < argc - i; j++)
    {
    // Add array to main accumulator
    for(unsigned int q = 0; q < da[j]->GetNumberOfComponents(); q++)
      {
      if(selection.size() == 0 || selection[source_index])
        {
        for(int k = 0; k < da[j]->GetNumberOfTuples(); k++)
          array->SetComponent(k, target_index, da[j]->GetComponent(k, q));
        target_index++;
        }
      source_index++;
      }
    }

  // Stick the array in the output
  array->SetName(out_arr_name.length() > 0 ? out_arr_name.c_str() : arr_name.c_str());
  if(flag_cell)
    ref->GetCellData()->AddArray(array);
  else
    ref->GetPointData()->AddArray(array);

  // Write the output
  WriteMesh<TMeshType>(ref, fnout.c_str(), flag_binary);
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

