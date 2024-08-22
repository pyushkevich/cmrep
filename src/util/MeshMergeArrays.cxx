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

int mesh_merge_arrays_usage()
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
  cout << "  -m <file>      Read mesh filenames from a file instead of command line" << endl;
  return -1;
}

struct Parameters {
  std::string arr_name, out_arr_name, fnout, fnref, fnsel, fnmeshfile;
  bool flag_cell = false;
  bool flag_binary = false;
  vector<string> source_meshes;
};


template <class TMeshType>
int MeshMergeArrays(const Parameters &p)
{
  // Read the reference mesh
  TMeshType *ref = ReadMesh<TMeshType>(p.fnref.c_str());

  // Read each of the input meshes
  std::vector<vtkDataArray *> da;
  unsigned int comp_total = 0;
  for(auto mesh_fn: p.source_meshes)
    {
    TMeshType *src = ReadMesh<TMeshType>(mesh_fn.c_str());

    // Get the corresponding array
    int idx = 0;
    vtkDataArray *arr = (p.flag_cell) 
      ? src->GetCellData()->GetArray(p.arr_name.c_str(), idx)
      : src->GetPointData()->GetArray(p.arr_name.c_str(), idx);

    // If no array, crap out
    if (arr == NULL)
      {
      cerr << "Missing array " << p.arr_name << " in mesh " << mesh_fn << endl;
      return -1;
      }

    // Get number of components
    unsigned int nc = arr->GetNumberOfComponents();
    cout << "Read " << mesh_fn << " with " << nc << " components." << std::endl;

    // Push back the array
    da.push_back(arr);
    comp_total += nc;
    }

  // Read the selection file
  vector<bool> selection;
  unsigned int comp_selected = 0;
  if(p.fnsel.size())
    {
    ifstream sfile(p.fnsel);
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
  array->SetNumberOfTuples(p.flag_cell ? ref->GetNumberOfCells() : ref->GetNumberOfPoints());

  // Read each of the input meshes and merge their arrays
  int target_index = 0, source_index = 0;
  for(auto *da_j : da)
    {
    // Add array to main accumulator
    for(unsigned int q = 0; q < da_j->GetNumberOfComponents(); q++)
      {
      if(selection.size() == 0 || selection[source_index])
        {
        for(int k = 0; k < da_j->GetNumberOfTuples(); k++)
          array->SetComponent(k, target_index, da_j->GetComponent(k, q));
        target_index++;
        }
      source_index++;
      }
    }

  // Stick the array in the output
  array->SetName(p.out_arr_name.length() > 0 ? p.out_arr_name.c_str() : p.arr_name.c_str());
  if(p.flag_cell)
    ref->GetCellData()->AddArray(array);
  else
    ref->GetPointData()->AddArray(array);

  // Write the output
  WriteMesh<TMeshType>(ref, p.fnout.c_str(), p.flag_binary);
  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  int i;
  Parameters p;
  
  // Read the optional parameters
  for (i = 1; i < argc; i++)
    {
    if (0 == strcmp(argv[i],"-o"))
      p.out_arr_name = argv[++i];
    else if (0 == strcmp(argv[i],"-r"))
      p.fnref = argv[++i];
    else if (0 == strcmp(argv[i],"-s"))
      p.fnsel = argv[++i];
    else if (0 == strcmp(argv[i],"-c"))
      p.flag_cell = true;
    else if (0 == strcmp(argv[i],"-B"))
      p.flag_binary = true;
    else if (0 == strcmp(argv[i],"-m"))
      p.fnmeshfile = argv[++i];
    else
      break;
    }

  // Check that there is enough input left
  int n_req_args = p.fnmeshfile.length() ? 2 : 3;
  if(i > argc - n_req_args)
    {
    cerr << "Not enough arguments provided on the command line, see usage" << endl;
    return mesh_merge_arrays_usage();
    }

  // Get the output filename
  p.fnout = argv[i++];
  
  // Get the input array
  p.arr_name = argv[i++];

  // Start of the input filenames
  int nin = argc - i;

  // Did the user provide an input file with mesh names?
  if(p.fnmeshfile.size())
    {
    if(nin > 0)
      {
      cerr << "Meshes can be specified via file (-m) or command line, but not both" << endl;
      return -1;
      }

    // Read the input file
    std::ifstream inputFile(p.fnmeshfile.c_str());
    if (!inputFile.is_open())
      {
      cerr << "Error opening file " << p.fnmeshfile << endl;
      return -1;
      }

    std::string line;
    while (std::getline(inputFile, line) && line.length() > 0)
      p.source_meshes.push_back(line);
    }
  else
    {
    for(int j = 0; j < nin; j++)
      p.source_meshes.push_back(argv[i + j]);
    }

  // Read the reference mesh
  if(p.fnref.length() == 0)
    p.fnref = p.source_meshes.front();

  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(p.fnref.c_str());
  reader->OpenVTKFile();
  reader->ReadHeader();
  bool is_ug = reader->IsFileUnstructuredGrid(), is_pd = reader->IsFilePolyData();
  reader->Delete();

  // Is this a polydata?
  if(is_ug)
    {
    return MeshMergeArrays<vtkUnstructuredGrid>(p);
    }
  else if(is_pd)
    {
    return MeshMergeArrays<vtkPolyData>(p);
    }
  else
    {
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}

