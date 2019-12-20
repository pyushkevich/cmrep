#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyDataReader.h>
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBYUWriter.h>
#include <vtkSTLWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

int usage()
{
  cout << "tetdumpattr - Dumps an attribute array in VTK tetra mesh" << endl;
  cout << "usage: " << endl;
  cout << "   dumpmeshattr [options] input.vtk array_name output.dat" << endl;
  cout << "options: " << endl;
  cout << "   -b    Dump array as binary" << endl;
  cout << "   -c    Request cell array (by default, point arrays are used" << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 3) return usage();
  
  // Read the appropriate mesh
  vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
  reader->SetFileName(argv[argc-3]);
  reader->Update();
  vtkUnstructuredGrid *grid = reader->GetOutput();

  // Read the options
  bool flagBinary = false;
  bool flagCell = false;
  for(size_t i = 1; i < argc-3; i++)
    {
    if(!strcmp(argv[i], "-b"))
      flagBinary = true;
    else if(!strcmp(argv[i], "-c"))
      flagCell = true;
    else 
      { 
      cerr << "Bad option" << argv[i] << endl;
      return usage();
      }
    }

  // Get the array
  vtkDataArray *array = NULL;
  if(flagCell)
    array = grid->GetCellData()->GetArray(argv[argc-2]);
  else
    array = grid->GetPointData()->GetArray(argv[argc-2]);

  // Check that the array exists
  if(!array)
    {
    cerr << "Array is missing in input data" << endl;
    return -1;
    }

  // Dump the array
  if(flagBinary)
    {
    FILE *f = fopen(argv[argc-1], "wb");
    fwrite(array->GetVoidPointer(0), array->GetDataTypeSize(), 
      array->GetNumberOfComponents() * array->GetNumberOfTuples(), f);
    fclose(f);
    }
  else
    {
    // Write the array values
    ofstream fout(argv[argc-1]);
    for(size_t iPoint = 0; iPoint < array->GetNumberOfTuples(); iPoint++)
      {
      for(size_t iComp = 0; iComp < array->GetNumberOfComponents(); iComp++)
        {
        fout << array->GetComponent(iPoint, iComp);
        if(iComp + 1 < array->GetNumberOfComponents())
          fout << " ";
        else 
          fout << endl;
        }
      }
    fout.close();
    }
}
