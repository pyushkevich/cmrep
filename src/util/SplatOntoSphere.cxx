#include <iostream>

#include <vtkPlatonicSolidSource.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkArrayCalculator.h>
#include <vtkCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>

using namespace std;

int usage() 
{
  cout << "sphere_splat: compute density of points on a sphere" << endl;
  cout << "usage:" << endl;
  cout << "  sphere_splat [options] inputs.vtk" << endl;
  cout << "options:" << endl;

  return -1;
}

int main(int argc, char *argv[])
{
  unsigned int n_tess = 5;

  vtkSmartPointer<vtkPlatonicSolidSource> ico = vtkPlatonicSolidSource::New();
  ico->SetSolidTypeToIcosahedron();
  ico->Update();

  vtkSmartPointer<vtkLoopSubdivisionFilter> subdivide = vtkSmartPointer<vtkLoopSubdivisionFilter>::New ();
  subdivide->SetNumberOfSubdivisions (n_tess);
  subdivide->SetInputConnection (ico->GetOutputPort ());

  subdivide->Update();

  // Normalize to unit sphere
  vtkSmartPointer<vtkArrayCalculator> calc = vtkArrayCalculator::New();
  calc->SetInputConnection(subdivide->GetOutputPort());
  calc->AddCoordinateVectorVariable("x");
  calc->SetCoordinateResults(1);
  calc->SetFunction("norm(x)");
  calc->Update();

  vtkSmartPointer<vtkPolyData> sphere = calc->GetPolyDataOutput();

  // Create a cell data array
  vtkSmartPointer<vtkDoubleArray> hits = vtkDoubleArray::New();
  hits->SetNumberOfComponents(1);
  hits->SetNumberOfTuples(sphere->GetNumberOfCells());
  hits->FillComponent(0, 0.0);
  hits->SetName("hits");
  sphere->GetCellData()->AddArray(hits);
  

  // Create a locator
  vtkSmartPointer<vtkCellLocator> loc = vtkCellLocator::New();
  loc->SetDataSet(sphere);
  loc->BuildLocator();

  // Process all of the input meshes
  for(int i = 2; i < argc;i++)
    {
    // Read mesh
    vtkSmartPointer<vtkPolyDataReader> reader = vtkPolyDataReader::New();
    reader->SetFileName(argv[i]);
    reader->Update();
    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();

    // Iterate over points
    for(int j = 0; j < mesh->GetNumberOfPoints(); j++)
      {
      double cp[3];
      vtkSmartPointer<vtkGenericCell> cell = vtkGenericCell::New();
      vtkIdType cellid;
      int subid;
      double dist2;

      // Run the locator
      loc->FindClosestPoint(mesh->GetPoint(j), cp, cell, cellid, subid, dist2);

      // Add a hit to the cell
      hits->SetTuple1(cellid, hits->GetTuple1(cellid) + 1);
      }

    cout << "." << flush;
    mesh->Delete();

    }

  cout << endl;

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
  writer->SetFileName(argv[1]);
  writer->SetInputData(sphere);
  writer->Update();

}
