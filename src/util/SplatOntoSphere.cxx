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
#include <vtkPointData.h>

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
  
  vtkSmartPointer<vtkDoubleArray> phits = vtkDoubleArray::New();
  phits->SetNumberOfComponents(1);
  phits->SetNumberOfTuples(sphere->GetNumberOfPoints());
  phits->FillComponent(0, 0.0);
  phits->SetName("phits");
  sphere->GetPointData()->AddArray(phits);
  
  vtkSmartPointer<vtkDoubleArray> theta = vtkDoubleArray::New();
  theta->SetNumberOfComponents(1);
  theta->SetNumberOfTuples(sphere->GetNumberOfPoints());
  theta->FillComponent(0, 0.0);
  theta->SetName("theta");
  sphere->GetPointData()->AddArray(theta);
  
  vtkSmartPointer<vtkDoubleArray> phi = vtkDoubleArray::New();
  phi->SetNumberOfComponents(1);
  phi->SetNumberOfTuples(sphere->GetNumberOfPoints());
  phi->FillComponent(0, 0.0);
  phi->SetName("phi");
  sphere->GetPointData()->AddArray(phi);
  
  // Another array for mean shape
  vtkSmartPointer<vtkDoubleArray> meanshape = vtkDoubleArray::New();
  meanshape->SetNumberOfComponents(3);
  meanshape->SetNumberOfTuples(sphere->GetNumberOfPoints());
  meanshape->FillComponent(0, 0.0);
  meanshape->FillComponent(1, 0.0);
  meanshape->FillComponent(2, 0.0);
  meanshape->SetName("meanshape");
  sphere->GetPointData()->AddArray(meanshape);

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

    // Get the phi and theta arrays
    vtkDataArray *phi = mesh->GetPointData()->GetArray("phi");
    vtkDataArray *theta = mesh->GetPointData()->GetArray("theta");

    // Iterate over points
    for(int j = 0; j < mesh->GetNumberOfPoints(); j++)
      {
      double p[3];
      double cp[3];
      vtkSmartPointer<vtkGenericCell> cell = vtkGenericCell::New();
      vtkIdType cellid;
      int subid;
      double dist2;

      // Compute the spherical coordinates
      double theta_x = theta->GetComponent(j,0), phi_x = phi->GetComponent(j,0);
      p[0] = sin(theta_x) * cos(phi_x);
      p[1] = sin(theta_x) * sin(phi_x);
      p[2] = cos(theta_x);

      // Run the locator
      loc->FindClosestPoint(p, cp, cell, cellid, subid, dist2);

      // Add a hit to the cell
      hits->SetTuple1(cellid, hits->GetTuple1(cellid) + 1);

      // Splat the coordinates onto the vertices of the triangle
      for(int k = 0; k < cell->GetNumberOfPoints(); k++)
        {
        vtkIdType m = cell->GetPointId(k);
        for(int d = 0; d < 3; d++)
          {
          meanshape->SetComponent(m, d, meanshape->GetComponent(m, d) + mesh->GetPoint(j)[d]); 
          phits->SetTuple1(m, 1 + phits->GetTuple1(m));
          }
        }
      }

    cout << "." << flush;
    }

  cout << endl;

  // Update the arrays
  for(int j = 0; j < sphere->GetNumberOfPoints(); j++)
    {
    double *x = sphere->GetPoint(j);
    theta->SetTuple1(j, acos(x[2] / (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])));
    phi->SetTuple1(j, atan2(x[1],x[0]));

    double n = phits->GetTuple1(j);
    if(n > 0)
      {
      for(int d = 0; d < 3; d++)
        {
        meanshape->SetComponent(j, d, meanshape->GetComponent(j, d) / n);
        }
      sphere->GetPoints()->SetPoint(j, 
        meanshape->GetComponent(j,0),
        meanshape->GetComponent(j,1),
        meanshape->GetComponent(j,2));
      }
    }

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
  writer->SetFileName(argv[1]);
  writer->SetInputData(sphere);
  writer->Update();

}
