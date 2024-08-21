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
  unsigned int n_tess = 3;

  /*
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
  
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
  writer->SetFileName("sphere_canon.vtk");
  writer->SetInputData(sphere);
  writer->Update();
  */

  vtkSmartPointer<vtkPolyDataReader> reader = vtkPolyDataReader::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  vtkSmartPointer<vtkPolyData> sphere = reader->GetOutput();

  // Process all of the input meshes
  for(int i = 2; i < argc;i++)
    {
    // Read mesh
    vtkSmartPointer<vtkPolyDataReader> reader = vtkPolyDataReader::New();
    reader->SetFileName(argv[i]);
    reader->Update();
    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();

    // Map the mesh onto the sphere
    vtkSmartPointer<vtkArrayCalculator> calcmesh = vtkArrayCalculator::New();
    calcmesh->SetInputData(mesh);
    calcmesh->AddScalarArrayName("theta");
    calcmesh->AddScalarArrayName("phi");
    calcmesh->SetCoordinateResults(1);
    calcmesh->SetFunction("sin(theta)*cos(phi)*iHat+sin(theta)*sin(phi)*jHat+cos(theta)*kHat");
    calcmesh->Update();

    vtkSmartPointer<vtkPolyData> mesh_sphere = calcmesh->GetPolyDataOutput();

    // Create a locator for the newly created sphere
    vtkSmartPointer<vtkCellLocator> loc = vtkCellLocator::New();
    loc->SetDataSet(mesh_sphere);
    loc->BuildLocator();

    // Updated sphere
    vtkSmartPointer<vtkPolyData> sphereCopy = vtkSmartPointer<vtkPolyData>::New();
    sphereCopy->DeepCopy(sphere);

    // Iterate over points in the landmark sphere
    for(int j = 0; j < sphere->GetNumberOfPoints(); j++)
      {
      double p[3];
      double cp[3];
      vtkSmartPointer<vtkGenericCell> cell = vtkGenericCell::New();
      vtkIdType cellid;
      int subid;
      double dist2;
      double pcoords[3], weights[3];

      // Compute the spherical coordinates
      loc->FindClosestPoint(sphere->GetPoint(j), cp, cell, cellid, subid, dist2);
      // cout << cp[0] << " " << cp[1] << " " << cp[2] << endl;
      cell->EvaluatePosition(cp, p, subid, pcoords, dist2, weights);
      cout << weights[0] << " " << weights[1] << " " << weights[2] << endl;

      // Splat the coordinates onto the vertices of the triangle
      sphereCopy->GetPoints()->SetPoint(j, 0, 0, 0);
      for(int k = 0; k < cell->GetNumberOfPoints(); k++)
        {
        vtkIdType m = cell->GetPointId(k);
        double xk = mesh->GetPoint(m)[0] * weights[k] + sphereCopy->GetPoint(j)[0];
        double yk = mesh->GetPoint(m)[1] * weights[k] + sphereCopy->GetPoint(j)[1];
        double zk = mesh->GetPoint(m)[2] * weights[k] + sphereCopy->GetPoint(j)[2];
        sphereCopy->GetPoints()->SetPoint(j, xk, yk, zk);
        }
      }

    char buffer[256];
    snprintf(buffer, 256, "mesh%04d.vtk", i);
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
    writer->SetFileName(buffer);
    writer->SetInputData(sphereCopy);
    writer->Update();
    
    cout << "." << flush;
    }

  cout << endl;
}
