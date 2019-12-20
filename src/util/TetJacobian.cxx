#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCell.h>
#include <vtkTetra.h>
#include <iostream>

using namespace std;

int usage()
{
  cout << 
    "tetjac - compute Jacobian between two tet-meshes\n"
    "usage:\n"
    "  tetjac mesh1.vtk mesh2.vtk output.vtk\n";
  return -1;
}


bool ComputeTetrahedralVolumes(vtkUnstructuredGrid *mesh, vtkFloatArray *cellarray, vtkFloatArray *pointarray, double &totalvol)
{
  totalvol = 0.0;
  cellarray->Allocate(mesh->GetNumberOfCells());
  pointarray->Allocate(mesh->GetNumberOfPoints());
  for(int j = 0; j < mesh->GetNumberOfPoints(); j++)
     pointarray->InsertNextValue(0.0);
  for(int j = 0; j < mesh->GetNumberOfCells(); j++)
    {
    vtkCell *cell =  mesh->GetCell(j);
    if(cell->GetCellType() == VTK_TETRA)
      {
      double p[4][3];
      mesh->GetPoint(cell->GetPointId(0), p[0]);
      mesh->GetPoint(cell->GetPointId(1), p[1]);
      mesh->GetPoint(cell->GetPointId(2), p[2]);
      mesh->GetPoint(cell->GetPointId(3), p[3]);
      double vol = vtkTetra::ComputeVolume(p[0], p[1], p[2], p[3]);
      cellarray->InsertNextValue(vol);
      totalvol += vol;
      for(size_t ivert = 0; ivert < 4; ivert++)
         pointarray->SetTuple1(cell->GetPointId(ivert), pointarray->GetTuple1(cell->GetPointId(ivert)) + vol/4.0);

      }
    else return false;
    }

  return true;
}

int main(int argc, char *argv[])
{
  // Process parameters
  if(argc != 4) return usage();
  char *fnFix = argv[1];
  char *fnMov = argv[2];
  char *fnOut = argv[3];

  // Load the meshes
  vtkUnstructuredGridReader *readFix = vtkUnstructuredGridReader::New();
  readFix->SetFileName(fnFix);
  readFix->Update();
  vtkUnstructuredGrid *fix = readFix->GetOutput();

  vtkUnstructuredGridReader *readMov = vtkUnstructuredGridReader::New();
  readMov->SetFileName(fnMov);
  readMov->Update();
  vtkUnstructuredGrid *mov = readMov->GetOutput();

  // Allocate the arrays
  vtkFloatArray *vc = vtkFloatArray::New();
  vtkFloatArray *vp = vtkFloatArray::New();
  vtkFloatArray *jc = vtkFloatArray::New();
  vtkFloatArray *jp = vtkFloatArray::New();
  vc->SetName("Cell Volume");
  vp->SetName("Point Volume");
  jc->SetName("Cell Jacobian");
  jp->SetName("Point Jacobian");

  // Compute the volumes of the two meshes
  double tvfix, tvmov;
  if(!ComputeTetrahedralVolumes(fix, vc, vp, tvfix) || !ComputeTetrahedralVolumes(mov, jc, jp, tvmov))
    {
    cerr << "Meshes contain non-tetrahedral elements. Aborting." << endl;
    return -1;
    }

  cout << "Fixed Mesh Volume: " << tvfix << endl;
  cout << "Moving Mesh Volume: " << tvmov << endl;

  // Compute the Jacobian by dividing by the fixed volume
  if(vc->GetNumberOfTuples() != jc->GetNumberOfTuples() 
    || vp->GetNumberOfTuples() != jp->GetNumberOfTuples())
    {
    cerr << "Mesh dimensions do not match. Aborting." << endl;
    return -1;
    }

  for(size_t i = 0; i < vc->GetNumberOfTuples(); i++)
    jc->SetTuple1(i, jc->GetTuple1(i) / vc->GetTuple1(i));

  for(size_t i = 0; i < vp->GetNumberOfTuples(); i++)
    jp->SetTuple1(i, jp->GetTuple1(i) / vp->GetTuple1(i));

  fix->GetCellData()->AddArray(vc);
  fix->GetCellData()->AddArray(jc);
  fix->GetPointData()->AddArray(vp);
  fix->GetPointData()->AddArray(jp);

  // Report change in volume
  printf("VOLUME_STATS: %12.8f %12.8f %12.8f %12.8f\n",
    tvfix, tvmov, tvmov/tvfix, (tvfix - tvmov) / tvfix);

  // Save mesh
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fnOut);
  writer->SetInputData(fix);
  writer->Update();
}
