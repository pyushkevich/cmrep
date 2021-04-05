#include <cstdio>
#include <iostream>
#include <vtkSmartPointer.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>
#include <vtkCellData.h>
#include "ReadWriteVTK.h"

using namespace std;

int usage()
{
  cout << "apply_remeshing_to_another: apply a remeshing operation to another mesh\n"
       << "usage:\n"
       << "  apply_remeshing_to_another [options] <reference> <remeshed> <target> <output>\n"
       << "parameters:\n"
       << "  reference:        A mesh was input to a remeshing operation\n"
       << "  remeshed:         A mesh that was output by a remeshing operation\n"
       << "  target:           A mesh that you want to be remeshed consistently\n"
       << "  output:           The output mesh\n";

  return -1;
}

int main(int argc, char *argv[])
{
  typedef vtkSmartPointer<vtkPolyData> PolyDataPtr;

  if(argc < 5)
    return usage();

  // Read the three meshes
  PolyDataPtr pd_ref = ReadVTKData(argv[argc-4]);
  PolyDataPtr pd_rem = ReadVTKData(argv[argc-3]);
  PolyDataPtr pd_trg = ReadVTKData(argv[argc-2]);

  // Create a locator
  pd_ref->BuildCells();
  vtkSmartPointer<vtkCellLocator> loc = vtkCellLocator::New();
  loc->SetDataSet(pd_ref);
  loc->BuildLocator();

  // Create output mesh
  PolyDataPtr pd_out = vtkPolyData::New();
  pd_out->SetPoints(vtkPoints::New());
  pd_out->GetPoints()->SetNumberOfPoints(pd_rem->GetNumberOfPoints());
  pd_out->GetCellData()->DeepCopy(pd_rem->GetCellData());

  // For every vertex in the remeshed mesh, find its barycentric coordinates in the reference mesh
  double tol2 = 1.0e-6;
  double pcoords[3];
  double weights[256];
  vtkSmartPointer<vtkGenericCell> cell_ptr = vtkGenericCell::New();
  for(unsigned int i = 0; i < pd_rem->GetNumberOfPoints(); i++)
    {
    double *x_rem = pd_rem->GetPoint(i);
    int indx = loc->FindCell(x_rem, tol2, cell_ptr, pcoords, weights);
    if(indx >= 0)
      {
      double x_new[] = {0., 0., 0.};
      for(unsigned int j = 0; j < cell_ptr->GetNumberOfPoints(); j++)
        {
        unsigned int v_j = cell_ptr->GetPointId(j);
        for(unsigned int d = 0; d < 3; d++)
          x_new[d] += weights[j] * pd_trg->GetPoint(v_j)[d];
        }
      pd_out->GetPoints()->SetPoint(i, x_new);
      }
    }

  // Save output mesh
  WriteVTKData(pd_out, argv[argc-1]);

  return 0;
}
