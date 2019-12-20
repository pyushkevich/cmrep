#include "TetraLib.h"

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkTetra.h>
#include <vtkCellData.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>

using namespace std;

typedef vnl_matrix_fixed<double, 4, 4> MatrixType;

void WriteTetgenOutputAsUnstructuredMesh(tetgenio &out, const char *fname)
{
  vtkPoints* outpoints = vtkPoints::New();
  vtkFloatArray* pointarray = vtkFloatArray::New();
  vtkFloatArray* cellarray  = vtkFloatArray::New();
  vtkUnstructuredGrid* outmesh = vtkUnstructuredGrid::New();

  // Add the points 
  outpoints->SetNumberOfPoints(out.numberofpoints);
  pointarray->SetName ("Point Volume");
  pointarray->Allocate(out.numberofpoints);
  cout <<  out.numberofpoints << " points in tetgen mesh" << endl;
  cout <<  out.numberoftetrahedronattributes << " cell arrays in tetgen mesh" << endl;
  cout <<  out.numberofpointattributes << " point arrays in tetgen mesh" << endl;
  for(size_t k = 0; k < out.numberofpoints; k++)
    {
    double pt[3];
    pt[0] = out.pointlist[ k*3 ];
    pt[1] = out.pointlist[ k*3 + 1 ];
    pt[2] = out.pointlist[ k*3 + 2 ];
    outpoints->SetPoint(k, pt[0], pt[1], pt[2]);
    pointarray->InsertNextValue(0.0);
    }
  outmesh->SetPoints(outpoints);
  
  // Debug
  /*
  for(size_t k = 0; k < outpoints->GetNumberOfPoints(); k++)
    {
    double *pt = outpoints->GetPoint(k);
    cout << " point " << k << " is " << pt << endl;
    }
  */

  // Allocate the mesh
  outmesh->Allocate(out.numberoftetrahedra);
  
  cellarray->SetName ("Cell Volume");
  cellarray->Allocate(out.numberoftetrahedra);
  // Add the cells
  for(size_t k = 0; k < out.numberoftetrahedra; k++)
    {
    vtkIdList* idlist = vtkIdList::New();
    unsigned int ids[4];
    ids[0] = out.tetrahedronlist[ k*4 ];
    ids[1] = out.tetrahedronlist[ k*4 + 1 ];
    ids[2] = out.tetrahedronlist[ k*4 + 2 ];
    ids[3] = out.tetrahedronlist[ k*4 + 3 ];
    
    idlist->InsertNextId (ids[0]);
    idlist->InsertNextId (ids[1]);
    idlist->InsertNextId (ids[2]);
    idlist->InsertNextId (ids[3]);
    outmesh->InsertNextCell (VTK_TETRA, idlist);

    // Calculate the volume
    double p[4][3];
    outmesh->GetPoint(ids[0], p[0]);
    outmesh->GetPoint(ids[1], p[1]);
    outmesh->GetPoint(ids[2], p[2]);
    outmesh->GetPoint(ids[3], p[3]);
    double vol = vtkTetra::ComputeVolume(p[0], p[1], p[2], p[3]);
    cellarray->InsertNextValue(vol);

    // Distribute the volume to all vertices equally
    for(size_t ivert = 0; ivert < 4; ivert++)
       pointarray->SetTuple1(ids[ivert], pointarray->GetTuple1(ids[ivert]) + vol/4.0);
    idlist->Delete();
    }
  if (outmesh->GetCellData())
    {
    cout << "Adding cell data.." << endl;
    outmesh->GetCellData()->AddArray(cellarray);
    }
    
  if (outmesh->GetPointData())
    {
    cout << "Adding point data.." << endl;
    outmesh->GetPointData()->AddArray(pointarray);
    }

  cout <<  outmesh->GetNumberOfPoints() << " points in output mesh" << endl;
  // Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
  /*
  out.save_nodes("vtkout");
  out.save_poly("vtkout");
  out.save_elements("vtkout");
  out.save_faces("vtkout");
  */

  // Write the vtk output
  vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fname);
  writer->SetInputData(outmesh);
  writer->Write();
  writer->Delete();
}
