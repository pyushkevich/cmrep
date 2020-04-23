#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkTriangle.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
//#include <itkOrientedRASImage.h>
#include <itkImage.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>
#include "ReadWriteVTK.h"

#define TRILIBRARY
#define ANSI_DECLARATORS

#include "triangle.c"


using namespace std;
using namespace itk;

typedef vnl_matrix_fixed<double, 4, 4> MatrixType;

int usage()
{
  cout << "contour2surf - Converts a VTK contour mesh to a VTK surface mesh" << endl;
  cout << "usage: " << endl;
  cout << "   contour2surf meshin.vtk meshout.vtk [triangleoptions]" << endl;
  cout << "triangleoptions: default is zpq32a160" << endl;
  return -1;
}

void WriteTriangleOutputAsPolyDataMesh(triangulateio &out, const char *fname)
{
  vtkPoints* outpoints = vtkPoints::New();
  vtkFloatArray* pointarray = vtkFloatArray::New();
  vtkFloatArray* cellarray  = vtkFloatArray::New();
  vtkPolyData* outmesh = vtkPolyData::New();

  // Add the points 
  outpoints->SetNumberOfPoints(out.numberofpoints);
  pointarray->SetName ("Point Area");
  pointarray->Allocate(out.numberofpoints);
  cout <<  out.numberofpoints << " points in triangle mesh" << endl;
  cout <<  out.numberoftriangleattributes << " cell arrays in triangle mesh" << endl;
  cout <<  out.numberofpointattributes << " point arrays in triangle mesh" << endl;
  for(int k = 0; k < out.numberofpoints; k++)
    {
    double pt[3];
    pt[0] = out.pointlist[ k*2 ];
    pt[1] = out.pointlist[ k*2 + 1 ];
    pt[2] = 0.0;
    outpoints->SetPoint(k, pt[0], pt[1], pt[2]);
    pointarray->InsertNextValue(0.0);
    }
  outmesh->SetPoints(outpoints);


  // Allocate the mesh
  outmesh->Allocate(out.numberoftriangles);

  cellarray->SetName ("Cell Area");
  cellarray->Allocate(out.numberoftriangles);
  // Add the cells 
  for(int k = 0; k < out.numberoftriangles; k++)
    {
    vtkIdList* idlist = vtkIdList::New();
    unsigned int ids[4];
    ids[0] = out.trianglelist[ k*3 ];
    ids[1] = out.trianglelist[ k*3 + 1 ];
    ids[2] = out.trianglelist[ k*3 + 2 ];

    idlist->InsertNextId (ids[0]);
    idlist->InsertNextId (ids[1]);
    idlist->InsertNextId (ids[2]);
    outmesh->InsertNextCell (VTK_TRIANGLE, idlist);

    // Calculate the area
    double p[3][3];
    outmesh->GetPoint(ids[0], p[0]);
    outmesh->GetPoint(ids[1], p[1]);
    outmesh->GetPoint(ids[2], p[2]);
    double area = vtkTriangle::TriangleArea(p[0], p[1], p[2]);
    cellarray->InsertNextValue(area);

    // Distribute the area
    for(size_t ivert = 0; ivert < 3; ivert++)
       pointarray->SetTuple1(ids[ivert], pointarray->GetTuple1(ids[ivert]) + area/3.0);
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
  WriteVTKData( outmesh, fname);
//    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//    writer->SetFileName(fname);
//    writer->SetInput(outmesh);
//    writer->Update();

}



int main(int argc, char **argv)
{
  
  struct triangulateio in, out, vor;


  // Check the parameters
  if(argc < 3) return usage();

  // triangle option 
  std::string optstr = std::string("");
  if (argc > 3)
    {
    for(int i = 3; i < argc; i++)
      optstr = optstr + std::string( argv[i] );
    }
  else
    optstr = optstr + std::string("zpq32a160");


  // Read the mesh
  vtkPolyData *mesh = ReadVTKData(argv[1]);

  // Cast mesh into triangle format
  // All indices start from 0.
  //in.firstnumber = 0; make sure to pass in -z option.

  // Enter the nodes
  in.numberofpoints = mesh->GetNumberOfPoints();
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

  for(int k = 0; k < mesh->GetNumberOfPoints(); k++)
    {
    double *pt = mesh->GetPoint(k);
    in.pointlist[ k*2 ] = pt[0];
    in.pointlist[ k*2 + 1 ] = pt[1];
    //cout << k << " " << pt[0] << " " << pt[1] << endl;
    }
  
  in.numberofpointattributes = 0;
  in.pointmarkerlist = NULL;

  // Enter the cells
  in.numberofsegments = mesh->GetNumberOfCells(); 
  in.segmentmarkerlist = NULL;
  in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  for(int k = 0; k < mesh->GetNumberOfCells(); k++)
    {
    vtkCell *cell = mesh->GetCell(k);
//    if(cell->GetCellType() != VTK_TRIANGLE)
//      throw("Wrong cell type");
    vtkIdType a0 = cell->GetPointId(0);
    vtkIdType a1 = cell->GetPointId(1);
    in.segmentlist[k*2] = a0; 
    in.segmentlist[k*2 + 1] = a1; 
    //cout << k << " " << a0 << " " << a1 << endl;

    }

    in.numberofholes = 0;
    in.numberofregions = 0;


  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out' */

  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
  out.pointattributelist = (REAL *) NULL;
  out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  out.triangleattributelist = (REAL *) NULL;
  out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
  /* Needed only if segments are output (-p or -c) and -P not used: */
  out.segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  out.segmentmarkerlist = (int *) NULL;
  out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */


  cout << "Options are: " << optstr << endl;
  triangulate((char *)optstr.c_str(), &in, &out, &vor);
  

  // Write to a mesh
  WriteTriangleOutputAsPolyDataMesh(out, argv[2]);

  return 0;
}
