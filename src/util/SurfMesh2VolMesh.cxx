#include <iostream>
#include <string>
#include <sstream>
#include <set>

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
//#include <itkOrientedRASImage.h>
#include <itkImage.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>
#include "ReadWriteVTK.h"
#include "TetraLib.h"

using namespace std;
using namespace itk;

typedef vnl_matrix_fixed<double, 4, 4> MatrixType;

int usage()
{
  cout << "surf2tet - Converts a VTK surface mesh to a VTK volume mesh" << endl;
  cout << "usage: " << endl;
  cout << "   surf2tet meshin.vtk meshout.vtk [tetgenoptions]" << endl;
  cout << "tetgenoptions: default is pq1.414a0.1" << endl;
  return -1;
}

int main(int argc, char **argv)
{
  
  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;


  // Check the parameters
  if(argc < 3) return usage();

  // tetgen option 
  std::string optstr = std::string("");
  if (argc > 3)
    {
    for(int i = 3; i < argc; i++)
      optstr = optstr + std::string( argv[i] );
    }
  else
    optstr = optstr + std::string("pq1.414a0.1");


  // Read the mesh
  vtkPolyData *mesh = ReadVTKData(argv[1]);

  // Clean the mesh
  vtkTriangleFilter *triangulator = vtkTriangleFilter::New();
  triangulator->SetInputData(mesh);
  triangulator->PassLinesOff();
  triangulator->PassVertsOff();
  triangulator->Update();
  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInputConnection(triangulator->GetOutputPort());
  cleaner->Update();
  vtkTriangleFilter *triangulator2 = vtkTriangleFilter::New();
  triangulator2->SetInputConnection(cleaner->GetOutputPort());
  triangulator2->PassLinesOff();
  triangulator2->PassVertsOff();
  triangulator2->Update();
  mesh = triangulator2->GetOutput();

  // Cast mesh into tetgen format
  // All indices start from 0.
  in.firstnumber = 0;

  // Enter the nodes
  in.numberofpoints = mesh->GetNumberOfPoints();
  in.pointlist = new REAL[in.numberofpoints * 3];
  for(size_t k = 0; k < mesh->GetNumberOfPoints(); k++)
    {
    double *pt = mesh->GetPoint(k);
    in.pointlist[ k*3 ] = pt[0];
    in.pointlist[ k*3 + 1 ] = pt[1];
    in.pointlist[ k*3 + 2 ] = pt[2];
    }
  
  // Enter the cells
  in.numberoffacets = mesh->GetNumberOfCells(); 
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  for(size_t k = 0; k < mesh->GetNumberOfCells(); k++)
    {
    vtkCell *cell = mesh->GetCell(k);
    if(cell->GetCellType() != VTK_TRIANGLE)
      throw("Wrong cell type");
    vtkIdType a0 = cell->GetPointId(0);
    vtkIdType a1 = cell->GetPointId(1);
    vtkIdType a2 = cell->GetPointId(2);
    f = &in.facetlist[k];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = a0;
    p->vertexlist[1] = a1;
    p->vertexlist[2] = a2;


    }

  // Output the PLC to files 'barin.node' and 'barin.poly'.
  /*
  in.save_nodes("vtkin");
  in.save_poly("vtkin");
  in.save_elements("vtkin");
  in.save_faces("vtkin");
  */

  //tetrahedralize("pq1.414a0.1", &in, &out);
  tetgenbehavior topt;
  topt.parse_commandline((char *) optstr.c_str());
  tetrahedralize(&topt, &in, &out);
  

  // Write to a mesh
  WriteTetgenOutputAsUnstructuredMesh(out, argv[2]);

  return 0;
}
