#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkThresholdPoints.h>
#include <vtkBoundingBox.h>
#include <vtkCellArray.h>
#include <vtkCellDataToPointData.h>
#include <vtkPolyData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBYUWriter.h>
#include <vtkSTLWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkMath.h>
#include <vtkTriangle.h>
#include <VTKMeshShortestDistance.h>
#include <VTKMeshHalfEdgeWrapper.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_cross.h>

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

using namespace std;

typedef std::pair<vtkIdType, vtkIdType> VertexPair;
typedef std::set< std::pair<vtkIdType, vtkIdType> > VertexPairSet;
typedef std::vector<VertexPair> VertexPairArray;

vtkPolyData *ReadVTKData(string fn)
{
  // Choose the reader based on extension
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(fn.c_str());
  reader->Update();
  return reader->GetOutput();
}

void WriteVTKData(vtkPolyData *data, string fn)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fn.c_str());
  writer->SetInput(data);
  writer->Update();
}
vtkPolyData *ReadVoronoiOutput(
  string fn,                  // Filename from which to read the points
  vtkBoundingBox *bbox,
  vtkPolyData *bnd,
  VertexPairArray &src,       // Output: for each cell in the VD, the pair of generators
  double xTol)      
{
  // Load the file
  ifstream fin(fn.c_str());

  // Load the numbers
  size_t nv, np, junk;
  
  // First two lines
  fin >> junk;
  fin >> nv; 

  vtkSelectEnclosedPoints *sel = vtkSelectEnclosedPoints::New();
  sel->SetTolerance(xTol);
  sel->Initialize(bnd);

  // Create an array of points
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(nv);

  // Create an array of in/out flags
  bool *ptin = new bool[nv];

  // Progress bar
  cout << "Selecting points inside mesh (n = " << nv << ")" << endl;
  cout << "|         |         |         |         |         |" << endl;
  size_t next_prog_mark = nv / 50;

  for(size_t i = 0; i < nv; i++)
    {
    double x,y,z;
    fin >> x;
    fin >> y;
    fin >> z;
    pts->SetPoint(i,x,y,z);

    // Is this point outside of the bounding box
    ptin[i] = bbox->ContainsPoint(x,y,z) && sel->IsInsideSurface(x,y,z);
    
    if(i >= next_prog_mark)
      {
      cout << "." << flush;
      next_prog_mark += nv / 50;
      }
    }
  cout << "." << endl;

  // Read the number of cells
  fin >> np;

  // Clear the src array
  src.clear();

  // Create the polygons 
  vtkCellArray *cells = vtkCellArray::New();
  for(size_t j = 0; j < np; j++)
    {
    bool isinf = false;
    bool isout = false;
    
    size_t m; fin >> m; m-=2;
    vtkIdType ip1, ip2; fin >> ip1; fin >> ip2;

    vtkIdType *ids = new vtkIdType[m];
    for(size_t k = 0; k < m; k++)
      {
      fin >> ids[k];

      // Is this point at infinity?
      if(ids[k] == 0) isinf = true; else ids[k]--;
      if(!ptin[ids[k]]) isout = true;
      }

    if(!isinf && !isout)
      {


      // add the cell
      cells->InsertNextCell(m, ids);

      // add the pair of generators
      src.push_back(make_pair(ip1, ip2)); // TODO: is this numbering 0-based?
      }

    delete ids;
    }

  // Create the vtk poly data
  vtkPolyData *poly = vtkPolyData::New();
  poly->SetPoints(pts);
  poly->SetPolys(cells);
  return poly;
}

inline vtkFloatingPointType TriangleArea(
  const vnl_vector_fixed<vtkFloatingPointType,3> &A, 
  const vnl_vector_fixed<vtkFloatingPointType,3> &B, 
  const vnl_vector_fixed<vtkFloatingPointType,3> &C)
{
  return 0.5 * vnl_cross_3d(B - A, C - A).magnitude();
}

double ComputeAverageEdgeLength(vtkPolyData *poly)
{
  double l = 0.0;
  size_t n = 0;
  
  vtkIdType nCells = poly->GetNumberOfCells();
  for(vtkIdType iCell = 0; iCell < nCells; iCell++)
    {
    // Get the points in this cell
    vtkIdType nPoints, *xPoints;
    poly->GetCellPoints(iCell, nPoints, xPoints);

    for(vtkIdType j = 0; j < nPoints; j++)
      {
      vtkIdType k = (j + 1) % nPoints;
      vnl_vector_fixed<vtkFloatingPointType,3> x1(poly->GetPoint(xPoints[j]));
      vnl_vector_fixed<vtkFloatingPointType,3> x2(poly->GetPoint(xPoints[k]));
      l += sqrt(dot_product(x1-x2,x1-x2));
      n++;
      }
    }

  return l / n;
}

void ComputeAreaElement(vtkPolyData *poly, vnl_vector<vtkFloatingPointType> &elt)
{
  // For each triangle in the polydata compute its area
  vtkIdType nCells = poly->GetNumberOfCells();
  
  // Initialize the area elt array
  elt.set_size(poly->GetNumberOfPoints());
  elt.fill(0.0);
  
  for(vtkIdType iCell = 0; iCell < nCells; iCell++)
    {
    // Get the points in this cell
    vtkIdType nPoints, *xPoints;
    poly->GetCellPoints(iCell, nPoints, xPoints);
    
    // Only triangles are admitted
    if(nPoints != 3)
      { 
      cerr << "Irregular face (n = " << nPoints << ") detected!" << endl;
      break;
      }

    // Get the three points
    vnl_vector_fixed<vtkFloatingPointType, 3> X0(poly->GetPoint(xPoints[0]));
    vnl_vector_fixed<vtkFloatingPointType, 3> X1(poly->GetPoint(xPoints[1]));
    vnl_vector_fixed<vtkFloatingPointType, 3> X2(poly->GetPoint(xPoints[2]));

    // Compute the area
    double xArea = TriangleArea(X0, X1, X2);
    if(xArea < 0)
      {
      cerr << "Negative area returned at cell " << iCell << endl;
      break;
      }

    // Add the area to all points
    elt[xPoints[0]] += xArea; elt[xPoints[1]] += xArea; elt[xPoints[2]] += xArea;
    }
}

int usage()
{
  cout << "Usage: " << endl;
  cout << "    cmrep_vskel [options] boundary.vtk output_skeleton.vtk" << endl;
  cout << "Parameters: " << endl;
  cout << "    boundary.vtk         Boundary mesh to skeletonize" << endl;
  cout << "    output.vtk           Where to output the skeleton" << endl;
  cout << "General Options:         " << endl;
  cout << "    -Q PATH              Path to the qvoronoi executable" << endl;
  cout << "Pruning Options:" << endl;
  cout << "    -e N                 Minimal number of mesh edges separating two generator" << endl;
  cout << "                         points of a VD face for it to be considered (try 2, 3)" << endl;
  cout << "    -p X.XX              Prune the mesh using factor X.XX (try 2.0). The " << endl;
  cout << "                         pruning algorithm deletes faces in the VD for " << endl;
  cout << "                         which the ratio of the geodesic distance between " << endl;
  cout << "                         the generating points and the euclidean distance " << endl;
  cout << "                         between these points is less than X.XX" << endl;
  cout << "    -c N                 Take at most N connected components of the skeleton" << endl;
  cout << "    -s mesh.vtk          Load a skeleton from mesh.vtk and compare to the output skeleton" << endl;
  cout << "    -g                   Compute full geodesic information. This is only useful for" << endl;
  cout << "                         debugging the pruning code." << endl;
  cout << "    -t                   Tolerance for the inside/outside search algorithm (default 1e-6)" << endl;
  cout << "                         Use lower values if holes appear in the skeleton" << endl;
  return -1;
}
  

int main(int argc, char *argv[])
{

  // Command line arguments
  string fnMesh, fnOutput, fnSkel, fnQVoronoi="qvoronoi";
  double xPrune = 2.0, xSearchTol = 1e-6;
  int nComp = 0, nDegrees = 0;
  bool flagGeodFull = false;

  // Check that there are at least three command line arguments
  if(argc < 3) return usage();
  fnMesh = argv[argc-2];
  fnOutput = argv[argc-1];
  
  // Parse the command line for options
  for(int iArg = 1; iArg < argc - 3; ++iArg)
    {
    string arg = argv[iArg];
    if(arg == "-p")
      {
      xPrune = atof(argv[++iArg]);
      }
    else if(arg == "-t")
      {
      xSearchTol = atof(argv[++iArg]);
      }
    else if(arg == "-c")
      {
      nComp = atoi(argv[++iArg]);
      }
    else if(arg == "-Q")
      {
      fnQVoronoi = argv[++iArg];
      }
    else if(arg == "-e")
      {
      nDegrees = atoi(argv[++iArg]);
      }
    else if(arg == "-s")
      {
      fnSkel = argv[++iArg];
      }
    else if(arg == "-g")
      {
      flagGeodFull = true;
      }
    else
      {
      cerr << "Bad option " << arg << endl;
      return -1;
      }
    }

  // Load the input mesh
  vtkPolyData *bndraw = ReadVTKData(fnMesh);
  bndraw->BuildLinks();
  bndraw->BuildCells();

  // The raw boundary must be triangulated and cleaned
  vtkTriangleFilter *fTriangle = vtkTriangleFilter::New();
  fTriangle->SetInput(bndraw);
  vtkCleanPolyData *fClean = vtkCleanPolyData::New();
  fClean->SetInput(fTriangle->GetOutput());
  fClean->Update();
  vtkPolyData *bnd = fClean->GetOutput();

  // Create a bounding box object
  bnd->GetPoints()->ComputeBounds();
  double *bbBnd = bnd->GetPoints()->GetBounds();
  printf("Bounding Box : %f %f %f %f %f %f\n", bbBnd[0], bbBnd[1], bbBnd[2], bbBnd[3], bbBnd[5], bbBnd[6]);
  vtkBoundingBox fBoundBox;
  fBoundBox.SetBounds(bbBnd);

  // Create a temporary file where to store the points
  char *fnPoints = tmpnam(NULL);
  cout << "Storing mesh point coordinates in " << fnPoints << endl;
  FILE *f = fopen(fnPoints, "wt");
  fprintf(f, "%d\n3\n", (int) bnd->GetNumberOfPoints());
  for(vtkIdType i = 0; i < bnd->GetNumberOfPoints(); i++)
    fprintf(f, "%f %f %f\n", bnd->GetPoint(i)[0], bnd->GetPoint(i)[1], bnd->GetPoint(i)[2]);
  fclose(f);

  // Call Qvoronoi 
  // qvoronoi p Fv < $WORK/$ID.boundary_points.txt > $WORK/$ID.voronoi.txt
  string fnVoronoiOutput = string(fnPoints) + "_voronoi.txt";
  char command[1024];
  sprintf(command, "%s p Fv < %s > %s", fnQVoronoi.c_str(), fnPoints, fnVoronoiOutput.c_str());
  cout << "Executing system command \"" << command << "\"" << endl;
  if(system(command) < 0)
    {
    cerr << "Call to QVoronoi failed" << endl;
    return -1;
    }

  // Array of generating points for each cell
  VertexPairArray pgen;

  // Load the file
  ifstream fin(fnVoronoiOutput.c_str());

  // Load the numbers
  size_t nv, np, junk;
  
  // First two lines
  fin >> junk;
  fin >> nv; 

  vtkSelectEnclosedPoints *sel = vtkSelectEnclosedPoints::New();
  sel->SetTolerance(xSearchTol);
  sel->Initialize(bnd);

  // Create an array of points
  vtkPoints *pts = vtkPoints::New();
  pts->SetNumberOfPoints(nv);

  // Create an array of in/out flags
  bool *ptin = new bool[nv];

  // Progress bar
  cout << "Selecting points inside mesh (n = " << nv << ")" << endl;
  cout << "|         |         |         |         |         |" << endl;
  size_t next_prog_mark = nv / 50;

  for(size_t i = 0; i < nv; i++)
    {
    double x,y,z;
    fin >> x;
    fin >> y;
    fin >> z;
    pts->SetPoint(i,x,y,z);

    // Is this point outside of the bounding box
    ptin[i] = fBoundBox.ContainsPoint(x,y,z) && sel->IsInsideSurface(x,y,z);
    
    if(i >= next_prog_mark)
      {
      cout << "." << flush;
      next_prog_mark += nv / 50;
      }
    }
  cout << "." << endl;

  // Read the number of cells
  fin >> np;

  // Progress bar
  cout << "Selecting faces using pruning criteria (n = " << np << ")" << endl;
  cout << "|         |         |         |         |         |" << endl;
  next_prog_mark = np / 50;

  // Create and configure Dijkstra's alg for geodesic distance
  VTKMeshHalfEdgeWrapper hewrap_geo(bnd);
  EuclideanDistanceMeshEdgeWeightFunction wfunc_geo;
  VTKMeshShortestDistance dijkstra_geo;
  dijkstra_geo.SetInputMesh(&hewrap_geo);
  dijkstra_geo.SetEdgeWeightFunction(&wfunc_geo);
  dijkstra_geo.ComputeGraph();

  // Create and configure Dijkstra's alg for edge counting
  VTKMeshHalfEdgeWrapper hewrap_edge(bnd);
  UnitLengthMeshEdgeWeightFunction wfunc_edge;
  VTKMeshShortestDistance dijkstra_edge;
  dijkstra_edge.SetInputMesh(&hewrap_edge);
  dijkstra_edge.SetEdgeWeightFunction(&wfunc_edge);
  dijkstra_edge.ComputeGraph();

  // Keep track of number pruned
  size_t npruned_geo=0, npruned_edge=0;

  // Create the polygons 
  vtkCellArray *cells = vtkCellArray::New();
  
  // Allocate the output radius data array
  vtkDoubleArray *daRad = vtkDoubleArray::New();
  daRad->SetNumberOfComponents(1);
  daRad->SetName("Radius");

  // Another array for prune strength
  vtkDoubleArray *daPrune = vtkDoubleArray::New();
  daPrune->SetNumberOfComponents(1);
  daPrune->SetName("Pruning Ratio");

  // Another array for prune strength
  vtkDoubleArray *daGeod = vtkDoubleArray::New();
  daGeod->SetNumberOfComponents(1);
  daGeod->SetName("Geodesic");

  for(size_t j = 0; j < np; j++)
    {
    bool isinf = false;
    bool isout = false;
    
    size_t m; fin >> m; m-=2;
    vtkIdType ip1, ip2; fin >> ip1; fin >> ip2;

    vtkIdType *ids = new vtkIdType[m];
    for(size_t k = 0; k < m; k++)
      {
      fin >> ids[k];

      // Is this point at infinity?
      if(ids[k] == 0) isinf = true; else ids[k]--;
      if(!ptin[ids[k]]) isout = true;
      }

    if(!isinf && !isout)
      {
      bool pruned = false;
      double r = 0, dgeo = 0;

      // Get the edge distance between generators
      dijkstra_edge.ComputeDistances(ip1, nDegrees);
      double elen = dijkstra_edge.GetVertexDistance(ip2);
      if(elen < nDegrees)
        { pruned = true; npruned_edge++; }
      else
        {
        // Get the Euclidean distance between generator points
        vnl_vector_fixed<vtkFloatingPointType,3> p1(bnd->GetPoint(ip1)); 
        vnl_vector_fixed<vtkFloatingPointType,3> p2(bnd->GetPoint(ip2)); 
        r = (p1 - p2).magnitude();
        
        // The geodesic distance between generators should exceed d * xPrune;
        dijkstra_geo.ComputeDistances(ip1, r * xPrune + 1);

        // Get the distance
        dgeo = dijkstra_geo.GetVertexDistance(ip2);

        // If the geodesic is too short, don't insert point
        if(dgeo < r * xPrune)
          { pruned = true; npruned_geo++; }
        }

      if(!pruned)
        {
        // add the cell
        cells->InsertNextCell(m, ids);
        daRad->InsertNextTuple(&r);
        daGeod->InsertNextTuple(&dgeo);
        double ratio = dgeo / r;
        daPrune->InsertNextTuple(&ratio);

        // add the pair of generators
        pgen.push_back(make_pair(ip1, ip2)); // TODO: is this numbering 0-based?
        }
      }

    if(j >= next_prog_mark)
      {
      cout << "." << flush;
      next_prog_mark += np / 50;
      }

    delete ids;
    }

  cout << "." << endl;
  cout << "Edge contraint pruned " << npruned_edge << " faces." << endl;
  cout << "Geodesic to Euclidean distance ratio contraint pruned " << npruned_geo << " faces." << endl;

  // Create the vtk poly data
  vtkPolyData *skel = vtkPolyData::New();
  skel->SetPoints(pts);
  skel->SetPolys(cells);
  skel->GetCellData()->AddArray(daRad);
  skel->GetCellData()->AddArray(daGeod);
  skel->GetCellData()->AddArray(daPrune);
  skel->BuildCells();
  skel->BuildLinks();

  // Drop the singleton points from the diagram (points not in any cell)
  vtkCleanPolyData *fltClean = vtkCleanPolyData::New();
  fltClean->SetInput(skel);
  fltClean->Update();

  // The output from the next branch
  vtkPolyData *polySave = fltClean->GetOutput();

  // Compute the connected components
  if(nComp > 0)
    {
    vtkPolyDataConnectivityFilter *fltConnect = vtkPolyDataConnectivityFilter::New();
    fltConnect->SetInput(fltClean->GetOutput());

    if(nComp == 1)
      fltConnect->SetExtractionModeToLargestRegion();
    else 
      {
      fltConnect->SetExtractionModeToSpecifiedRegions();
      fltConnect->InitializeSpecifiedRegionList();
      for(int rr = 0; rr < nComp; rr++)
        fltConnect->AddSpecifiedRegion(rr);
      }

    fltConnect->ScalarConnectivityOff();
    fltConnect->Update();
    polySave = fltConnect->GetOutput();
    cout << "Connected component constraint pruned " 
      << skel->GetNumberOfCells() - polySave->GetNumberOfCells() 
      << " faces." << endl; 
    }

  // Convert the cell data to point data
  vtkCellDataToPointData *c2p = vtkCellDataToPointData::New();
  c2p->SetInput(polySave);
  c2p->PassCellDataOn();
  c2p->Update();
  vtkPolyData *final = c2p->GetPolyDataOutput();

  // Compute mean, median thickness?
  double int_area = 0, int_thick = 0;
  vtkDataArray *finalRad = final->GetCellData()->GetArray("Radius");
  for(int i = 0; i < final->GetNumberOfCells(); i++)
    {
    double r = finalRad->GetTuple1(i);
    vtkCell *c = final->GetCell(i);
    vnl_vector_fixed<vtkFloatingPointType, 3> p1(final->GetPoint(c->GetPointId(0)));
    vnl_vector_fixed<vtkFloatingPointType, 3> p2(final->GetPoint(c->GetPointId(1)));
    vnl_vector_fixed<vtkFloatingPointType, 3> p3(final->GetPoint(c->GetPointId(2)));
    double a = fabs(TriangleArea(p1, p2, p3));
    int_area += a;
    int_thick += r * a;
    }
  cout << "Surface area: " << int_area << endl;
  cout << "Mean thickness: " << int_thick / int_area << endl;
    
  WriteVTKData(c2p->GetPolyDataOutput(), fnOutput);
}

