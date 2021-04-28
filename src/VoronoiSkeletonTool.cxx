#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <SparseMatrix.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkThresholdPoints.h>
#include <vtkBoundingBox.h>
#include <vtkCellArray.h>
#include <vtkCellDataToPointData.h>
#include <vtkPolyData.h>
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
#include <vtkPolyDataWriter.h>
#include <vtkMath.h>
#include <vtkTriangle.h>
#include <VTKMeshShortestDistance.h>
#include <VTKMeshHalfEdgeWrapper.h>
#include <vtkQuadricClustering.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_cross.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

#include <itkOrientedRASImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkVectorImage.h>

#include "util/ReadWriteVTK.h"

using namespace std;

typedef std::pair<vtkIdType, vtkIdType> VertexPair;
typedef std::set< std::pair<vtkIdType, vtkIdType> > VertexPairSet;
typedef std::vector<VertexPair> VertexPairArray;

void WriteVTKData(vtkUnstructuredGrid *data, string fn)
{
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fn.c_str());
  writer->SetInputData(data);
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

    delete[] ids;
    }

  // Create the vtk poly data
  vtkPolyData *poly = vtkPolyData::New();
  poly->SetPoints(pts);
  poly->SetPolys(cells);
  return poly;
}

inline double TriangleArea(
  const vnl_vector_fixed<double,3> &A, 
  const vnl_vector_fixed<double,3> &B, 
  const vnl_vector_fixed<double,3> &C)
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
      vnl_vector_fixed<double,3> x1(poly->GetPoint(xPoints[j]));
      vnl_vector_fixed<double,3> x2(poly->GetPoint(xPoints[k]));
      l += sqrt(dot_product(x1-x2,x1-x2));
      n++;
      }
    }

  return l / n;
}

void ComputeAreaElement(vtkPolyData *poly, vnl_vector<double> &elt)
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
    vnl_vector_fixed<double, 3> X0(poly->GetPoint(xPoints[0]));
    vnl_vector_fixed<double, 3> X1(poly->GetPoint(xPoints[1]));
    vnl_vector_fixed<double, 3> X2(poly->GetPoint(xPoints[2]));

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
  cout << "    -z <level> <mode>    Subdivide input mesh prior to skeletonization" << endl;
  cout << "                         mode is either 'loop' or 'linear'" << endl;
  cout << "Pruning Options:" << endl;
  cout << "    -e N                 Minimal number of mesh edges separating two generator" << endl;
  cout << "                         points of a VD face for it to be considered (try 2, 3)" << endl;
  cout << "    -p X.XX              Prune the mesh using factor X.XX (try 2.0). The " << endl;
  cout << "                         pruning algorithm deletes faces in the VD for " << endl;
  cout << "                         which the ratio of the geodesic distance between " << endl;
  cout << "                         the generating points and the euclidean distance " << endl;
  cout << "                         between these points is less than X.XX" << endl;
  cout << "    -c N                 Take at most N connected components of the skeleton" << endl;
  cout << "    -g                   Compute full geodesic information. This is only useful for" << endl;
  cout << "                         debugging the pruning code." << endl;
  cout << "    -t                   Tolerance for the inside/outside search algorithm (default 1e-6)" << endl;
  cout << "                         Use lower values if holes appear in the skeleton. Set to zero to" << endl;
  cout << "                         disable pruning of outside vertices" << endl;
  cout << "Output Options: " << endl;;
  cout << "    -s mesh.vtk          Load a skeleton from mesh.vtk and compare to the output skeleton" << endl;
  cout << "    -R N xyz.mat d.mat   Generate N random samples from the skeleton and save their coordiantes" << endl;
  cout << "                         to xyz.mat and geodesic distances to d.mat" << endl;
  cout << "    -T name.vtk          Generate thickness map on the boundary. The thickness is the distance" << endl;
  cout << "                         from each boundary point to the closest pruned skeleton point" << endl;
  cout << "    -I in.nii thickness.nii depth.nii " << endl; 
  cout << "                         Generate thickness map in an image. Input is a binary image. " << endl;
  cout << "                         Output 1 is a thickness image; Output 2 is a depth map" << endl;
  cout << "    -q n_bins            Postprocess skeleton with VTK's quadric clustering filter" << endl;
  cout << "                         The effect is to reduce the number of vertices in the skeleton" << endl;
  cout << "                         Parameter n_bins is the number of bins in each dimension" << endl;
  cout << "                         A good value for n_bins is 20-50" << endl;
  cout << "    -d mesh.vtk          Generate a Delaunay tetrahedralization of the input point set, with" << endl;
  cout << "                         the pruned parts of the skeleton excluded. Use with tetfill to generate" << endl;
  cout << "                         a thickness map in image space (different from -I output, this is " << endl;
  cout << "                         distance from the skeleton vertices to the boundary generator points" << endl;
  cout << "Other Options: " << endl;
  cout << "    -S image array mode  Sample from image 'image' and store as array 'array'" << endl;
  cout << "                         mode is one of 'mean', 'max'. Must be used together with -T command" << endl;
  return -1;
}
  
// Sampling info
enum SamplingMode {
  MEAN, MAX 
};

struct SamplingStruct {
  string fnImage, array;
  SamplingMode mode;
};

enum SubMode { LINEAR = 0, LOOP };

int main(int argc, char *argv[])
{

  // Command line arguments
  string fnMesh, fnOutput, fnSkel, fnQVoronoi="qvoronoi", fnOutThickness;
  string fnImgRef, fnImgThick, fnImgDepth;
  string fnTetraMesh;
  double xPrune = 2.0, xSearchTol = 1e-6;
  int nComp = 0, nDegrees = 0, nRandSamp = 0, nBins = 0;
  bool flagGeodFull = false;
  int subLevel = 0;
  SubMode subMode = LINEAR;

  // Sampling stuff
  vector<SamplingStruct> sampling;

  // Check that there are at least three command line arguments
  if(argc < 3) return usage();
  fnMesh = argv[argc-2];
  fnOutput = argv[argc-1];
  string fnXYZ, fnDist;
  
  // Parse the command line for options
  for(int iArg = 1; iArg < argc - 3; ++iArg)
    {
    string arg = argv[iArg];
    if(arg == "-p")
      {
      xPrune = atof(argv[++iArg]);
      }
    else if(arg == "-z")
      {
      subLevel = atoi(argv[++iArg]);
      std::string mode = argv[++iArg];
      subMode = mode == "loop" ? LOOP : LINEAR;
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
    else if(arg == "-T")
      {
      fnOutThickness = argv[++iArg];
      }
    else if(arg == "-q")
      {
      nBins = atoi(argv[++iArg]);
      }
    else if(arg == "-I")
      {
      fnImgRef = argv[++iArg];
      fnImgThick = argv[++iArg];
      fnImgDepth = argv[++iArg];
      }
    else if(arg == "-S")
      {
      SamplingStruct ss;
      ss.fnImage = argv[++iArg];
      ss.array = argv[++iArg];
      string mode = argv[++iArg];
      if(mode == "mean") ss.mode = MEAN;
      else if(mode == "max") ss.mode = MAX;
      else
        cerr << "Bad mode for -S: " << ss.mode << endl;
      sampling.push_back(ss);
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
    else if(arg == "-R")
      {
      nRandSamp = atoi(argv[++iArg]);
      fnXYZ = argv[++iArg];
      fnDist = argv[++iArg];
      }
    else if(arg == "-d")
      {
      fnTetraMesh = argv[++iArg];
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
  fTriangle->SetInputData(bndraw);
  vtkCleanPolyData *fClean = vtkCleanPolyData::New();
  fClean->SetInputConnection(fTriangle->GetOutputPort());
  fClean->SetTolerance(1e-4);
  fClean->Update();
  vtkPolyData *bnd = fClean->GetOutput();

  // Do we want subdivision?
  if(subLevel > 0)
    {
    if(subMode == LINEAR)
      {
      vtkLinearSubdivisionFilter *fSub = vtkLinearSubdivisionFilter::New();
      fSub->SetInputConnection(fClean->GetOutputPort());
      fSub->SetNumberOfSubdivisions(subLevel);
      fSub->Update();
      bnd = fSub->GetOutput();
      }
    else
      {
      vtkLoopSubdivisionFilter *fSub = vtkLoopSubdivisionFilter::New();
      fSub->SetInputConnection(fClean->GetOutputPort());
      fSub->SetNumberOfSubdivisions(subLevel);
      fSub->Update();
      bnd = fSub->GetOutput();
      }
    }

  // Create a bounding box object
  bnd->GetPoints()->ComputeBounds();
  double *bbBnd = bnd->GetPoints()->GetBounds();
  printf("Bounding Box : %f %f %f %f %f %f\n", bbBnd[0], bbBnd[1], bbBnd[2], bbBnd[3], bbBnd[5], bbBnd[6]);
  vtkBoundingBox fBoundBox;
  fBoundBox.SetBounds(bbBnd);

  // Create a temporary file where to store the points
#ifndef WIN32
  char fnTemplate[] = "/tmp/voronoi_pts.XXXXXX";
  char *fnPoints = mktemp(fnTemplate);
#else
  char *fnPoints = tmpnam(NULL);
#endif
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

  // Create a point array to hold radius information. The radius of the MIB
  // is computed as distance from vertices of Voronoi cells to generating
  // points.
  vtkDoubleArray *daPointRadius = vtkDoubleArray::New();
  daPointRadius->SetNumberOfComponents(1);
  daPointRadius->SetNumberOfTuples(nv);
  daPointRadius->SetName("VoronoiRadius");
  daPointRadius->FillComponent(0, NAN);

  // Set of tetrahedra corresponding to the Voronoi triangulation
  typedef std::set<size_t> Hedron;
  std::vector<Hedron> hedra(nv);

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
    if(xSearchTol > 0)
      ptin[i] = fBoundBox.ContainsPoint(x,y,z) && sel->IsInsideSurface(x,y,z);
    else
      ptin[i] = fBoundBox.ContainsPoint(x,y,z);
    
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
        vnl_vector_fixed<double,3> p1(bnd->GetPoint(ip1)); 
        vnl_vector_fixed<double,3> p2(bnd->GetPoint(ip2)); 
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


        // For each vertex of the cell, add the generating points to its (tetra)hedra
        for(size_t k = 0; k < m; k++)
          {
          hedra[ids[k]].insert(ip1);
          hedra[ids[k]].insert(ip2);
          }

        // Compute the radius at each vertex of the cell
        for(size_t k = 0; k < m; k++)
          {
          size_t v = ids[k];
          vnl_vector_fixed<double, 3> Vk(pts->GetPoint(v));
          vnl_vector_fixed<double, 3> G(bnd->GetPoint(ip1));
          double r = (Vk - G).magnitude();
          double r_old = daPointRadius->GetComponent(v, 0);
          if(::isnan(r_old))
            daPointRadius->SetComponent(v, 0, r);
          }
        }
      }

    if(j >= next_prog_mark)
      {
      cout << "." << flush;
      next_prog_mark += np / 50;
      }

    delete[] ids;
    }

  cout << "." << endl;
  cout << "Edge contraint pruned " << npruned_edge << " faces." << endl;
  cout << "Geodesic to Euclidean distance ratio contraint (" << xPrune << ") pruned " << npruned_geo << " faces." << endl;

  // Clean up files
  if (fin.is_open()) 
    {
    fin.close();
    remove(fnVoronoiOutput.c_str());
    } 
  remove(fnPoints);

  // Did we get tetrahedra?
  if(fnTetraMesh.size()) 
    {
    vtkUnstructuredGrid *testtet = vtkUnstructuredGrid::New();
    testtet->SetPoints(bnd->GetPoints());

    vtkDoubleArray *tetrad = vtkDoubleArray::New();
    tetrad->SetNumberOfComponents(1);
    tetrad->SetName("VoronoiRadius");
    testtet->GetCellData()->AddArray(tetrad);

    vtkDoubleArray *tetra_ctr = vtkDoubleArray::New();
    tetra_ctr->SetNumberOfComponents(3);
    tetra_ctr->SetName("VoronoiCenter");
    testtet->GetCellData()->AddArray(tetra_ctr);

    for(unsigned int i = 0; i < nv; i++)
      {
      if(hedra[i].size() == 4)
        {
        vtkIdType tetids[4];
        Hedron::const_iterator hit = hedra[i].begin();
        for(unsigned int j = 0; j < 4; j++, hit++)
          tetids[j] = *hit;
        testtet->InsertNextCell(VTK_TETRA, 4, tetids);
        tetrad->InsertNextTuple1(daPointRadius->GetComponent(i, 0));

        double *ctr = pts->GetPoint(i);
        tetra_ctr->InsertNextTuple3(ctr[0], ctr[1], ctr[2]);
        }
      }

    // Write the tetrahedra
    WriteVTKData(testtet, fnTetraMesh.c_str());
    }

  // Create the vtk poly data
  vtkPolyData *skel = vtkPolyData::New();
  skel->SetPoints(pts);
  skel->SetPolys(cells);
  skel->GetCellData()->AddArray(daRad);
  skel->GetCellData()->AddArray(daGeod);
  skel->GetCellData()->AddArray(daPrune);
  skel->GetPointData()->AddArray(daPointRadius);
  skel->BuildCells();
  skel->BuildLinks();

  // Drop the singleton points from the diagram (points not in any cell)
  vtkCleanPolyData *fltClean = vtkCleanPolyData::New();
  fltClean->SetInputData(skel);
  fltClean->Update();
  cout << "Clean filter: trimmed " 
    << skel->GetNumberOfPoints() << " vertices to "
    << fltClean->GetOutput()->GetNumberOfPoints() << endl;


  // The output from the next branch
  vtkPolyData *polySave = fltClean->GetOutput();

  // Compute the connected components
  if(nComp > 0)
    {
    vtkPolyDataConnectivityFilter *fltConnect = vtkPolyDataConnectivityFilter::New();
    fltConnect->SetInputConnection(fltClean->GetOutputPort());

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

    // Don't see why this is necessary, but Connectivity filter does not remove points,
    // just faces. So we need another clear filter
    vtkCleanPolyData *fltWhy = vtkCleanPolyData::New();
    fltWhy->SetInputConnection(fltConnect->GetOutputPort());
    fltWhy->Update();

    cout << "Connected component constraint pruned " 
      << polySave->GetNumberOfCells() - fltWhy->GetOutput()->GetNumberOfCells() << " faces and "
      << polySave->GetNumberOfPoints() - fltWhy->GetOutput()->GetNumberOfPoints() << " points." << endl; 
    polySave = fltWhy->GetOutput();
    }

  // Convert the cell data to point data
  vtkCellDataToPointData *c2p = vtkCellDataToPointData::New();
  c2p->SetInputData(polySave);
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
    if(c->GetNumberOfPoints() == 3)
      {
      vnl_vector_fixed<double, 3> p1(final->GetPoint(c->GetPointId(0)));
      vnl_vector_fixed<double, 3> p2(final->GetPoint(c->GetPointId(1)));
      vnl_vector_fixed<double, 3> p3(final->GetPoint(c->GetPointId(2)));
      double a = fabs(TriangleArea(p1, p2, p3));
      int_area += a;
      int_thick += r * a;
      }
    }
  cout << "Surface area: " << int_area << endl;
  cout << "Mean thickness: " << int_thick / int_area << endl;
    
  vtkPolyData *skelfinal = c2p->GetPolyDataOutput();

  // Quadric clustering
  if(nBins > 0)
    {
    // Calculate appropriate bin size
    double *bbBnd = skelfinal->GetPoints()->GetBounds();
    vtkBoundingBox fbb;
    fbb.SetBounds(bbBnd);
    double binsize = fbb.GetMaxLength() / nBins;

    vtkQuadricClustering *fCluster = vtkQuadricClustering::New();
    fCluster->SetNumberOfDivisions(
      ceil(fbb.GetLength(0) / binsize), 
      ceil(fbb.GetLength(1) / binsize),
      ceil(fbb.GetLength(2) / binsize));
    fCluster->SetInputData(c2p->GetPolyDataOutput());
    fCluster->SetCopyCellData(1);
    fCluster->Update();

    printf("QuadClustering (%d x %d x %d blocks) :\n", 
      fCluster->GetNumberOfXDivisions(),
      fCluster->GetNumberOfYDivisions(),
      fCluster->GetNumberOfZDivisions());
    printf("  Input mesh: %d points, %d cells\n", 
      (int) skelfinal->GetNumberOfPoints(),  
      (int) skelfinal->GetNumberOfCells());
    printf("  Output mesh: %d points, %d cells\n", 
      (int) fCluster->GetOutput()->GetNumberOfPoints(),  
      (int) fCluster->GetOutput()->GetNumberOfCells());

    // Convert cell data to point data again
    vtkCellDataToPointData *c2p = vtkCellDataToPointData::New();
    c2p->SetInputConnection(fCluster->GetOutputPort());
    c2p->PassCellDataOn();
    c2p->Update();
    skelfinal = c2p->GetPolyDataOutput();
    }

  // Write the skeleton out
  WriteVTKData(skelfinal, fnOutput);

  // Generate thickness map
  if(fnOutThickness.length())
    {
    // Initialize thickness array
    vtkDoubleArray *daRad = vtkDoubleArray::New();
    daRad->SetNumberOfComponents(1);
    daRad->SetName("Thickness");

    // Create locator for finding closest points
    vtkCellLocator *loc = vtkCellLocator::New();
    loc->SetDataSet(skelfinal);
    loc->CacheCellBoundsOn();
    loc->BuildLocator();

    // Vector of sampling images and arrays
    typedef itk::VectorImage<double, 3> SampleImageType;
    vector<SampleImageType::Pointer> imgSam(sampling.size(), nullptr);
    vector<vtkDoubleArray *> daSam(sampling.size(), NULL);

    // Interpolators
    typedef itk::LinearInterpolateImageFunction<SampleImageType, double> SampleInterpType;
    vector<SampleInterpType::Pointer> interp(sampling.size(), nullptr);

    // Vector of step sizes for each image
    vector<double> stepsize(sampling.size(), 0);

    // Create and add all the sampling arrays 
    for(int i = 0; i < sampling.size(); i++)
      {
      // Load the image that needs to be sampled
      typedef itk::ImageFileReader<SampleImageType> SampleReaderType;
      SampleReaderType::Pointer reader = SampleReaderType::New();
      reader->SetFileName(sampling[i].fnImage.c_str());
      reader->Update();
      imgSam[i] = reader->GetOutput(); 

      // Interpolator
      interp[i] = SampleInterpType::New();
      interp[i]->SetInputImage(imgSam[i]);

      // Create an output array
      int ncomp = imgSam[i]->GetNumberOfComponentsPerPixel();
      daSam[i] = vtkDoubleArray::New();
      daSam[i]->SetNumberOfComponents(ncomp);
      daSam[i]->SetNumberOfTuples(bndraw->GetNumberOfPoints());
      daSam[i]->SetName(sampling[i].array.c_str());

      // Compute step size - smallest of spacing values
      double minspc = 
        std::min(
          std::min(imgSam[i]->GetSpacing()[0], imgSam[i]->GetSpacing()[1]),
          imgSam[i]->GetSpacing()[2]);

      // Arbitrary - set to 1/4 voxel extent
      stepsize[i] = minspc / 4;

      // Add the array to the boundary mesh
      bndraw->GetPointData()->AddArray(daSam[i]);
      }

    // Compute thickness values
    for(size_t i = 0; i < (size_t) bndraw->GetNumberOfPoints(); i++)
      {
      double xs[3], d2, d;
      int subid;
      vtkIdType cellid;

      loc->FindClosestPoint(
        bndraw->GetPoint(i), xs, cellid, subid, d2);

      d = sqrt(d2);
      daRad->InsertNextTuple(&d);

      // Sample along the interval between point i and point xs
      for(int j = 0; j < imgSam.size(); j++)
        {
        int n_steps = std::max((int) (0.5 + d / stepsize[j]), 5);

        // Current sample
        SampleImageType::PixelType sample, sample_sum, sample_max;

        // Ncomp
        int ncomp = imgSam[j]->GetNumberOfComponentsPerPixel();

        for(int k = 0; k < n_steps; k++)
          {
          // Interpolate a point here
          double t = k / (n_steps - 1.0);
          itk::Point<double, 3> pt;
          pt[0] = bndraw->GetPoint(i)[0] * t + xs[0] * (1-t);
          pt[1] = bndraw->GetPoint(i)[1] * t + xs[1] * (1-t);
          pt[2] = bndraw->GetPoint(i)[2] * t + xs[2] * (1-t);

          // Map point to DICOM coordinates from NIFTI coordinates - hard coded
          pt[0] *= -1; pt[1] *= -1;

          // Map to voxel continuous index
          SampleInterpType::ContinuousIndexType cix;
          imgSam[j]->TransformPhysicalPointToContinuousIndex(pt, cix);
          sample = interp[j]->EvaluateAtContinuousIndex(cix);

          // If first point, store
          if(k == 0)
            {
            sample_sum = sample;
            sample_max = sample;
            }
          else
            {
            for(int c = 0; c < ncomp; c++)
              {
              sample_sum[c] += sample[c];
              sample_max[c] = std::max(sample_max[c], sample[c]);
              }
            }
          }

        // Store the information
        for(int c = 0; c < ncomp; c++)
          {
          if(sampling[j].mode == MEAN)
            {
            daSam[j]->SetComponent(i, c, sample_sum[c] / n_steps);
            }
          else
            {
            daSam[j]->SetComponent(i, c, sample_max[c]);
            }
          }
        }
      }

    // Add array to boundary data
    bndraw->GetPointData()->AddArray(daRad);

    // Write thickness map
    WriteVTKData(bndraw, fnOutThickness);
    } 

  if(fnImgRef.length() && fnImgThick.length())
    {
    // Read the reference image
    typedef itk::OrientedRASImage<short,3> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fnImgRef.c_str());
    try 
      {
      reader->Update();
      } 
    catch (itk::ExceptionObject &exc)
      {
      cerr << "Unable to read reference image: " << exc << endl;
      return -1;
      }
    ImageType::Pointer ref = reader->GetOutput();

    // Create the depth and thickness maps
    typedef itk::Image<float, 3> FloatImageType;
    FloatImageType::Pointer iout[2];
    for(size_t i = 0; i < 2; i++)
      {
      iout[i] = FloatImageType::New();
      iout[i]->SetRegions(ref->GetBufferedRegion());
      iout[i]->CopyInformation(ref);
      iout[i]->Allocate();
      iout[i]->FillBuffer(0.0);
      }
    
    // Create locator for finding closest points
    vtkCellLocator *loc = vtkCellLocator::New();
    loc->SetDataSet(skelfinal);
    loc->CacheCellBoundsOn();
    loc->BuildLocator();

    // Compute distance at each 1 voxel in the reference image
    typedef itk::ImageRegionConstIteratorWithIndex<ImageType> RefIterator;
    typedef itk::ImageRegionIterator<FloatImageType> OutIterator;
    RefIterator rit(ref, ref->GetBufferedRegion());
    OutIterator itthick(iout[0], ref->GetBufferedRegion());
    OutIterator itdepth(iout[1], ref->GetBufferedRegion());
    for(; !rit.IsAtEnd(); ++rit, ++itthick, ++itdepth)
      {
      if(rit.Get())
        {
        // Get the RAS coordinate of the voxel
        ImageType::PointType p;
        ref->TransformIndexToRASPhysicalPoint(rit.GetIndex(), p);
        double pt[3]; pt[0] = p[0]; pt[1] = p[1]; pt[2] = p[2];

        // Find closest point on the skeleton
        double xs[3], d2, d, t;
        int subid;
        vtkIdType cellid;
        
        loc->FindClosestPoint(pt, xs, cellid, subid, d2);
        d = sqrt(d2);

        // Get thickness at that cell
        t = finalRad->GetTuple1(cellid);

        // Set thickness, depth
        itthick.Set(t);
        itdepth.Set(d / t);
        /*
          {
          cout << "Index: " << rit.GetIndex() << endl;
          cout << "Point: " << p << endl;
          cout << "Distance: " << d << endl;
          cout << "Thickness: " << t << endl;
          cout << "Matching cell: " << cellid << endl;
          cout << "Corr.Point: " << xs[0] << " " << xs[1] << " " << xs[2] << endl;
          }
        */
        }
      }

    // Write output images
    typedef itk::ImageFileWriter<FloatImageType> WriterType;
    for(size_t i = 0; i < 2; i++)
      {
      WriterType::Pointer writer = WriterType::New();
      writer->SetInput(iout[i]);
      writer->SetFileName(i==0 ? fnImgThick : fnImgDepth);
      writer->Update();
      }
    }

  // Generate random samples
  if(nRandSamp > 0)
    {
    // Convert the mesh to triangles
    vtkTriangleFilter *fltTri = vtkTriangleFilter::New();
    fltTri->SetInputData(skelfinal);

    // Clean the mesh
    vtkCleanPolyData *fltClean = vtkCleanPolyData::New();
    fltClean->SetInputConnection(fltTri->GetOutputPort());
    fltClean->Update();
    vtkPolyData *xMesh = fltClean->GetOutput();
    xMesh->BuildCells();
    xMesh->BuildLinks();

    // Compute all the edges
    typedef ImmutableSparseArray<double> ImmutableArr;
    typedef ImmutableArr::VNLSourceType MutableArr;

    // Initialize the adjacency matrix
    size_t np = xMesh->GetNumberOfPoints();
    MutableArr M0(np, np);

    // Traverse all the cells in the VTK mesh, recording all the available
    // edges in the mesh. 
    for(unsigned int iCell = 0; iCell < (unsigned int) xMesh->GetNumberOfCells(); iCell++)
      {
      // Get the points for this cell
      vtkIdType nPoints, *xPoints;
      xMesh->GetCellPoints(iCell, nPoints, xPoints);

      // Walk around the list of points
      for(unsigned int j = 0; j < (unsigned int) nPoints; j++)
        {
        // Get the head and the tail of the current half-edge
        unsigned int iTail = xPoints[j], iHead = xPoints[(j+1) % nPoints];

        // Compute the distance
        vnl_vector<double> x1(xMesh->GetPoint(iTail), 3);
        vnl_vector<double> x2(xMesh->GetPoint(iHead), 3);
        double dist = (x1 - x2).magnitude();

        // Add the edge to the list
        M0(iTail,iHead) = dist;
        M0(iHead,iTail) = dist;
        }
      }

    // Build a sparse matrix from the edge data
    ImmutableArr M;
    M.SetFromVNL(M0); 

    // Generate a random subset of landmarks in the image
    size_t *xLandmarks = new size_t[np];
    for(int i = 0; i < (int) np; i++) xLandmarks[i] = i;
    random_shuffle(xLandmarks, xLandmarks+np);

    // Compute distances between landmarks
    typedef DijkstraShortestPath<double> Dijkstra;
    unsigned int *row_index = new unsigned int[M.GetNumberOfRows()];
    for(size_t q = 0; q < M.GetNumberOfRows(); q++)
      row_index[q] = M.GetRowIndex()[q];
    unsigned int *col_index = new unsigned int[M.GetNumberOfSparseValues()];
    for(size_t q = 0; q < M.GetNumberOfSparseValues(); q++)
      col_index[q] = M.GetColIndex()[q];

    Dijkstra dijk(np, row_index, col_index, M.GetSparseData());

    // Initialize the distance matrix
    vnl_matrix<double> mDist(nRandSamp, nRandSamp, 0.0);

    // Using the brute force approach
    for(int i = 0; i < (int) nRandSamp; i++)
      {
      // Compute all pairs shortest paths
      dijk.ComputePathsFromSource(xLandmarks[i]);
      const double *xDist = dijk.GetDistanceArray();

      // Get the landmark-to-landmark distances
      for(int j = 0; j < (int) nRandSamp; j++)
        {
        double d = xDist[xLandmarks[j]];
        mDist[i][j] = d * d;
        }

      cout << ".";
      if( ((i+1) % 64) == 0 || i + 1 == (int) nRandSamp )
        cout << " n = " << i+1 << endl;
      else
        cout << flush;
      }

    // Save the distance matrix
    ofstream of1(fnDist.c_str());
    of1 << mDist << endl;
    of1.close();
    // vnl_matlab_filewrite exporter(argv[3]);
    // exporter.write(mDist, "dist");

    // Save the coordinates of the landmarks
    vnl_matrix<double> mCoord(nRandSamp, 3, 0.0);
    for(int i = 0; i < (int) nRandSamp; i++)
      {
      double *xPoint = xMesh->GetPoint(xLandmarks[i]);
      mCoord[i][0] = xPoint[0];
      mCoord[i][1] = xPoint[1];
      mCoord[i][2] = xPoint[2];
      }

    // Save the distance matrix
    ofstream of2(fnXYZ.c_str());
    of2 << mCoord << endl;
    of2.close();
    // vnl_matlab_filewrite exporter2(argv[2]);
    // exporter2.write(mCoord, "xyz");

    delete[] xLandmarks;
    delete[] row_index;
    delete[] col_index;
    }
}

