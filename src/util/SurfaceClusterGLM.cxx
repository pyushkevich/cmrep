#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkClipDataSet.h>
#include <vtkThreshold.h>
#include <vtkConnectivityFilter.h>
#include "vtkPointData.h"
#include "ReadWriteVTK.h"
#include "vtkTriangleFilter.h"
#include "vtkBandedPolyDataContourFilter.h"
#include "vtkClipPolyData.h"
#include "vtkThresholdPoints.h"
#include "vtkPolyDataConnectivityFilter.h"
// #include "vtkUnstructuredGridToPolyDataFilter.h"
#include "vtkMergeCells.h"
#include "vtkGeometryFilter.h"
#include "vtkTriangle.h"
#include "vtkSmartPointer.h"
#include "vnl/vnl_file_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_rank.h"
#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>

#include <vtksys/SystemTools.hxx>


#include <exception>
#include <string>
#include <cstdarg>

extern "C" {
  double tnc ( double t, double df, double delta, int *ifault );
}

class MCException : public std::exception
{
public:
  MCException(const char *fmt, ...)
    {
    char buffer[1024];
    va_list parg;
    va_start(parg, fmt);
    vsprintf(buffer, fmt, parg);
    va_end(parg);
    message=buffer;
    }

  virtual ~MCException() throw() {}

  virtual const char *what() const throw()
    { return message.c_str(); }

private:
  std::string message;
};


using namespace std;

const char *usage_text = 
  "meshglm - a statistical tool for VTK meshes usage:\n"
  "  meshglm [switches]\n\n"
  "required switches:\n"
  "  -m / --mesh mesh_in.vtk mesh_out.vtk      \n"
  "                 Specify input and output meshes. The -m option can be\n"
  "                 repeated multiple times for different anatomical regions.\n"
  "                 Input mesh should provide the values of the dependent \n"
  "                 variable for each subject in a single attribute array.\n"
  "  -a / --array NAME                         \n"
  "                 Name of the data array that holds the dependent variable.\n"
  "                 The attribute array must have as many components as there\n"
  "                 are subjects in the GLM.\n"
  "  -g / --glm design.txt contrast.txt        \n"
  "                 GLM specification. The design matrix is k x m, where k is \n"
  "                 the number of subjects and m is the number of independent \n"
  "                 variables. The contrast vector is 1 x m.\n\n"
  "optional switches:\n"
  "  -c / --cell    \n"
  "                 Run GLM on cell data instead of point data\n"
  "  -p / --perm N  \n"
  "                 Perform permutation testing with N random permutations. \n"
  "                 This will generate a point-wise (or cell-wise) corrected\n"
  "                 p-value. \n"
  "  -s / --stat C|T|P                        \n"
  "                 The statistic on which permutation tests are performed:\n"
  "                 C: contrast; T: t-statistic; P: p-value (default: T).\n"
  "  -t / --threshold T                       \n"
  "                 The threshold to define clusters. If omitted, no cluster\n"
  "                 analysis will be performed. The threshold applies to the\n"
  "                 statistic specified by the -s flag. It must be positive.\n"
  "  -d / --diffusion T \n"
  "                 Amount of diffusion (spatial smoothing) to apply to data\n"
  "                 before statistical analysis. A simple implementation of \n"
  "                 the heat equation is used. Very roughly, this corresponds\n"
  "                 to smoothing the data with a Gaussian with standard \n"
  "                 deviation sqrt(2T)*L, where L is the average distance \n"
  "                 between neighboring vertices.\n"
  "  --delta-t DT\n"
  "                 Specify a different step size for diffusion (default: 0.01).\n"
  "  -e / --edges   \n"
  "                 Generate separate output files containing cluster edges.\n"
  "                 The file names are derived from the -m output parameters.\n"
  "  -M / --missing \n"
  "                 Input data (per-vertex observations) have missing data encoded\n"
  "                 as NaNs. The algorithm will exclude these observations from the\n"
  "                 GLM and compute the GLM with only non-missing observations\n"
  "                 NOTE: on some operating systems NaNs can only be read from VTK\n"
  "                 meshes in binary (not ASCII) format.\n"
  "  -B / --binary  Write output meshes in binary format. This is turned on automatically\n"
  "                 when the -M flag is supplied\n"
  "  -f / --flane   Use the Freedman-Lane procedure for permutation testing\n"
  "  -n / --nuissance nuissance.txt\n"
  "                 Optional Freedman-Lane specification. The nuissance vector is 1 x m, where\n"
  "                 m is the number of independent variables. Nuissance variables should \n"
  "                 be set to 1. Default specification assumes 0 variables in the \n"
  "                 contrast vector are nuissance variables. \n"
  "  -X <array>     Exclusion array. Should have the same dimensionality as the array \n"
  "                 passed to -a. For every point/subject where the exclusion array is \n"
  "                 non-zero, the corresponding value in the main array will be set to NaN\n"
  "  -T             Apply triangulation filter to input meshes (polydata only). This is\n"
  "                 in cases when input data has non-triangle faces\n"
  "  -R rowmask.txt Row mask. The file should contain as many numbers as there are rows in\n"
  "                 the design matrix, with 1 meaning the subject should be included in the \n"
  "                 analysis, 0 meaning excluded. This makes it easier to run analysis on \n"
  "                 different subsets of subjects without generating separate input meshes \n";

int usage()
{
  cout << usage_text; 
  return -1;
}

struct Cluster
{
  // Surface area of the cluster
  double area;

  // Power of the thresholded data over the cluster
  double power;

  // Cluster size in nodes (?)
  size_t n;

  // p-values
  double pArea, pPower;

  // t-value
  double tvalue;
 
  // Dummy constructor
  Cluster() : area(0.0), power(0.0), n(0), tvalue(0.0) {};
};

typedef std::vector<Cluster> ClusterArray;

// Templated IO classes to read mesh data
template <class TMeshType>
TMeshType * ReadMesh(const char *fname)
{ return NULL; }

template <>
vtkUnstructuredGrid *ReadMesh<>(const char *fname)
{
  vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}

template <>
vtkPolyData *ReadMesh<>(const char *fname)
{
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}

template <class TMeshType>
void WriteMesh(TMeshType *mesh, const char *fname, bool force_binary = false)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname, bool force_binary)
{
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(force_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname, bool force_binary)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(force_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

enum Domain { POINT, CELL };

bool operator < (const vnl_vector_fixed<double,3> &a, const vnl_vector_fixed<double, 3> &b)
{
  return(
    (a[0] < b[0]) ||
    (a[0] == b[0] && a[1] < b[1]) ||
    (a[0] == b[0] && a[1] == b[1] && a[2] < b[2])); 
}

vtkFloatArray *
AddArrayToMesh(
  vtkDataSet *t, 
  Domain dom, 
  const char *name, 
  unsigned int n_comp, 
  float defval, 
  bool make_active = false)
{
  vtkFloatArray * array = vtkFloatArray::New();
  array->SetNumberOfComponents(n_comp);
  array->SetName(name);

  if(dom == POINT)
    {
    array->SetNumberOfTuples(t->GetNumberOfPoints());
    t->GetPointData()->AddArray(array);
    if(make_active)
      t->GetPointData()->SetActiveScalars(name);
    }
  else
    {
    array->SetNumberOfTuples(t->GetNumberOfCells());
    t->GetCellData()->AddArray(array);
    if(make_active)
      t->GetCellData()->SetActiveScalars(name);
    }

  for(unsigned int k = 0; k < n_comp; k++)
    array->FillComponent(k, defval);

  return array;
}

vtkFloatArray *
AddArrayToMesh(vtkDataSet *t, Domain dom, const string &name, 
  unsigned int n_comp, float defval, bool make_active = false)
{ return AddArrayToMesh(t, dom, name.c_str(), n_comp, defval, make_active); }

vtkDataArray * GetArrayFromMesh(vtkDataSet *t, Domain dom, const char *name)
{
  vtkDataArray * data;
  if(dom == POINT)
    data = t->GetPointData()->GetArray(name);
  else
    data = t->GetCellData()->GetArray(name);
  return data;
}
vtkDataArray * GetArrayFromMesh(vtkDataSet *t, Domain dom, const string &name)
{ return GetArrayFromMesh(t, dom, name.c_str()); }

vtkDataArray * GetScalarsFromMesh(vtkDataSet *t, Domain dom)
{
  vtkDataArray * data;
  if(dom == POINT)
    data = t->GetPointData()->GetScalars();
  else
    data = t->GetCellData()->GetScalars();
  return data;
}

template <class TMeshType>
void
MergeMeshPieces(vtkDataSet *p1, vtkDataSet *p2, TMeshType *out)
{
  // Create a hashtable for coordinates
  typedef vnl_vector_fixed<double, 3> Vec;
  typedef std::map<Vec, vtkIdType> Hash;
  Hash hash;

  vtkPoints *pts = vtkPoints::New();
  out->SetPoints(pts);
  out->Allocate(p1->GetNumberOfCells() + p2->GetNumberOfCells());
  
  // Insert points into the output mesh from the first input
  for(vtkIdType i = 0; i < p1->GetNumberOfPoints(); i++)
    {
    hash[Vec(p1->GetPoint(i))] = out->GetPoints()->InsertNextPoint(p1->GetPoint(i));
    }
     
  // Insert points from the second input, unless they are duplicate
  for(vtkIdType i = 0; i < p2->GetNumberOfPoints(); i++)
    {
    Vec p(p2->GetPoint(i));
    typename Hash::iterator it = hash.find(p);
    if(it == hash.end())
      hash[p] = out->GetPoints()->InsertNextPoint(p.data_block());
    }

  // Insert cells from first array
  for(vtkIdType i = 0; i < p1->GetNumberOfCells(); i++)
    {
    vtkCell *cell = p1->GetCell(i);
    out->InsertNextCell(cell->GetCellType(), cell->GetNumberOfPoints(), cell->GetPointIds()->GetPointer(0));
    }

  // Insert cells from second array
  for(vtkIdType i = 0; i < p2->GetNumberOfCells(); i++)
    {
    vtkCell *cell = p2->GetCell(i);
    vtkIdType ids[1024];
    for(int j = 0; j < cell->GetNumberOfPoints(); j++)
      {
      Vec p(p2->GetPoint(cell->GetPointId(j)));
      ids[j] = hash[p];
      }
    out->InsertNextCell(cell->GetCellType(), cell->GetNumberOfPoints(), ids);
    }

  // Merge all the cell arrays
  for(int k = 0; k < p1->GetCellData()->GetNumberOfArrays(); k++)
    {
    int index;
    vtkDataArray *src1 = p1->GetCellData()->GetArray(k);
    vtkDataArray *src2 = p2->GetCellData()->GetArray(src1->GetName(), index);

    if(index < 0 || src1->GetNumberOfComponents() != src2->GetNumberOfComponents())
      continue;

    vtkFloatArray * arr = AddArrayToMesh(out, CELL, src1->GetName(), src1->GetNumberOfComponents(), 0);
  
    for(vtkIdType i = 0; i < p1->GetNumberOfCells(); i++)
      for(int m = 0; m < src1->GetNumberOfComponents(); m++)
        arr->SetComponent(i, m, src1->GetComponent(i, m));

    for(vtkIdType i = 0; i < p2->GetNumberOfCells(); i++)
      for(int m = 0; m < src1->GetNumberOfComponents(); m++)
        arr->SetComponent(i + p1->GetNumberOfCells(), m, src2->GetComponent(i, m));
    }

  // Merge all the point arrays
  for(int k = 0; k < p1->GetPointData()->GetNumberOfArrays(); k++)
    {
    int index;
    vtkDataArray *src1 = p1->GetPointData()->GetArray(k);
    vtkDataArray *src2 = p2->GetPointData()->GetArray(src1->GetName(), index);

    if(index < 0 || src1->GetNumberOfComponents() != src2->GetNumberOfComponents())
      continue;

    vtkFloatArray * arr = AddArrayToMesh(out, POINT, src1->GetName(), src1->GetNumberOfComponents(), 0);

    for(vtkIdType i = 0; i < p1->GetNumberOfPoints(); i++)
      for(int m = 0; m < src1->GetNumberOfComponents(); m++)
        arr->SetComponent(i, m, src1->GetComponent(i, m));

    for(vtkIdType i = 0; i < p2->GetNumberOfPoints(); i++)
      {
      Vec p(p2->GetPoint(i));
      vtkIdType ip = hash[p];
      for(int m = 0; m < src1->GetNumberOfComponents(); m++)
        arr->SetComponent(ip, m, src2->GetComponent(i, m));
      }
    }

  // Make sure the normals are preserved
  if(p1->GetPointData()->GetNormals() && p2->GetPointData()->GetNormals())
    out->GetPointData()->SetActiveNormals(p1->GetPointData()->GetNormals()->GetName());

  if(p1->GetCellData()->GetNormals() && p2->GetCellData()->GetNormals())
    out->GetPointData()->SetActiveNormals(p1->GetCellData()->GetNormals()->GetName());
}

class Edge : public std::pair<vtkIdType, vtkIdType>
{
public:
  Edge(vtkIdType a, vtkIdType b) : std::pair<vtkIdType,vtkIdType>(min(a,b), max(a,b)) {}
};

void PointDataDiffusion(vtkDataSet *mesh, double time, double dt, const char *array)
{
  // Diffusion simulates heat equation, dF/dt = -Laplacian(F), for t = time
  // We use the most basic approximation of the laplacian L(F) = [Sum_{j\in N(i)} F(j) - F(i)] / |N(i)|

  // Create a set of all edges in the mesh
  typedef set<Edge> EdgeSet;
  EdgeSet edges;

  // Report
  printf("Performing diffusion on point data (t = %f, delta_t = %f)\n", time, dt);

  // Source array
  vtkDataArray *f = mesh->GetPointData()->GetArray(array);

  // Get all edges into the edge set
  for(int i = 0; i < mesh->GetNumberOfCells(); i++)
    {
    vtkCell *cell = mesh->GetCell(i);
    vtkIdType *p = cell->GetPointIds()->GetPointer(0);
    if(cell->GetCellType() == VTK_TRIANGLE)
      {
      edges.insert(Edge(p[0], p[1]));
      edges.insert(Edge(p[1], p[2]));
      edges.insert(Edge(p[0], p[2]));
      }
    else if(cell->GetCellType() == VTK_TETRA)
      {
      edges.insert(Edge(p[0], p[1]));
      edges.insert(Edge(p[0], p[2]));
      edges.insert(Edge(p[0], p[3]));
      edges.insert(Edge(p[1], p[2]));
      edges.insert(Edge(p[1], p[3]));
      edges.insert(Edge(p[2], p[3]));
      }
    }

  // Count the number of neighbors of each vertex
  std::vector<int> nbr(mesh->GetNumberOfPoints(), 0);
  for(EdgeSet::iterator it = edges.begin(); it!=edges.end(); ++it)
    {
    nbr[it->first]++; 
    nbr[it->second]++;
    }

  // Create an array for the updates
  vtkFloatArray *f_upd = vtkFloatArray::New();
  f_upd->SetNumberOfComponents(f->GetNumberOfComponents());
  f_upd->SetNumberOfTuples(f->GetNumberOfTuples());

  // Copy f to f_upd
  for(int i = 0; i < f->GetNumberOfTuples(); i++)
    for(int j = 0; j < f->GetNumberOfComponents(); j++)
      f_upd->SetComponent(i, j, f->GetComponent(i, j));

  // Iterate
  unsigned int jt = 0;
  for(double t = 0; t < time - dt/2; t+=dt)
    {
    // Update f_upd
    for(EdgeSet::iterator it = edges.begin(); it!=edges.end(); ++it)
      {
      double wa = dt / nbr[it->first];
      double wb = dt / nbr[it->second];
      int ia = it->first;
      int ib = it->second;

      for(int j = 0; j < f->GetNumberOfComponents(); j++)
        {
        double val_a_j = f->GetComponent(ia,j);
        double val_b_j = f->GetComponent(ib,j);

        if(!vnl_math::isnan(val_b_j))
          f_upd->SetComponent(ia, j, f_upd->GetComponent(ia, j) + wa * (val_b_j - val_a_j));

        if(!vnl_math::isnan(val_a_j))
          f_upd->SetComponent(ib, j, f_upd->GetComponent(ib, j) + wb * (val_a_j - val_b_j));
        }
      }

    // Copy f_upd to f
    for(int i = 0; i < f->GetNumberOfTuples(); i++)
      for(int j = 0; j < f->GetNumberOfComponents(); j++)
        f->SetComponent(i, j, f_upd->GetComponent(i, j));

    cout << "." << flush;
    if((++jt) % 100 == 0 || t+dt >= time - 0.5 * dt)
      cout << " t = " << t+dt << endl;
    }
}

void CellDataDiffusion(vtkDataSet *mesh, double time, double dt, const char *array)
{
  // Diffusion, but between cells. This is really pretty ad hoc now
  
  // The links have to be built, but this function is only available through
  // classes vtkPolyData, vtkUnstructuredGrid. So we have to cast
  vtkPolyData *mesh_pd = dynamic_cast<vtkPolyData *>(mesh);
  vtkUnstructuredGrid *mesh_ug = dynamic_cast<vtkUnstructuredGrid *>(mesh);
  if (mesh_pd)
    mesh_pd->BuildLinks();
  else if (mesh_ug)
    mesh_ug->BuildLinks();
  else 
    throw MCException("Unexpected mesh type in CellDataDiffusion");

  // Create a set of all edges in the mesh. These are pairs of adjacent cells that
  // share an edge
  typedef set<Edge> EdgeSet;
  EdgeSet edges;

  // Report
  printf("Performing diffusion on cell data (t = %f, delta_t = %f)\n", time, dt);

  // Get all edges into the edge set
  for(int i = 0; i < mesh->GetNumberOfCells(); i++)
    {
    vtkCell *cell = mesh->GetCell(i);
    if(cell->GetCellType() == VTK_TETRA)
      {
      for(int j = 0; j < cell->GetNumberOfFaces(); j++)
        {
        vtkSmartPointer<vtkIdList> nbr = vtkIdList::New();
        vtkCell *face = cell->GetFace(j);
        mesh->GetCellNeighbors(i, face->GetPointIds(), nbr);
        for(int k = 0; k < nbr->GetNumberOfIds(); k++)
          edges.insert(Edge(i, nbr->GetId(k)));
        }
      }
    else if(cell->GetCellType() == VTK_TRIANGLE)
      {
      for(int j = 0; j < cell->GetNumberOfEdges(); j++)
        {
        vtkSmartPointer<vtkIdList> nbr = vtkIdList::New();
        vtkCell *edge = cell->GetEdge(j);
        mesh->GetCellNeighbors(i, edge->GetPointIds(), nbr);
        for(int k = 0; k < nbr->GetNumberOfIds(); k++)
          edges.insert(Edge(i, nbr->GetId(k)));
        }
      }
    else throw MCException("Wrong cell type in CellDataDiffusion");
    }

  printf("There are %d pairs of adjacent cells\n", (int) edges.size());

  // Count the number of neighbors of each cell
  std::vector<int> nbr(mesh->GetNumberOfCells(), 0);
  for(EdgeSet::iterator it = edges.begin(); it!=edges.end(); ++it)
    {
    nbr[it->first]++; 
    nbr[it->second]++;
    }

  // Create an array for the updates
  vtkDataArray *f = mesh->GetCellData()->GetArray(array);
  vtkFloatArray *f_upd = vtkFloatArray::New();
  f_upd->SetNumberOfComponents(f->GetNumberOfComponents());
  f_upd->SetNumberOfTuples(f->GetNumberOfTuples());

  // Copy f to f_upd
  for(int i = 0; i < f->GetNumberOfTuples(); i++)
    for(int j = 0; j < f->GetNumberOfComponents(); j++)
      f_upd->SetComponent(i, j, f->GetComponent(i, j));

  // Iterate
  unsigned int jt = 0;
  for(double t = 0; t < time - dt/2; t+=dt)
    {
    // Update f_upd
    for(EdgeSet::iterator it = edges.begin(); it!=edges.end(); ++it)
      {
      double wa = dt / nbr[it->first];
      double wb = dt / nbr[it->second];
      int ia = it->first;
      int ib = it->second;

      for(int j = 0; j < f->GetNumberOfComponents(); j++)
        {
        double val_a_j = f->GetComponent(ia,j);
        double val_b_j = f->GetComponent(ib,j);

        if(!vnl_math::isnan(val_b_j))
          f_upd->SetComponent(ia, j, f_upd->GetComponent(ia, j) + wa * (val_b_j - val_a_j));

        if(!vnl_math::isnan(val_a_j))
          f_upd->SetComponent(ib, j, f_upd->GetComponent(ib, j) + wb * (val_a_j - val_b_j));

        /*
        f_upd->SetComponent(ia, j, 
          f_upd->GetComponent(ia, j) + wa * (f->GetComponent(ib,j) - f->GetComponent(ia, j)));
        f_upd->SetComponent(ib, j, 
          f_upd->GetComponent(ib, j) + wb * (f->GetComponent(ia,j) - f->GetComponent(ib, j)));
          */
        }
      }

    // Copy f_upd to f
    for(int i = 0; i < f->GetNumberOfTuples(); i++)
      for(int j = 0; j < f->GetNumberOfComponents(); j++)
        f->SetComponent(i, j, f_upd->GetComponent(i, j));

    cout << "." << flush;
    if((++jt) % 100 == 0 || t+dt >= time - 0.5 * dt)
      cout << " t = " << t+dt << endl;
    }
}

template <class TMeshType>
ClusterArray ComputeClusters(
  TMeshType *mesh,
  Domain dom,
  double thresh, 
  TMeshType **mout = NULL)
{ ClusterArray ca(1) ; return ca; }

class ClusterComputer
{
public:
  ClusterComputer(vtkDataSet *mesh, Domain dom, double thresh);
  ClusterArray ComputeClusters(bool build_combo_mesh);
  vtkUnstructuredGrid *GetFullMesh()
    { return mout; }

private:
  vtkSmartPointer<vtkConnectivityFilter> fConnect;
  vtkSmartPointer<vtkClipDataSet> fClip;
  vtkSmartPointer<vtkThreshold> fThresh, fThreshInv;
  vtkSmartPointer<vtkGeometryFilter> fGeometry;

  vtkSmartPointer<vtkDataSet> mesh;
  vtkSmartPointer<vtkUnstructuredGrid> mout;

  double thresh;
  Domain dom;
  char *scalar_name;
};

ClusterComputer
::ClusterComputer(vtkDataSet *mesh, Domain dom, double thresh)
{
  this->thresh = thresh;
  this->dom = dom;
  this->mesh = mesh;

  // The connectivity filter
  fConnect = vtkConnectivityFilter::New();
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();

  // Generate pipeline
  if (dom == POINT)
    {
    fClip = vtkClipDataSet::New();
    fClip->SetInputData(mesh);
    fClip->SetInputArrayToProcess(0, 0, 0, 
      vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS); 
    fClip->SetValue(thresh);
    fConnect->SetInputConnection(fClip->GetOutputPort());
    } 
  else
    {  
    fThresh = vtkThreshold::New();
    fThresh->SetInputData(mesh);
    fThresh->SetInputArrayToProcess(0, 0, 0, 
      vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
    fThresh->ThresholdByUpper(thresh);
    fConnect->SetInputConnection(fThresh->GetOutputPort());

    fThreshInv = vtkThreshold::New();
    fThreshInv->SetInputData(mesh);
    fThreshInv->SetInputArrayToProcess(0, 0, 0, 
      vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
    fThreshInv->ThresholdByLower(thresh - 1.e-6);
    } 

  // Initialize the output object
  mout = vtkUnstructuredGrid::New();

  // Copy the names of the scalars
  scalar_name = GetScalarsFromMesh(mesh, dom)->GetName();
}

double TetraVolume(vtkDataSet *mesh, 
  vtkIdType a0, vtkIdType a1, vtkIdType a2, vtkIdType a3)
{
  return vtkTetra::ComputeVolume(
    mesh->GetPoint(a0),
    mesh->GetPoint(a1),
    mesh->GetPoint(a2),
    mesh->GetPoint(a3));
}

double ComputeCellVolume(vtkCell *cell)
{
  vtkSmartPointer<vtkIdList> ids;
  if(cell->GetCellType() == VTK_TETRA)
    {
    ids = cell->GetPointIds();
    }
  else if(cell->GetCellType() == VTK_WEDGE)
    {
    ids = vtkIdList::New();
    vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
    cell->Triangulate(0, ids, pts);
    }
  else throw MCException("Volume calculation only supported for Tetra and Wedge cells");

  double vol = 0;
  for(int k = 0; k < ids->GetNumberOfIds(); k+=4)
    {
    vol += fabs(vtkTetra::ComputeVolume(
        cell->GetPoints()->GetPoint(ids->GetId(k)),
        cell->GetPoints()->GetPoint(ids->GetId(k+1)),
        cell->GetPoints()->GetPoint(ids->GetId(k+2)),
        cell->GetPoints()->GetPoint(ids->GetId(k+3))));
    }
  return vol;
}

ClusterArray
ClusterComputer::ComputeClusters(bool build_combo_mesh)
{
  // Do we need clipped output?
  if(dom == POINT)
    fClip->SetGenerateClippedOutput(build_combo_mesh);
 
  // Perform clipping and connectivity analysis
  if(dom == POINT)
    fClip->Modified();
  else
    fThresh->Modified();

  fConnect->Update();
  vtkUnstructuredGrid *p = dynamic_cast<vtkUnstructuredGrid *>(fConnect->GetOutput());
  int nelt = (dom == POINT) ?  p->GetNumberOfPoints() : p->GetNumberOfCells();

  // Leave if there are no clusters
  if(nelt == 0)
    return ClusterArray(0);

  // Create a temporary array to store cell area element
  vector<double> daArea(nelt, 0.0);

  // Compute the area of each triangle in the cluster set
  for(int k = 0; k < p->GetNumberOfCells(); k++)
    {
    vtkCell *cell = p->GetCell(k);
    if(cell->GetCellType() == VTK_TRIANGLE)
      {
      // Compute the area of the triangle
      double p0[3], p1[3], p2[3];
      vtkIdType a0 = cell->GetPointId(0), a1 = cell->GetPointId(1), a2 = cell->GetPointId(2);
      p->GetPoint(a0, p0); p->GetPoint(a1, p1); p->GetPoint(a2, p2); 
      double area = vtkTriangle::TriangleArea(p0, p1, p2);

      if(dom == POINT)
        {  
        // Split the volume between neighbors
        daArea[a0] += area / 3.0; daArea[a1] += area / 3.0; daArea[a2] += area / 3.0;
        }
      else
        {
        // No need to split, working with cells
        daArea[k] = area;
        }
      }
    else if(cell->GetCellType() == VTK_TETRA || cell->GetCellType() == VTK_WEDGE)
      {
      double vol = ComputeCellVolume(cell);
      if(dom == POINT)
        {
        for(int m = 0; m < cell->GetNumberOfPoints(); m++)
          daArea[cell->GetPointId(m)] += vol / cell->GetNumberOfPoints();
        }
      else
        {
        daArea[k] = vol;
        }
      }
    else
      throw MCException("Wrong cell type, should be VTK_TRIANGLE, VTK_TETRA or VTK_WEDGE");

    // Compute the area of the triangle
    double p0[3], p1[3], p2[3];
    vtkIdType a0 = cell->GetPointId(0), a1 = cell->GetPointId(1), a2 = cell->GetPointId(2);
    p->GetPoint(a0, p0); p->GetPoint(a1, p1); p->GetPoint(a2, p2); 
    double area = vtkTriangle::TriangleArea(p0, p1, p2);

    if(dom == POINT)
      {  
      // Split the volume between neighbors
      daArea[a0] += area / 3.0; daArea[a1] += area / 3.0; daArea[a2] += area / 3.0;
      }
    else
      {
      // No need to split, working with cells
      daArea[k] = area;
      }
    }

  // Get the input statistic
  vtkDataArray *daStat = GetArrayFromMesh(p, dom, scalar_name);

  // Get the region ID
  vtkDataArray *daRegion = GetScalarsFromMesh(p, dom);

  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
  for(int i = 0; i < nelt; i++)
    {
    size_t region = (size_t) (daRegion->GetTuple1(i));
    double x = daStat->GetTuple1(i);
    double a = daArea[i];
    ca[region].n++;
    ca[region].area += a;
    ca[region].power += a * x;
    ca[region].tvalue += x;
    }

  // Get the output if needed
  if(build_combo_mesh)
    {
    if(dom == POINT)
      {
      // Get the stuff that got cut out
      fClip->InsideOutOn();
      fClip->Update();
      vtkUnstructuredGrid *fout = fClip->GetOutput();
      fClip->InsideOutOff();

      // Assign each cell in cut away piece a ClusterId of -1
      AddArrayToMesh(fout, CELL, "RegionId", 1, -1.0);

      // Assign a RegionId to each cell in the thresholded part
      vtkFloatArray *prid = AddArrayToMesh(p, CELL, "RegionId", 1, 0.0);

      for(int i = 0; i < p->GetNumberOfCells(); i++)
        {
        vtkCell *cell = p->GetCell(i);
        double c1 = daRegion->GetTuple1(cell->GetPointId(0));
        double c2 = daRegion->GetTuple1(cell->GetPointId(1));
        double c3 = daRegion->GetTuple1(cell->GetPointId(2));
        if(c1 != c2 || c1 != c3 || c2 != c3)
          throw MCException("Inconsistent RegionID mapping");

        prid->SetTuple1(i, c1);
        }

      // Merge the two arrays
      MergeMeshPieces<vtkUnstructuredGrid>(p, fout, mout);
      mout->GetPointData()->SetActiveScalars(scalar_name);
      cout << "SC = " <<  mout->GetPointData()->GetScalars() << endl;
      }
    else
      {
      // Apply an inverted threshold
      fThreshInv->Update();
      vtkUnstructuredGrid *cut = fThreshInv->GetOutput();

      // Each cell in p already has a region id. Just need to add regionids to the cut
      AddArrayToMesh(cut, CELL, "RegionId", 1, -1);

      // Merge the pieces
      MergeMeshPieces<vtkUnstructuredGrid>(p, cut, mout);
      mout->GetCellData()->SetActiveScalars(scalar_name);
      }
    }

  // Return the cluster array
  return ca;
}

template <>
ClusterArray ComputeClusters(
  vtkPolyData *mesh,
  Domain dom,
  double thresh, 
  vtkPolyData **mout)
{
  // This is the thresholded out / clipped out region
  vtkSmartPointer<vtkDataSet> fin = NULL, fout = NULL; 

  if (dom == POINT)
	{
    // Initialize mesh
    vtkSmartPointer<vtkClipPolyData> fContour = vtkClipPolyData::New();

    // Clip the data field at the threshold
    fContour->SetInputData(mesh);
    fContour->SetValue(thresh);

    // Generate the background too 
    if(mout)
      fContour->GenerateClippedOutputOn();

    // Run the filter
    fContour->Update();
    fin = fContour->GetOutput();

    if(mout)
      fout = fContour->GetClippedOutput();
    } 
  else
    {  
    // Initialize mesh
    vtkSmartPointer<vtkThreshold> fThresh = vtkThreshold::New();

    // Threshold the cell data field
    fThresh->SetInputData(mesh);
    fThresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
    fThresh->ThresholdByUpper(thresh);
    fThresh->Update();
    fin = fThresh->GetOutput();
    } 

  // Get the connected components
  vtkSmartPointer<vtkConnectivityFilter> fConnect = vtkConnectivityFilter::New();
  fConnect->SetInputData(fin);
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  fConnect->Update();

  // Convert the result to polydata (why?)
  vtkSmartPointer<vtkGeometryFilter> geo = vtkGeometryFilter::New();
  geo->SetInputData(fConnect->GetOutput());
  geo->Update();
  
  vtkSmartPointer<vtkPolyData> p = geo->GetOutput();

  // Number of elements to process
  int nelt = (dom == POINT) ? 
    p->GetNumberOfPoints() : p->GetNumberOfCells();

  // Create output data arrays for computing area element
  vtkFloatArray * daArea = AddArrayToMesh(p, dom, "area_element", 1, 0.0);

  // Compute the area of each triangle in the cluster set
  for(int k = 0; k < p->GetNumberOfCells(); k++)
    {
    vtkCell *cell = p->GetCell(k);
    if(cell->GetCellType() != VTK_TRIANGLE)
      throw MCException("Wrong cell type, should be VTK_TRIANGLE");

    // Compute the area of the triangle
    double p0[3], p1[3], p2[3];
    vtkIdType a0 = cell->GetPointId(0), a1 = cell->GetPointId(1), a2 = cell->GetPointId(2);
    p->GetPoint(a0, p0); p->GetPoint(a1, p1); p->GetPoint(a2, p2); 
    double area = vtkTriangle::TriangleArea(p0, p1, p2);

    if(dom == POINT)
      {  
      // Split the volume between neighbors
      daArea->SetTuple1(a0, area / 3.0 + daArea->GetTuple1(a0));
      daArea->SetTuple1(a1, area / 3.0 + daArea->GetTuple1(a1));
      daArea->SetTuple1(a2, area / 3.0 + daArea->GetTuple1(a2));
      }
    else
      {
      // No need to split, working with cells
      daArea->SetTuple1(k, area);
      }
    }

  // The scalars in fin are the input statistic
  vtkDataArray * daStat = GetScalarsFromMesh(fin, dom);

  // The scalars in p are the region IDs
  vtkDataArray * daRegion = GetScalarsFromMesh(p, dom);

  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
  for(int i = 0; i < nelt; i++)
    {
    size_t region = (size_t) (daRegion->GetTuple1(i));
    double x = daStat->GetTuple1(i);
    double a = daArea->GetTuple1(i);
    ca[region].n++;
    ca[region].area += a;
    ca[region].power += a * x;
    ca[region].tvalue += x;
    }

  // Get the output if needed
  if(mout != NULL)
    {
    if(dom == POINT)
      {
      // Assign each cell in cut away piece a ClusterId of -1
      AddArrayToMesh(fout, CELL, "RegionId", 1, -1.0);

      // Assign a RegionId to each cell in the thresholded part
      vtkFloatArray * prid = AddArrayToMesh(p, CELL, "RegionId", 1, 0.0);

      for(int i = 0; i < p->GetNumberOfCells(); i++)
        {
        vtkCell *cell = p->GetCell(i);
        double c1 = daRegion->GetTuple1(cell->GetPointId(0));
        double c2 = daRegion->GetTuple1(cell->GetPointId(1));
        double c3 = daRegion->GetTuple1(cell->GetPointId(2));
        if(c1 != c2 || c1 != c3 || c2 != c3)
          throw MCException("Inconsistent RegionID mapping");

        prid->SetTuple1(i, c1);
        }

      // Merge the two arrays
      *mout = vtkPolyData::New();
      MergeMeshPieces(p, fout, *mout);
      (*mout)->GetPointData()->SetActiveScalars(daStat->GetName());
      }
    else
      {
      // Apply an inverted threshold
      vtkThreshold *fThreshInv = vtkThreshold::New();

      // Threshold the cell data field
      fThreshInv->SetInputData(mesh);
      fThreshInv->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
      fThreshInv->ThresholdByLower(thresh - 1.e-6);
      fThreshInv->Update();
      vtkUnstructuredGrid *cut = fThreshInv->GetOutput();
      fThreshInv->Delete();

      // Each cell in p already has a region id. Just need to add regionids to the cut
      AddArrayToMesh(cut, CELL, "RegionId", 1, -1);

      // Merge the pieces
      *mout = vtkPolyData::New();
      MergeMeshPieces(p, cut, *mout);
      (*mout)->GetCellData()->SetActiveScalars(daStat->GetName());
      }
    }

  // Return the cluster array
  return ca;
}

template <>
ClusterArray ComputeClusters(
  vtkUnstructuredGrid *mesh,
  Domain dom,
  double thresh, 
  vtkUnstructuredGrid **mout )
{
  vtkClipDataSet *fContour = NULL;
  vtkThreshold *fThresh = NULL;

  vtkUnstructuredGrid *f;
  if(dom == POINT)
    {  
    // Initialize mesh
    fContour = vtkClipDataSet::New();

    // Clip the data field at the threshold
    fContour->SetInputData(mesh);
    fContour->SetValue(thresh);
    fContour->Update();
    f = fContour->GetOutput();
    }
  else
    {  
     // Initialize mesh
     fThresh = vtkThreshold::New();

     // Threshold the cell data field
     fThresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
     fThresh->SetInputData(mesh);
     fThresh->ThresholdByUpper(thresh);
     fThresh->Update();
     f = fThresh->GetOutput();
  }

  vtkDataSetAttributes *fdata = (dom == POINT) ? 
    (vtkDataSetAttributes *) f->GetPointData() : 
    (vtkDataSetAttributes *) f->GetCellData();

  // Get the connected components
  vtkConnectivityFilter * fConnect = vtkConnectivityFilter::New();
  fConnect->SetInputData(f);
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  fConnect->Update();
  vtkUnstructuredGrid *p = dynamic_cast<vtkUnstructuredGrid *>(fConnect->GetOutput());

  vtkDataSetAttributes *pdata = (dom == POINT) ? 
    (vtkDataSetAttributes *) p->GetPointData() : 
    (vtkDataSetAttributes *) p->GetCellData();

  // Number of elements to process
  int nelt = (dom == POINT) ? 
    p->GetNumberOfPoints() : p->GetNumberOfCells();

  // Create output data arrays for computing volume element
  vtkFloatArray *daArea = vtkFloatArray::New();
  daArea->SetName("volume_element");
  daArea->SetNumberOfComponents(1);
  daArea->SetNumberOfTuples(nelt);
  daArea->FillComponent(0, 0.0);

  // Compute the volume of each tetra in the cluster set
  for(int k = 0; k < p->GetNumberOfCells(); k++)
    {
    vtkCell *cell = p->GetCell(k);

    // Do this, otherwise output is meaningless
    if(cell->GetCellType() != VTK_TETRA)
      throw MCException("Wrong cell type");

    // Compute the volume of the tetra
    vtkIdType a0 = cell->GetPointId(0);
    vtkIdType a1 = cell->GetPointId(1);
    vtkIdType a2 = cell->GetPointId(2);
    vtkIdType a3 = cell->GetPointId(3);
    double p0[3], p1[3], p2[3], p3[3];
    p->GetPoint(a0, p0);
    p->GetPoint(a1, p1);
    p->GetPoint(a2, p2);
    p->GetPoint(a3, p3);

    double area = vtkTetra::ComputeVolume(p0, p1, p2, p3);

    if(dom == POINT)
      {  
      // Split the volume between neighbors
      daArea->SetTuple1(a0, abs(area) / 4.0 + daArea->GetTuple1(a0));
      daArea->SetTuple1(a1, abs(area) / 4.0 + daArea->GetTuple1(a1));
      daArea->SetTuple1(a2, abs(area) / 4.0 + daArea->GetTuple1(a2));
      daArea->SetTuple1(a3, abs(area) / 4.0 + daArea->GetTuple1(a3));
      }
    else
      {
      // No need to split, working with cells
      daArea->SetTuple1(k, abs(area));
      }
    }

  vtkDataArray *daRegion = pdata->GetScalars();
  vtkDataArray *daData = fdata->GetScalars();

  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
  for(int i = 0; i < nelt; i++)
    {
    size_t region = (size_t) (daRegion->GetTuple1(i));
    double x = daData->GetTuple1(i);
    double a = daArea->GetTuple1(i);
    ca[region].n++;
    ca[region].area += a;
    ca[region].power += a * x;
    ca[region].tvalue += x;
    }

  // Get the output if needed
  if(mout != NULL)
    {
    pdata->AddArray(daArea);
    *mout = p;
    }
  else
    {
    // Delete the intermediates
    daArea->Delete();
    fConnect->Delete();
    if (dom==POINT)
       fContour->Delete();
    else
       fThresh->Delete();
    }

  // Return the cluster array
  return ca;
}

class GeneralLinearModel
{
public:
  GeneralLinearModel(
    vtkDataArray *data, vtkDataArray *contrast, 
    vtkDataArray *tstat, vtkDataArray *pval, vtkDataArray *beta, vtkDataArray *residual, vtkDataArray *dfarray,
    const vnl_matrix<double> &design_matrix, 
    const vnl_matrix<double> &contrast_vector);

  void Compute(const vnl_matrix<double> &Yperm, bool need_t, bool need_p, bool need_res);
  void ComputeWithMissingData(const vnl_matrix<double> &Yperm, bool need_t, bool need_p, bool need_res);
  //int compute_rank(vnl_matrix <double> D); 
  double PtoT(double p);
  vnl_matrix<double> Y; // data matrix for permutations

private:
  // Pointers to the data arrays
  vtkDataArray *data, *contrast, *tstat, *pval, *beta, *residual, *dfarray;

  // Copies of the matrix, other data
  vnl_matrix<double> X,cv;
  int nelt, rank, df;
};

GeneralLinearModel::GeneralLinearModel(
    vtkDataArray *data, vtkDataArray *contrast, 
    vtkDataArray *tstat, vtkDataArray *pval, vtkDataArray *beta,
    vtkDataArray *residual, vtkDataArray *dfarray,
    const vnl_matrix<double> &design_matrix, 
    const vnl_matrix<double> &contrast_vector)
{
  // Copy input
  this->data = data;
  this->contrast = contrast;
  this->tstat = tstat;
  this->pval = pval;
  this->beta = beta;
  this->X = design_matrix;
  this->cv = contrast_vector;
  this->residual = residual;
  this->dfarray = dfarray;

  // Initialize other stuff
  nelt = data->GetNumberOfTuples();

  // Rank and degrees of freedom
  rank = vnl_rank(X.transpose() * X, vnl_rank_row);
  df = X.rows() - rank;
 
  // Y matrix
  Y.set_size(X.rows(), nelt);
  for(size_t i = 0; i < X.rows(); i++)
    for(int j = 0; j < nelt; j++)
      Y(i,j) = data->GetComponent(j,i); 
  
}

void
GeneralLinearModel::Compute(const vnl_matrix<double> &Yperm, bool need_t, bool need_p, bool need_res)
{

  // Shuffle the rows according to current permutan_resr (size_t i = 0; i < permutation.size(); i++)
  //  for (size_t j = 0; j < X.cols(); j++)
  //    	Xperm[i][j] = X[permutation[i]][j];

  // Standard GLM calculation
  vnl_matrix<double> A = vnl_matrix_inverse<double>(X.transpose() * X).pinverse(rank);
  vnl_matrix<double> bhat = (A * X.transpose()) * Y;
  
  // Compute the contrast
  vnl_matrix<double> res = cv * bhat;

  // Copy the betas and the contrast to the output vector
  for (int j = 0; j < nelt; j++)
    {
    contrast->SetTuple1(j, res(0,j));
    beta->SetTuple(j, bhat.get_column(j).data_block());
    }

  // The rest only if we need the t-stat
  if(need_t)
    {
    vnl_matrix<double> errmat = Y - X * bhat;
    
    // Residual variance
    vnl_matrix<double> resvar(1, nelt);
    resvar.fill(0.0);
    for(int j = 0; j < nelt; j++)
       {
       for(size_t i = 0; i < X.rows(); i++)
          resvar(0, j) += errmat(i,j) * errmat(i,j);
       if(j==0){
	 std::cout << "res_mag" << resvar(0,j) << std::endl;
	std::cout << "df" << df << std::endl;
	}
       resvar(0, j) = resvar(0, j) / (double) df;
       }

    //if we need the residual, copy to the output vector
    if(need_res)
    {	
      for (int j = 0; j< nelt; j++)
      {
      residual->SetTuple(j, errmat.get_column(j).data_block());
      }
    }

    // Compute t-stat / p-value
    vnl_matrix<double> den(1, nelt);
    den = (cv * (A * cv.transpose())) * resvar;
    for (int j = 0; j < nelt; j++)
      {
      double t = den(0, j) > 0 ? res(0, j) / sqrt(den(0, j)) : 0.0;
      tstat->SetTuple1(j, t);
      dfarray->SetTuple1(j, (double) df);
      if(need_p)
        {
        int dummy;
        pval->SetTuple1(j, 1.0 - tnc(t, df, 0.0, &dummy)); 
        }
      }
    }

}

void
GeneralLinearModel::ComputeWithMissingData(const vnl_matrix<double> &Yperm, bool need_t, bool need_p, bool need_res)
{
  // Shuffle the rows according to current permutation
  //for (size_t i = 0; i < permutation.size(); i++)
  //  for (size_t j = 0; j < X.cols(); j++)
  //    Xperm[i][j] = X[permutation[i]][j];

  // Number of subjects
  int ns = Yperm.rows();

  // Matrix A
  vnl_matrix<double> A, Xj, AXjT;
  vnl_vector<double> Yj;
  double cv_A_cvT;

  // Iterate over all the elements
  for (int j = 0; j < nelt; j++)
    {
    
    bool same_pattern = true;
    int n_nans = 0;

    // Determine which elements have missing data (nan) and if that matches the last element we processed
    for(int k = 0; k < ns; k++)
      {
      if(vnl_math::isnan(Yperm(k, j)))
        {
        n_nans++;
        if(j == 0 || !vnl_math::isnan(Yperm(k,j-1)))
          same_pattern = false;
        }
      else
        {
        if(j == 0 || vnl_math::isnan(Yperm(k,j-1)))
          same_pattern = false;
        }
      }

    // If the pattern is changed, we need to compute the new A matrix
    if(!same_pattern)
      {
      // Copy rows from Xperm into Xj
      Xj.set_size(ns - n_nans, X.cols());
      for(int k = 0, p = 0; k < ns; k++)
        if(!vnl_math::isnan(Yperm(k,j)))
          Xj.set_row(p++, X.get_row(k));

      
      // Compute A - should be Xj?
      A = vnl_matrix_inverse<double>(Xj.transpose() * Xj).pinverse(rank);
      AXjT = A * Xj.transpose();

      // Compute cv * A * cv^T
      cv_A_cvT = (cv * (A * cv.transpose()))(0,0);

      // Resize Y
      Yj.set_size(ns - n_nans);
      }

    // Copy the Y vector's elements into Yk
    for(int k = 0, p = 0; k < ns; k++)
      if(!vnl_math::isnan(Yperm(k,j)))
        Yj[p++] = Yperm(k,j);
    
    // Compute the estimated betas
    vnl_vector<double> beta_j = AXjT * Yj;

    // Store in output array
    beta->SetTuple(j, beta_j.data_block());

    // Compute the contrast - con_j same as res
    vnl_vector<double> con_j = cv * beta_j;
    contrast->SetTuple1(j, con_j[0]);

    //Compute the degrees of freedom
    int df_j = Xj.rows()-rank;

    // The rest only if we need the t-stat
    if(need_t)
      {
      vnl_vector<double> Xj_beta = Xj * beta_j;
      vnl_vector<double> res_j = Yj - Xj * beta_j;

      double resvar = res_j.squared_magnitude() / (double) df_j;

      //if needed, copy the residuals to the output array
      if(need_res)
        {
        vnl_vector<double> res_j_ns;
        res_j_ns.set_size(ns);
        for(int k = 0, p = 0; k < ns; k++)
          if(!vnl_math::isnan(Yperm(k,j)))
            res_j_ns[k] = res_j[p++];
          else
            res_j_ns[k] = NAN;
        residual->SetTuple(j,res_j_ns.data_block());
        }

      // Compute t-stat / p-value
      double den = cv_A_cvT * resvar;
      double t = (den > 0) ? con_j[0] / sqrt(den) : 0.0;
      tstat->SetTuple1(j, t);
      dfarray->SetTuple1(j, (double) df_j);
      if(need_p)
        {
        int dummy;
        pval->SetTuple1(j, 1.0 - tnc(t, df_j, 0.0, &dummy));
        }
      }
    }
}


//int
//GeneralLinearModel::compute_rank(vnl_matrix<double> D)
//{

//  int n = D.rows();
//  int m = D.cols();
//  const double eps = 1e-9;
//  int rank = 0;
//  vector<bool> row_selected(n,false);
//  for ( int i = 0; i < m; i++){
//	int j;
//	for ( j = 0; j < n; j++){
//		if( !row_selected[j] && abs(D[j][i] > eps))
//			break;
//	}
//	if( j != n){
//		++rank;
//		row_selected[j] = true;
//		for ( int p = i+1; p<m; ++p)
//			D[j][p] /= D[j][i];
//		for ( int k = 0; k < n; ++k) {
//			if ( k != j && abs(D[k][i]) > eps){
//				for(int p = i+1; p<m; ++p)
//					D[k][p] -= D[j][p]*D[k][i];
//			}
//		}
//	}
//}
//  return rank;

//}

double
GeneralLinearModel::PtoT(double p)
{
  // We do a binary search for the t
  int dummy;
  double t0 = 0, t1 = 1e10; 
  double p0 = 1.0 - tnc(t0, df, 0.0, &dummy);
  double p1 = 1.0 - tnc(t1, df, 0.0, &dummy);
  double t = 0, pt = 0;

  if(p > p0 || p < p1)
    throw MCException("Can not find t-stat for p threshold %f", p);

  for(int k = 0; k < 64; k++)
    {
    t = 0.5 * (t0 + t1);
    pt = 1.0 - tnc(t, df, 0.0, &dummy);
    if(p > pt)
      { t1 = t; p1 = pt; }
    else
      { t0 = t; p0 = pt; }
    if(fabs(p - pt) < 1e-7)
      break;
    }
  return t;
}

/*
template <class TMeshType>
void GeneralLinearModel(
  const vnl_matrix<double> &origmat, 
  const vnl_matrix<double> &origcon, 
  const vector<int> &indiv_labels, 
  TMeshType *pd, 
  const char *var, 
  std::string VOIttest, 
  Domain dom)
{
  // What data are we working on ?
  vtkDataSetAttributes *dsa = (dom == POINT) ? 
    (vtkDataSetAttributes *) pd->GetPointData() : 
    (vtkDataSetAttributes *) pd->GetCellData();
  
  // Get relevant arrays
  vtkDataArray *data = dsa->GetArray(var);
  vtkDataArray *ttest = dsa->GetArray(VOIttest.c_str());

  int nelt = data->GetNumberOfTuples();

  // Set up design matrix
  vnl_matrix<double> mat;
  mat.set_size(origmat.rows(), origmat.cols());

  // Shuffle the rows according to current permutation
  for (size_t i = 0; i < indiv_labels.size(); i++)
    for (size_t j = 0; j < origmat.cols(); j++)
      mat[i][j] = origmat[indiv_labels[i]][j];

  // copy contrast matrix
  vnl_matrix<double> con = origcon;

  // Load all images into a Y matrix (can we do this)
  vnl_matrix<double> Y(mat.rows(), nelt);
  for(size_t i = 0; i < mat.rows(); i++)
    for(size_t j = 0; j < nelt; j++)
      Y(i,j) = data->GetComponent(j,i); 

  // Compute degrees of freedom
  unsigned int rank = vnl_rank(mat, vnl_rank_row);
  unsigned int df = Y.rows() - rank;

  // Standard GLM calculation
  vnl_matrix<double> A = vnl_matrix_inverse<double>(mat.transpose() * mat).pinverse(rank);
  vnl_matrix<double> bhat = (A * mat.transpose()) * Y;
  
  // Compute the contrast
  vnl_matrix<double> res = con * bhat;
  vnl_matrix<double> errmat = Y - mat * bhat;

  // Residual variance
  vnl_matrix<double> resvar(1, nelt);
  resvar.fill(0.0);
  for(size_t j = 0; j < nelt; j++)
     {
     for(size_t i = 0; i < mat.rows(); i++)
        resvar(0, j) += errmat(i,j)*errmat(i,j);
     resvar(0,j) = resvar(0,j)/(double)df;
     }

  // t-value
  vnl_matrix<double> tmat(1, nelt);
  vnl_matrix<double> den(1, nelt);

  den = (con * (A * con.transpose())) * resvar;
  tmat.fill(0.0);
  for (size_t j = 0; j < nelt; j++)
    if ( den(0,j) > 0 )
      { 
      tmat(0,j) = res(0,j) / sqrt(den(0,j)); 
      } 

  // Write the output 
  for(size_t j = 0; j < nelt; j++)
    ttest->SetTuple1(j, tmat(0,j)); 
}
*/

enum ThresholdType { CONTRAST, TSTAT, PVALUE };

struct Parameters
{
  // Input and output meshes
  vector<string> fn_mesh_input, fn_mesh_output;

  // Number of permutations
  size_t np;
  
  // Array
  string array_name, exclusion_array_name;

  // Domain
  Domain dom;

  // GLM parameters
  string glm_design, glm_contrast;

  // Freedman-Lane parameters
  string fl_nuissance;

  // Row mask
  string row_mask;

  // T statistic threshold
  double threshold;

  // Diffusion amount
  double diffusion, delta_t;

  // Whether tubes are being used
  bool flag_edges;

  // Whether the data has missing observations encoded as nans
  bool flag_missing_data;

  // Whether to write output in binary format
  bool flag_write_binary;

  // Whether to use the freedman-lane procedure
  bool flag_freedman_lane;

  // Triangulate?
  bool flag_triangulate;

  // What is the threshold on
  ThresholdType ttype;

  Parameters()
    {
    np = 0; dom = POINT; threshold = 0; diffusion = 0; flag_edges = false; ttype = TSTAT; delta_t = 0.01; 
    flag_missing_data = false; flag_write_binary = false; flag_freedman_lane = false; flag_triangulate = false;
    }
};

class IdHash : public std::vector<vtkIdType>
{
public:
  IdHash(unsigned int dim, vtkIdType *src) : std::vector<vtkIdType>(dim, 0) 
    {
    for(unsigned int i = 0; i < dim; i++)
      (*this)[i] = src[i];
    std::sort(this->begin(), this->end());
    }

  IdHash() : std::vector<vtkIdType>() {}

  bool operator < (const IdHash &b)
    {
    if(this->size() != b.size())
      throw MCException("2D and 3D cells mixed");

    for(unsigned int k = 0; k < size(); k++)
      {
      if((*this)[k] < b[k]) return true;
      if((*this)[k] > b[k]) return false;
      }
    return false;
    }
};

template <class TMeshType>
int AppendOpenEdgesFromMesh(
  TMeshType *src, vtkDataArray *regarr, int region, TMeshType *trg)
{
  // Hashtable
  typedef IdHash HashEntry;
  typedef std::map<HashEntry, int> Hash;
  Hash hash;
  int nadded = 0;

  // Traverse each cell
  for(int i = 0; i < src->GetNumberOfCells(); i++)
    {
    if(regarr->GetTuple1(i) == region)
      {
      vtkCell *cell = src->GetCell(i);

      // Edges or faces?
      bool edges = cell->GetCellDimension() == 2;

      // Get all edges/faces in the cell
      int ne = edges ? cell->GetNumberOfEdges() : cell->GetNumberOfFaces();
      for(int k = 0; k < ne; k++)
        {
        vtkCell *edge = edges ? cell->GetEdge(k) : cell->GetFace(k);
        HashEntry he(edge->GetNumberOfPoints(), edge->GetPointIds()->GetPointer(0));
        if(hash.find(he) == hash.end())
          hash[he] = 1;
        else
          hash[he]++;
        }
      }
    }

  // Travese again
  for(int i = 0; i < src->GetNumberOfCells(); i++)
    {
    if(regarr->GetTuple1(i) == region)
      {
      vtkCell *cell = src->GetCell(i);

      // Edges or faces?
      bool edges = cell->GetCellDimension() == 2;

      // Get all edges/faces in the cell
      int ne = edges ? cell->GetNumberOfEdges() : cell->GetNumberOfFaces();
      for(int k = 0; k < ne; k++)
        {
        vtkCell *edge = edges ? cell->GetEdge(k) : cell->GetFace(k);
        HashEntry he(edge->GetNumberOfPoints(), edge->GetPointIds()->GetPointer(0));

        if(hash[he] == 1)
          {
          trg->InsertNextCell(edge->GetCellType(), 
            edge->GetNumberOfPoints(), edge->GetPointIds()->GetPointer(0));
          nadded++;
          }
        }
      }
    }

  return nadded;
}

template<class TMeshType>
void ConvertToMeshType(vtkUnstructuredGrid *src, vtkSmartPointer<TMeshType> &trg)
{}

template <>
void ConvertToMeshType(vtkUnstructuredGrid *src, vtkSmartPointer<vtkPolyData> &trg)
{
  vtkSmartPointer<vtkGeometryFilter> geo = vtkGeometryFilter::New();
  geo->SetInputData(src);
  geo->Update();
  trg = geo->GetOutput();
}

template <>
void ConvertToMeshType(vtkUnstructuredGrid *src, vtkSmartPointer<vtkUnstructuredGrid> &trg)
{
  trg = src;
}

template <class TMeshType> 
vtkSmartPointer<TMeshType> triangulate(TMeshType *mesh)
{
  throw MCException("Triangulation not supported");
}

template <>
vtkSmartPointer<vtkPolyData> triangulate(vtkPolyData *mesh)
{
  vtkSmartPointer<vtkTriangleFilter> filter = vtkSmartPointer<vtkTriangleFilter>::New();
  filter->SetInputData(mesh);
  filter->PassLinesOff();
  filter->PassVertsOff();
  filter->Update();
  vtkSmartPointer<vtkPolyData> result = filter->GetOutput();
  return result;
}

vnl_matrix<double> read_matrix(
    const std::string &fn, const std::string &desc, vnl_matrix<double> &mat)
{
  try
  {
    mat = vnl_file_matrix<double>(fn.c_str());
  }
  catch(...)
  {
    throw MCException("Error reading %s from %s", desc.c_str(), fn.c_str());
  }
}

vnl_matrix<double> read_vector(
    const std::string &fn, const std::string &desc, vnl_vector<double> &vec)
{
  try
  {
    vnl_matrix<double> mat = vnl_file_matrix<double>(fn.c_str());
    vec = mat.flatten_row_major();
  }
  catch(...)
  {
    throw MCException("Error reading %s from %s", desc.c_str(), fn.c_str());
  }
}

template <class TMeshType>
int meshcluster(Parameters &p, bool isPolyData)
{
  // Load design matrix, contrast vector and row mask from files
  vnl_matrix<double> mat, con;
  vnl_vector<double> row_mask;

  read_matrix(p.glm_design, "design matrix", mat);
  read_matrix(p.glm_contrast, "contrast vector", con);
  if(p.row_mask.size())
    read_vector(p.row_mask, "row mask", row_mask);
  else
    row_mask = vnl_vector<double>(mat.rows(), 1.0);

  if(p.row_mask.size() != mat.rows())
    throw MCException("Row mask size does not match design matrix dimensions");

  // Check the design matrix for missing values. Currently, we simply throw out the
  // rows with missing (NaN) values. In the future, we should implement a GLM with
  // missing data, or just move the whole package to R.
  // Note: this actually does not work!!
  std::vector<int> rows_kept;
  unsigned int n_rows_before_cleanup = mat.rows();
  for(unsigned int i = 0; i < mat.rows(); i++)
    {
    bool keep_row = row_mask[i] > 0.0;
    for(unsigned int j = 0; j < mat.cols(); j++)
      {
      if(vnl_math::isnan(mat(i,j)))
        {
        std::cout << "Discarding row " << i << " due to missing data" << std::endl;
        keep_row = false;
        }
      }
    if(keep_row)
      rows_kept.push_back(i);
    }

  // Now trim the design matrix
  if(rows_kept.size() < mat.rows())
    {
    vnl_matrix<double> mat_keep(rows_kept.size(), mat.cols());
    for(unsigned int i = 0; i < rows_kept.size(); i++)
      mat_keep.set_row(i, mat.get_row(rows_kept[i]));
    mat = mat_keep;
    }

  // Make sure matrices are compatible
  if(con.columns() != mat.columns())
    throw MCException("Matrix and contrast vector must have same number of columns");

  if(con.rows() != 1)
    throw MCException("Contrast vector must have one row");

  // Create the shuffle arrays
  vector<int> true_order, permutation;
  for(size_t i = 0; i < mat.rows(); i++)
    true_order.push_back(i);
  permutation = true_order;

  // If the threshold is negative, throw up
  if(p.threshold < 0)
    throw MCException("Negative threshold provided. This is not supported. You should invert the contrast instead");

  // Names for the statistics arrays
  string an_contrast = "Contrast", an_tstat = "T-statistic";
  string an_pval = "P-value (uncorr)", an_pcorr = "P-value (FWER corrected)";
  string an_beta = "Beta", an_res = "Residuals", an_betan = "Beta_Nuissance", an_dfarray="DOF";

  // Create an array of GLM objects
  vector<GeneralLinearModel *> glm;

  // Create an array of GLM objects for the reduced model and initialize flm variables
  vector <GeneralLinearModel *> glm_reduced;
  vnl_matrix <double> nuiss, nuiss_mat, nuiss_con;

  if(p.flag_freedman_lane)
    {
    // If using freedman_lane procedure for permutation testing
    if(p.fl_nuissance.size() == 0)
      {
      // check if string is empty, use zero entries in contrast to determine nuissance variables
      nuiss.set_size(1, mat.cols());
      nuiss.fill(0);
      for(unsigned int i = 0; i < con.cols(); i++)
        if(con[0][i] == 0)
          nuiss[0][i] = 1;
      }
    else
      {
      // if the user specified the nuissance vector
      read_matrix(p.fl_nuissance, "nuissance vector", nuiss);

      if(nuiss.columns() != mat.columns())
        throw MCException("Matrix and nuissance vector must have same number of columns");

      if(nuiss.rows() != 1)
        throw MCException("Nuissance vector must have one row");
      }

    // Create a design matrix containing only the nuissance variables
    nuiss_mat.set_size(mat.rows(), nuiss.array_one_norm());
    for(unsigned int i = 0, p = 0; i < nuiss.cols(); i++)
      {
      if(nuiss[0][i] == 1)
        nuiss_mat.set_column(p++, mat.get_column(i));
      }

    // Create a contarst matrix of the same size for consistency (not used in computation)
    nuiss_con.set_size(1, nuiss.array_one_norm());
    nuiss_con.fill(1);
    }

  // Read the meshes and initialize statistics arrays
  typedef vtkSmartPointer<TMeshType> MeshPointer;
  vector<MeshPointer> mesh;
  for(size_t i = 0; i < p.fn_mesh_input.size(); i++)
    {
    // Read mesh and optionally triangulate
    cout << "Reading mesh " << p.fn_mesh_input[i] << endl;
    vtkSmartPointer<TMeshType> mi = ReadMesh<TMeshType>(p.fn_mesh_input[i].c_str());
    if(p.flag_triangulate)
      {
      vtkSmartPointer<TMeshType> mi_tri = triangulate(mi.GetPointer());
      mesh.push_back(mi_tri);
      }
    else
      {
      mesh.push_back(mi);
      }

    // Get the data array from the mesh and error check
    vtkDataArray * data = GetArrayFromMesh(mesh[i], p.dom, p.array_name);
    if(!data)
      throw MCException("Array %s is missing in mesh %s",
                        p.array_name.c_str(), p.fn_mesh_input[i].c_str());
    if(data->GetNumberOfComponents() != (int) n_rows_before_cleanup)
      throw MCException("Wrong number of components (%d) in array %s in mesh %s. Should be %d.",
                        data->GetNumberOfComponents(), p.array_name.c_str(),
                        p.fn_mesh_input[i].c_str(), mat.rows());

    // Check if there is an exclusion array
    if(p.exclusion_array_name.length())
      {
      vtkDataArray *excl = GetArrayFromMesh(mesh[i], p.dom, p.exclusion_array_name);
      if(!excl)
        throw MCException("Array %s is missing in mesh %s",
                          p.exclusion_array_name.c_str(), p.fn_mesh_input[i].c_str());

      if(excl->GetNumberOfComponents() != (int) n_rows_before_cleanup)
        throw MCException("Wrong number of components (%d) in array %s in mesh %s. Should be %d.",
                          data->GetNumberOfComponents(), p.exclusion_array_name.c_str(),
                          p.fn_mesh_input[i].c_str(), mat.rows());

      if(data->GetNumberOfTuples() != excl->GetNumberOfTuples())
        throw MCException("Wrong number of tuples (%d) in array %s in mesh %s. Should be %d.",
                          data->GetNumberOfComponents(), p.exclusion_array_name.c_str(),
                          p.fn_mesh_input[i].c_str(), data->GetNumberOfTuples());

      unsigned int n_excl = 0;
      for(unsigned int j = 0; j < data->GetNumberOfTuples(); j++)
        for(int k = 0; k < data->GetNumberOfComponents(); k++)
          if(excl->GetComponent(j, k) != 0.0)
            {
            data->SetComponent(j, k, std::nan(0));
            ++n_excl;
            }

      cout << "  Converted " << n_excl << " values to NaN based on exclusion array" << endl;
      }

    // Clean the data if necessary
    if(rows_kept.size() < n_rows_before_cleanup)
      {
      vtkDataArray *rawData = data;
      char newName[1024];
      sprintf(newName,"MissingRows_%s", p.array_name.c_str());
      rawData->SetName(newName);

      data = AddArrayToMesh(mesh[i], p.dom, p.array_name, rows_kept.size(), 0, false);
      for(int i = 0; i < data->GetNumberOfTuples(); i++)
        {
        for(int j = 0; j < data->GetNumberOfComponents(); j++)
          {
          data->SetComponent(i, j, rawData->GetComponent(i, rows_kept[j]));
          }
        }
      }

    // Do diffusion if needed
    if(p.diffusion > 0)
      {
      if(p.dom == POINT)
        PointDataDiffusion(mesh[i], p.diffusion, p.delta_t, p.array_name.c_str());
      else
        CellDataDiffusion(mesh[i], p.diffusion, p.delta_t, p.array_name.c_str());
      }

    // Add the statistics array for output. The actual permutation testing always uses
    // CONTRAST or TSTAT, never the PVALUE (as that would waste time computing tcdf)
    vtkFloatArray * contrast = AddArrayToMesh(mesh[i], p.dom, an_contrast, 1, 0, p.ttype == CONTRAST);
    vtkFloatArray * tstat = AddArrayToMesh(mesh[i], p.dom, an_tstat, 1, 0, p.ttype != CONTRAST);
    vtkFloatArray * dfarray = AddArrayToMesh(mesh[i], p.dom, an_dfarray, 1, 0, false);
    vtkFloatArray * pval = AddArrayToMesh(mesh[i], p.dom, an_pval, 1, 0, false);
    vtkFloatArray * beta = AddArrayToMesh(mesh[i], p.dom, an_beta, mat.cols(), 0, false);
    vtkFloatArray * residual {0} ;
    vtkFloatArray * beta_nuiss {0};

    // Create new GLM
    glm.push_back(new GeneralLinearModel(data, contrast, tstat, pval, beta, residual, dfarray, mat, con));

    if(p.flag_freedman_lane)
      residual = AddArrayToMesh(mesh[i], p.dom, an_res, mat.rows(), NAN, false);
    beta_nuiss = AddArrayToMesh(mesh[i], p.dom, an_betan, nuiss_mat.cols(), 0, false);
    glm_reduced.push_back(new GeneralLinearModel(data, contrast, tstat, pval, beta_nuiss, residual, dfarray, nuiss_mat, nuiss_con));
    }

  // If the threshold is on the p-value, convert it to a threshold on T-statistic
  if(p.ttype == PVALUE)
    {
    p.ttype = TSTAT;
    if(p.threshold > 0)
      {
      double tthresh = glm.front()->PtoT(p.threshold);
      printf("P-value threshold %f converted to t-stat threshold %f\n", p.threshold, tthresh);
      p.threshold = tthresh;
      }
    }

  // Create an array of cluster computers
  vector<ClusterComputer *> clustcomp;

  if(p.np > 0)
    {
    printf("Executing GLM on %d random permutations\n", (int) p.np);
    if(p.threshold > 0)
      {
      for(size_t i = 0; i < mesh.size(); i++)
        clustcomp.push_back(new ClusterComputer(mesh[i], p.dom, p.threshold));
      }
    }
  // For Freedman-Lane procedure, need to perform GLM with the reduced model containing only nuissance variables
  if(p.flag_freedman_lane){

    for(size_t i = 0; i < mesh.size(); i++)
      {
      if(p.flag_missing_data)
        glm_reduced[i]->ComputeWithMissingData(glm_reduced[i]->Y,true, false, true );
      else
        glm_reduced[i]->Compute(glm_reduced[i]->Y,true, false, true);
      }
    }

  // Run permutation analysis
  vector<double> hArea(p.np), hPower(p.np), hStat(p.np);

  for(size_t ip = 0; ip < p.np; ip++)
    {

    // Initialize the histogram at zero
    hArea[ip] = 0; hPower[ip] = 0; hStat[ip] = 0;

    // Apply a random permutation to the labels array and update P matrix
    random_shuffle(permutation.begin(), permutation.end());
    
    // Build up the histogram of cluster areas (and powers)
    for(size_t i = 0; i < mesh.size(); i++)
      {
      vnl_matrix <double> Ytrue = glm[i]->Y;
      vnl_matrix<double> Yperm = Ytrue;

      // Shuffle the data according to the current permutation

      if(p.flag_freedman_lane) // Permute Y using nuissance model
        {
        // Read residual array from mesh
        vtkDataArray * res = GetArrayFromMesh(mesh[i], p.dom, an_res);
        vtkDataArray * beta_n = GetArrayFromMesh(mesh[i], p.dom, an_betan);
        vnl_matrix <double> res_mat;
        res_mat.set_size(Ytrue.rows(), Ytrue.cols());

        // Shuffle the residuals according to the current permutation
        for(size_t i = 0; i < res_mat.rows(); i++)
          for(int j = 0; j < res_mat.cols(); j++)
            res_mat(i,j) = res->GetComponent(j,i);

        vnl_matrix <double> res_matperm = res_mat;

        for (size_t i = 0; i < permutation.size(); i++)
          res_matperm.set_row(i,res_mat.get_row(permutation[i]));

        // Read the nuissance model beta values from the mesh
        vnl_matrix <double> beta_nuiss;
        beta_nuiss.set_size(nuiss_mat.cols(), Ytrue.cols());
        for(size_t i = 0; i < beta_nuiss.rows(); i++)
          for(int j = 0; j < beta_nuiss.cols(); j++)
            beta_nuiss(i,j) = beta_n->GetComponent(j,i);

        // Compute the permuted values

        Yperm = res_matperm + nuiss_mat*beta_nuiss;
        }
      else
        {
        for (size_t i = 0; i < permutation.size(); i++)
          Yperm.set_row(i,Ytrue.get_row(permutation[i]));
        }

      if(p.flag_missing_data){
        glm[i]->ComputeWithMissingData(Yperm, p.ttype != CONTRAST, p.ttype == PVALUE, false);
        }
      else{
        glm[i]->Compute(Yperm, p.ttype != CONTRAST, p.ttype == PVALUE, false);
        }

      // Generate a list of clusters (based on current scalars)
      if(p.threshold > 0)
        {
        ClusterArray ca = clustcomp[i]->ComputeClusters(false);

        // Now find the largest cluster
        for(size_t c = 0; c < ca.size(); c++)
          {
          if(ca[c].area > hArea[ip]) hArea[ip] = ca[c].area;
          if(ca[c].power > hPower[ip]) hPower[ip] = ca[c].power;
          }
        }

      // Likewise, find the largest scalar
      vtkDataArray *stat = GetScalarsFromMesh(mesh[i], p.dom);
      hStat[ip] = stat->GetMaxNorm();
      }

    // Some fancy formatting
    cout << "." << flush;
    if((ip+1) % 100 == 0 || (ip+1) == p.np)
      cout << " " << ip+1 << endl;
    } // end of permutation loop

  // Sort the histograms
  sort(hArea.begin(), hArea.end());
  sort(hPower.begin(), hPower.end());
  sort(hStat.begin(), hStat.end());

  // Going back to the original meshes, assign a cluster p-value to each mesh
  for(size_t i = 0; i < mesh.size(); i++)
    {
    vtkSmartPointer<TMeshType> mout;

    // Compute GLM (get all arrays, this is for keeps)
    if(p.flag_missing_data)
      glm[i]->ComputeWithMissingData(glm[i]->Y, true, true, false);
    else
      glm[i]->Compute(glm[i]->Y, true, true, false);

    // The rest is done only if permutations > 0
    if(p.np > 0)
      {
      if(p.threshold > 0)
        {
        // Generate true clusters (with full output)
        ClusterArray ca = clustcomp[i]->ComputeClusters(true);
        printf("MESH %s HAS %d CLUSTERS \n", p.fn_mesh_input[i].c_str(), (int) ca.size());

        if(ca.size() > 0)
          {
          ConvertToMeshType<TMeshType>(clustcomp[i]->GetFullMesh(), mout);

          // Assign a p-value to each cluster
          for(size_t c = 0; c < ca.size(); c++)
            {
            // Brute force search in the histogram :(
            size_t zArea = 0, zPower = 0;
            while(zArea < p.np && hArea[zArea] < ca[c].area) zArea++;
            while(zPower < p.np && hPower[zPower] < ca[c].power) zPower++;
            ca[c].pArea = 1.0 - zArea * 1.0 / p.np;
            ca[c].pPower = 1.0 - zPower * 1.0 / p.np;
            bool sig = (ca[c].pArea <= 0.05 || ca[c].pPower <= 0.05);
            printf("Cluster %03d:  AvgT = %6f; Area = %6f (p = %6f);  Power = %6f (p = %6f); %s\n",
                   (int) c, ca[c].tvalue/ca[c].n, ca[c].area, ca[c].pArea, ca[c].power, ca[c].pPower,
                   sig ? "***" : "");
            }

          // Create output mesh arrays for p-values
          string snArea = "p-cluster-area";
          vtkFloatArray * aArea = AddArrayToMesh(mout, CELL, snArea.c_str(), 1, 0);

          string snPower = "p-cluster-power";
          vtkFloatArray * aPower = AddArrayToMesh(mout, CELL, snPower.c_str(), 1, 0);

          // Get the cluster ID array
          vtkDataArray * aClusterId = mout->GetCellData()->GetArray("RegionId");

          // Set the mesh arrays' p-values
          for(int ic = 0; ic < mout->GetNumberOfCells(); ic++)
            {
            int r = (int) aClusterId->GetTuple1(ic);
            if(r >= 0)
              {
              // Set the cell array values
              aArea->SetTuple1(ic, ca[r].pArea);
              aPower->SetTuple1(ic, ca[r].pPower);
              }
            else
              {
              aArea->SetTuple1(ic, 1.0);
              aPower->SetTuple1(ic, 1.0);
              }
            }

          // If cluster edges requested, compute them
          if(p.flag_edges)
            {
            // Generate output
            TMeshType *emesh = TMeshType::New();
            emesh->SetPoints(mout->GetPoints());
            emesh->Allocate();

            // Set up arrays in the mesh
            vtkFloatArray * daArea = AddArrayToMesh(emesh, CELL, snArea.c_str(), 1, 0);
            vtkFloatArray * daPower = AddArrayToMesh(emesh, CELL, snPower.c_str(), 1, 0);
            vtkFloatArray * daClusterId = AddArrayToMesh(emesh, CELL, "RegionId", 1, 0);

            // Repeat for each cluster
            for(size_t c = 0; c < ca.size(); c++)
              {
              int added = AppendOpenEdgesFromMesh<TMeshType> (mout, aClusterId, c, emesh);
              for(int q = 0; q < added; q++)
                {
                daArea->InsertNextTuple1(ca[c].pArea);
                daPower->InsertNextTuple1(ca[c].pPower);
                daClusterId->InsertNextTuple1(c);
                }
              }

            // Figure out the name to save this to
            string path = vtksys::SystemTools::GetFilenamePath(p.fn_mesh_output[i]);
            string fnbase = vtksys::SystemTools::GetFilenameWithoutExtension(p.fn_mesh_output[i]);
            string ext = vtksys::SystemTools::GetFilenameExtension(p.fn_mesh_output[i]);

            string fnedge = p.fn_mesh_output[i];
            fnedge.replace(fnedge.length() - ext.length(), ext.length(), string("_edges") + ext);
            WriteMesh<TMeshType>(emesh, fnedge.c_str(), p.flag_write_binary);
            }
          }
        else
          {
          // No clusters found
          mout = mesh[i];
          }
        }
      else
        {
        // No threshold no clustering
        mout = mesh[i];
        }

      // Compute the corrected point-wise/cell-wise p-value for mout
      // printf("hStat: [%f %f %f %f %f]\n", hStat[0], hStat[p.np/4], hStat[p.np/2], hStat[3*p.np/4], hStat[p.np-1]);
      vtkDataArray * stat = GetScalarsFromMesh(mout, p.dom);
      vtkFloatArray * pcorr = AddArrayToMesh(mout, p.dom, an_pcorr, 1, 0);
      for(int j = 0; j < pcorr->GetNumberOfTuples(); j++)
        {
        double su = stat->GetTuple1(j);
        double pc = 1.0 - (lower_bound(hStat.begin(), hStat.end(), su) - hStat.begin()) / ((double) p.np);
        pcorr->SetTuple1(j, pc);
        }
      }
    else
      {
      // No permutation - just save the mesh on which GLM ran
      mout = mesh[i];
      }

    // Save the output mesh
    cout << "Saving output mesh.." << endl;
    WriteMesh<TMeshType>(mout, p.fn_mesh_output[i].c_str(), p.flag_write_binary);
    }
  return EXIT_SUCCESS;

}

int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();

  try
    {

    // Read command line
    Parameters p;

    for (int i = 1; i < argc; i++)
      {
      string arg = argv[i];

      if(arg == "-m" || arg == "--mesh")
        {
        p.fn_mesh_input.push_back(argv[++i]);
        p.fn_mesh_output.push_back(argv[++i]);
        }
      else if(arg == "-a" || arg == "--array")
        {
        p.array_name = argv[++i];
        }
      else if(arg == "-X")
        {
        p.exclusion_array_name = argv[++i];
        }
      else if(arg == "-p" || arg == "--perm")
        {
        p.np = atoi(argv[++i]);
        }
      else if(arg == "-g" || arg == "--glm")
        {
        p.glm_design = argv[++i];
        p.glm_contrast = argv[++i];
        }
      else if(arg == "-f" || arg == "--flane")
        {
        p.flag_freedman_lane = true;
        p.fl_nuissance = "";
        }
      else if(arg == "-n" || arg == "--nuiss")
        {
        p.fl_nuissance = argv[++i];
        }
      else if (arg == "-T")
        {
        p.flag_triangulate = true;
        }
      else if(arg == "-c" || arg == "--cell")
        {
        p.dom = CELL;
        }
      else if(arg == "-t" || arg == "--thresh")
        {
        p.threshold = atof(argv[++i]);
        }
      else if(arg == "-d" || arg == "--diffusion")
        {
        p.diffusion = atof(argv[++i]);
        }
      else if(arg == "--delta-t")
        {
        p.delta_t = atof(argv[++i]);
        }
      else if(arg == "-e" || arg == "--edges")
        {
        p.flag_edges = true;
        }
      else if(arg == "-M" || arg == "--missing")
        {
        p.flag_missing_data = true;
        p.flag_write_binary = true;
        }
      else if(arg == "-B" || arg == "--binary")
        {
        p.flag_write_binary = true;
        }
      else if(arg == "-s" || arg == "--stat")
        {
        string stat = argv[++i];
        if(stat == "T")
          p.ttype = TSTAT;
        else if(stat == "C")
          p.ttype = CONTRAST;
        else if(stat == "P")
          p.ttype = PVALUE;
        else
          throw MCException("Unknown parameter to --stat command (%s)", stat.c_str());
        }
      else throw MCException("Unknown command line switch %s", arg.c_str());
      }

    // Try reading the first mesh
    if (p.fn_mesh_input.size() == 0)
      throw MCException("No meshes specified on command line");

    // Read the meshes
    vtkDataReader *reader = vtkDataReader::New();
    reader->SetFileName(p.fn_mesh_input[0].c_str());
    reader->OpenVTKFile();
    reader->ReadHeader();

    if(reader->IsFileUnstructuredGrid())
      {
      reader->Delete();
      return meshcluster<vtkUnstructuredGrid>(p, false);
      }
    else if(reader->IsFilePolyData())
      {
      reader->Delete();
      return meshcluster<vtkPolyData>(p, true);
      }
    else
      {
      reader->Delete();
      throw MCException("Unsupported VTK data type in input file");
      }
    }
  catch (MCException &exc)
    {
    cerr << "Exception caught: " << exc.what() << endl; 
    return -1;
    }
}
