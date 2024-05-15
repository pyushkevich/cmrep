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
#include "VTKMeshHalfEdgeWrapper.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <random>
#include <deque>

#include <vtksys/SystemTools.hxx>


#include <exception>
#include <string>
#include <cstdarg>
#include <thread>
#include <mutex>
#include <functional>

const double glm_NAN = std::numeric_limits<double>::quiet_NaN();

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
    vsnprintf(buffer, 1024, fmt, parg );
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
  "  -M / --missing <min_obs|min_obs_fraction>\n"
  "                 Input data (per-vertex observations) have missing data encoded\n"
  "                 as NaNs. The algorithm will exclude these observations from the\n"
  "                 GLM and compute the GLM with only non-missing observations\n"
  "                 You must specify the minimum number of non-NaN observations at\n"
  "                 a vertex that are required to perform GLM at that vertex (e.g., 100),\n"
  "                 or a minimum fraction of observations that must be not NaN (e.g., 0.5).\n"
  "                 NOTE: on some operating systems NaNs can only be read from VTK\n"
  "                 meshes in binary (not ASCII) format. \n"
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
  "                 different subsets of subjects without generating separate input meshes \n"
  "  --tfce <del_h> Perform threshold-free cluster enhancement (TFCE). This will compute a \n"
  "                 TFCE map (Smith and Nichols, 2009) using the statistic selected with -s \n"
  "                 and step size delta. Parameter del_h (step size in cluster height) should be\n"
  "                 order of magnitude smaller than the range of the statistic, e.g., 0.1 for t-maps\n"
  "  --tfce-h <val> Set the value of the exponent applied to the cluster height, default 2.0\n"
  "  --tfce-e <val> Set the value of the exponent applied to the cluster extent, default 0.5\n"
  "  --threads <n>  Set the number of threads used in parallel\n"
  "  -z             Standardize the dependent variable and predictors before running GLM\n"
  "                 This way standardized betas can be obtained from the model\n";

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
TMeshType * ReadMesh(const char *)
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
  printf("Before diffusion: f[1000,10]=%8.4f\n", f->GetComponent(1000,10));
  unsigned int ncomp = f->GetNumberOfComponents();
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
      float *val_a = f_upd->GetPointer(ia * ncomp);
      float *val_b = f_upd->GetPointer(ib * ncomp);
      for(int j = 0; j < ncomp; j++)
        {
        if(!vnl_math::isnan(val_b[j]))
          val_a[j] += wa * (val_b[j] - val_a[j]);

        if(!vnl_math::isnan(val_a[j]))
          val_b[j] += wb * (val_a[j] - val_b[j]);
        }
      /*
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
      */
      }

    // Copy f_upd to f
    for(int i = 0; i < f->GetNumberOfTuples(); i++)
      for(int j = 0; j < f->GetNumberOfComponents(); j++)
        f->SetComponent(i, j, f_upd->GetComponent(i, j));

    cout << "." << flush;
    if((++jt) % 100 == 0 || t+dt >= time - 0.5 * dt)
      cout << " t = " << t+dt << endl;
    }
  printf("After diffusion: f[1000,10]=%8.4f\n", f->GetComponent(1000,10));
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
  printf("Before diffusion: f[1000,10]=%8.4f\n", f->GetComponent(1000,10));

  unsigned int ncomp = f->GetNumberOfComponents();
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
      float *val_a = f_upd->GetPointer(ia * ncomp);
      float *val_b = f_upd->GetPointer(ib * ncomp);
      for(int j = 0; j < ncomp; j++)
        {
        if(!vnl_math::isnan(val_b[j]))
          val_a[j] += wa * (val_b[j] - val_a[j]);

        if(!vnl_math::isnan(val_a[j]))
          val_b[j] += wb * (val_a[j] - val_b[j]);
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
  printf("After diffusion: f[1000,10]=%8.4f\n", f->GetComponent(1000,10));

}

template <class TDataset>
class TFCEComputer
{
public:
  TFCEComputer(TDataset *mesh, const char *arr_name, Domain dom, double dh, double E, double H)
  {
  }

  void Compute(vtkDataArray *tfce)
  {
    std::cerr << "TFCE not implemented for this mesh type" << std::endl;
    tfce->Fill(0.0);
  }
};

template <>
class TFCEComputer<vtkPolyData>
{
public:
  TFCEComputer(vtkPolyData *mesh, const char *arr_name, Domain dom, double dh, double E, double H)
    : mesh(mesh), he(mesh), dom(dom), dh(dh), E(E), H(H),
      triangle_queue(he.GetNumberOfFaces())
  {
    // Get the right array
    arr = dom == POINT ? mesh->GetPointData()->GetArray(arr_name) : mesh->GetCellData()->GetArray(arr_name);

    // Triangles are assigned a cluster index
    tri_cluster_idx.set_size(he.GetNumberOfFaces());

    // Initialize the per-triangle tfce map
    cluster_f.set_size(he.GetNumberOfFaces());
    tri_tfce.set_size(he.GetNumberOfFaces());

    // Compute the area of each triangle
    tri_area.set_size(he.GetNumberOfFaces());
    for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
      {
      // Compute the area of the triangle
      double p0[3], p1[3], p2[3];
      vtkCell *cell = mesh->GetCell(i);
      vtkIdType a0 = cell->GetPointId(0), a1 = cell->GetPointId(1), a2 = cell->GetPointId(2);
      mesh->GetPoint(a0, p0); mesh->GetPoint(a1, p1); mesh->GetPoint(a2, p2);
      tri_area[i] = std::abs(vtkTriangle::TriangleArea(p0, p1, p2));
      }

    // Initialize the partial area array
    tri_supra_thresh_area.set_size(he.GetNumberOfFaces());

    // Initialize the queue
    triangle_queue.clear();
  }

  void Compute(vtkDataArray *tfce)
  {
    // Reset the per-triangle TFCE map
    tri_tfce.fill(0.0);

    // Iterate over the levels
    for(double h = dh; true; h += dh)
      {
      // Mark triangles as visitable (0) or not visitable (1) based on whether they
      // are partially or fully above the threshold h
      tri_cluster_idx.fill(-1);
      unsigned int n_above = 0;

      // Initialize the suprathreshold area array to the full area array
      tri_supra_thresh_area.copy_in(tri_area.data_block());

      if(dom == POINT)
        {
        // Iterate over the half-edges. If one of the vertices in a half-edge is above
        // the threshold, then the face associated with this half-edge should be
        // assigned 0, i.e., it will be assigned to a cluster.
        // Also compute the suprathreshold area of each face. This can be computed by
        // scaling the area of the triangle by the proportion of each edge that is
        // above the threshold.
        for(unsigned int i = 0; i < he.GetNumberOfHalfEdges(); i++)
          {
          int v1 = he.GetHalfEdgeVertex(i), v0 = he.GetHalfEdgeTailVertex(i);
          int f = he.GetHalfEdgeFace(i);
          double h0 = arr->GetTuple1(v0), h1 = arr->GetTuple1(v1);
          if(h0 >= h || h1 >= h)
            {
            tri_cluster_idx[he.GetHalfEdgeFace(i)] = 0;
            n_above++;
            if(h0 < h)
              {
              double t = (h - h0) / (h1 - h0);
              tri_supra_thresh_area[f] *= 1 - t;
              }
            else if(h1 < h)
              {
              double t = (h - h0) / (h1 - h0);
              tri_supra_thresh_area[f] *= t;
              }
            }
          }
        }
      else
        {
        // Just check each face - no need to check vertices
        for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
          if(arr->GetTuple1(i) >= h)
            {
            tri_cluster_idx[i] = 0;
            n_above++;
            }
        }

      // Zero out the area of faced that are not visited
      for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
        if(tri_cluster_idx[i] < 0)
          tri_supra_thresh_area[i] = 0.0;

      // Check if any triangles are above the threshold, if not, quit
      if(n_above == 0)
        break;

      // Initialize the per-cluster function f (h^H * area^E) for each cluster
      cluster_f.fill(0.0);

      // Now we need to perform fast marching over the triangles, grouping them into
      // clusters and accumulating cluster areas. For this, we need some kind of a heap
      // to keep track of all the "active" triangles, i.e., triangles that could be
      // added to the current cluster.
      int curr_cluster = 1;
      double area = 0;
      triangle_queue.clear();
      for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
        {
        // If a triangle is not on the candidate list, skip it, it has already been visited
        // or it should not be visited (-1)
        if(tri_cluster_idx[i] != 0) continue;

        // Add the current vertex to the queue
        triangle_queue.push_back(i);
        while(triangle_queue.size() > 0)
          {
          // Look at the front of the queue, mark as part of current cluster, add area
          unsigned int f = triangle_queue.front();
          tri_cluster_idx[f] = curr_cluster;
          area += tri_supra_thresh_area[f];

          // Look at the neighbors of this triangle, add them to the queue if not visited
          int he_start = he.GetFaceHalfEdge(f), he_curr = he_start;
          for(unsigned int j = 0; j < 3; j++)
            {
            unsigned int f_nbr = he.GetHalfEdgeFace(he.GetHalfEdgeOpposite(he_curr));
            if(f_nbr == f)
              throw std::string("wrong logic");
            if(f_nbr != VTKMeshHalfEdgeWrapper::NA && tri_cluster_idx[f_nbr] == 0)
              {
              tri_cluster_idx[f_nbr] = curr_cluster;
              triangle_queue.push_back(f_nbr);
              }

            he_curr = he.GetHalfEdgeNext(he_curr);
            }

          // Pop the current triangle off the queue
          triangle_queue.pop_front();
          }

        // Record the area of this cluster -> send it to the vertices.
        cluster_f[curr_cluster] = std::pow(h, H) * std::pow(area, E);
        // printf("h = %8.6f  cluster %03d  area = %8.6f  f = %8.6f\n", h, curr_cluster, area, cluster_f[curr_cluster]);

        // Increase the cluster number, reset the area
        curr_cluster++;
        area = 0.0;
        }

      // Iterate over the triangles. Each triangle should be assigned a cluster,
      // except for triangles that are below the current height h, which will be
      // given cluster value -1. We initialize the cluster array to 0, which means
      // no triangle has been assigned yet
      for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
        if(tri_cluster_idx[i] >= 0)
          tri_tfce[i] += cluster_f[tri_cluster_idx[i]];

      // Generate an output mesh
      /*
      char buffer[256];
      snprintf(buffer, 256, "tfve_cluster_index_%6.4f", h);
      vtkSmartPointer<vtkFloatArray> a_idx = AddArrayToMesh(mesh, CELL, buffer, 1, 0.0);
      snprintf(buffer, 256, "tfve_cluster_parea_%6.4f", h);
      vtkSmartPointer<vtkFloatArray> a_area = AddArrayToMesh(mesh, CELL, buffer, 1, 0.0);
      snprintf(buffer, 256, "tfve_cluster_aratio_%6.4f", h);
      vtkSmartPointer<vtkFloatArray> a_area_ratio = AddArrayToMesh(mesh, CELL, buffer, 1, 0.0);
      snprintf(buffer, 256, "tfve_cluster_f_%6.4f", h);
      vtkSmartPointer<vtkFloatArray> a_f = AddArrayToMesh(mesh, CELL, buffer, 1, 0.0);
      for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
        {
        a_idx->SetTuple1(i, tri_cluster_idx[i]);
        a_area->SetTuple1(i, tri_supra_thresh_area[i]);
        a_area_ratio->SetTuple1(i, tri_supra_thresh_area[i] / tri_area[i]);
        a_f->SetTuple1(i, cluster_f[tri_cluster_idx[i]]);
        }
        */
      }

    // The final TFCE array is on cells - it does not really make so much sense to
    // project back to vertex space based on how we compute it
    tfce->SetNumberOfComponents(1);
    tfce->SetNumberOfTuples(he.GetNumberOfFaces());
    for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
      tfce->SetTuple1(i, tri_tfce[i]);

    /*
    // Assign the tfce value to vertices or cells depending on what the use wants
    if(dom == POINT)
      {
      tfce->SetNumberOfComponents(1);
      tfce->SetNumberOfTuples(he.GetNumberOfVertices());
      tfce->Fill(0.0);
      for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
        {
        // Split the tFCE between the vertices
        vtkCell *cell = mesh->GetCell(i);
        for(unsigned int j = 0; j < 3; j++)
          {
          auto pid = cell->GetPointId(j);
          tfce->SetTuple1(pid, tfce->GetTuple1(pid) + tri_tfce[i] / 3.0);
          }
        }
      }
    else
      {
      tfce->SetNumberOfComponents(1);
      tfce->SetNumberOfTuples(he.GetNumberOfFaces());
      for(unsigned int i = 0; i < he.GetNumberOfFaces(); i++)
        tfce->SetTuple1(i, tri_tfce[i]);
      }
    */
  }

private:
  vtkPolyData *mesh;
  VTKMeshHalfEdgeWrapper he;
  Domain dom;
  double dh, E, H;

  // Active data array
  vtkDataArray *arr;

  // For every triangle in the mesh, we need to keep track of a cluster index
  // assigned to it.
  vnl_vector<int> tri_cluster_idx;
  vnl_vector<double> tri_area, tri_supra_thresh_area;

  // For every vertex, associate it with a cluster
  vnl_vector<double> cluster_f, tri_tfce;

  // A queue used to manage the current cluster
  std::deque<int> triangle_queue;
};

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
  vtkSmartPointer<vtkFloatArray> scalar_backup;

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
    fThresh->SetUpperThreshold(thresh);
    fThresh->SetThresholdFunction(vtkThreshold::THRESHOLD_UPPER);
    fConnect->SetInputConnection(fThresh->GetOutputPort());

    fThreshInv = vtkThreshold::New();
    fThreshInv->SetInputData(mesh);
    fThreshInv->SetInputArrayToProcess(0, 0, 0, 
      vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
    fThreshInv->SetLowerThreshold(thresh);
    fThreshInv->SetThresholdFunction(vtkThreshold::THRESHOLD_LOWER);
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
  ClusterArray ca(0);
  if(nelt > 0)
    {
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
    ca = ClusterArray(fConnect->GetNumberOfExtractedRegions());
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
    }

  // Return the cluster array
  return ca;
}

class GeneralLinearModel
{
public:
  GeneralLinearModel(
    vtkFloatArray *data, vtkFloatArray *contrast,
    vtkFloatArray *tstat, vtkFloatArray *pval, vtkFloatArray *beta, vtkFloatArray *residual, vtkFloatArray *dfarray,
    const vnl_matrix<double> &design_matrix,
    const vnl_matrix<double> &contrast_vector);

  // Make a deep copy of another GLM
  GeneralLinearModel(const GeneralLinearModel *src);

  void Compute(const vnl_matrix<double> &Yperm, bool need_t, bool need_p, bool need_res);
  void ComputeWithMissingData(const vnl_matrix<double> &Yperm, bool need_t, bool need_p, bool need_res);
  //int compute_rank(vnl_matrix <double> D); 
  double PtoT(double p);
  vnl_matrix<double> Y; // data matrix for permutations

private:
  void CommonInit();
  vtkSmartPointer<vtkFloatArray> DeepCopyArray(vtkFloatArray *src);

  // Pointers to the data arrays
  vtkSmartPointer<vtkFloatArray> data, contrast, tstat, pval, beta, residual, dfarray;

  // Copies of the matrix, other data
  vnl_matrix<double> X,cv;
  int nelt, rank, df;
};

GeneralLinearModel::GeneralLinearModel(
    vtkFloatArray *data, vtkFloatArray *contrast,
    vtkFloatArray *tstat, vtkFloatArray *pval, vtkFloatArray *beta,
    vtkFloatArray *residual, vtkFloatArray *dfarray,
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
  this->CommonInit();
  }

GeneralLinearModel::GeneralLinearModel(const GeneralLinearModel *src)
{
  // Copy input
  this->data = DeepCopyArray(src->data);
  this->contrast = DeepCopyArray(src->contrast);
  this->tstat = DeepCopyArray(src->tstat);
  this->pval = DeepCopyArray(src->pval);
  this->beta = DeepCopyArray(src->beta);
  this->residual = DeepCopyArray(src->residual);
  this->dfarray = DeepCopyArray(src->dfarray);
  this->X = src->X;
  this->cv = src->cv;
  this->CommonInit();
}

vtkSmartPointer<vtkFloatArray>
GeneralLinearModel
::DeepCopyArray(vtkFloatArray *src)
{
  if(src)
    {
    vtkNew<vtkFloatArray> copy;
    copy->DeepCopy(src);
    return copy;
    }
  else return nullptr;
}

void
GeneralLinearModel::Compute(const vnl_matrix<double> &Yperm, bool need_t, bool need_p, bool need_res)
  {

  // Shuffle the rows according to current permutan_resr (size_t i = 0; i < permutation.size(); i++)
  //  for (size_t j = 0; j < X.cols(); j++)
  //    	Xperm[i][j] = X[permutation[i]][j];

  // Standard GLM calculation
  vnl_matrix<double> A = vnl_matrix_inverse<double>(X.transpose() * X).pinverse(rank);
  vnl_matrix<double> bhat = (A * X.transpose()) * Yperm;
  
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
    vnl_matrix<double> errmat = Yperm - X * bhat;
    
    // Residual variance
    vnl_matrix<double> resvar(1, nelt);
    resvar.fill(0.0);
    for(int j = 0; j < nelt; j++)
      {
      for(size_t i = 0; i < X.rows(); i++)
        resvar(0, j) += errmat(i,j) * errmat(i,j);
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

template <typename TDataArray> void set_tuple_to_nan(TDataArray *arr, int j)
{
  for(int k = 0; k < arr->GetNumberOfComponents(); k++)
    arr->SetComponent(j, k, glm_NAN);
}

template <typename TDataArray> void set_tuple_to_zero(TDataArray *arr, int j)
{
  for(int k = 0; k < arr->GetNumberOfComponents(); k++)
    arr->SetComponent(j, k, 0);
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

    // If the number of nans equals the number of observations, we cannot do any stats and we
    // should just assign NaN to everything in the output
    if(n_nans == ns)
      {
      set_tuple_to_nan(beta.GetPointer(), j);
      set_tuple_to_nan(contrast.GetPointer(), j);
      if(need_t)
        {
        set_tuple_to_nan(tstat.GetPointer(), j);
        dfarray->SetTuple1(j, 0);
        if(residual)
          set_tuple_to_nan(residual.GetPointer(), j);

        if(need_p)
          set_tuple_to_nan(pval.GetPointer(), j);
        }
      continue;
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

void GeneralLinearModel::CommonInit()
{
  // Initialize other stuff
  nelt = data->GetNumberOfTuples();

  // Rank and degrees of freedom
  // rank = vnl_rank(X.transpose() * X, vnl_rank_row);
  rank = X.columns();
  df = X.rows() - rank;

  // Y matrix
  Y.set_size(X.rows(), nelt);
  for(size_t i = 0; i < X.rows(); i++)
    for(int j = 0; j < nelt; j++)
      Y(i,j) = data->GetComponent(j,i);
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

  // Minimum number of valid (non-missing) observations that must be present
  // at a vertex for it to be included in the GLM. Specified either as a number
  // of observations or as a fraction.
  double min_valid_obs;

  // Whether to write output in binary format
  bool flag_write_binary;

  // Whether to use the freedman-lane procedure
  bool flag_freedman_lane;

  // Triangulate?
  bool flag_triangulate;

  // What is the threshold on
  ThresholdType ttype;

  // TFCE config
  double tfce_delta_h = 0.0, tfce_H = 2.0, tfce_E = 0.5;

  // Z-transform
  bool flag_z_transform = false;

  // Max threads
  int max_threads;

  Parameters()
    {
    np = 0; dom = POINT; threshold = 0; diffusion = 0; flag_edges = false; ttype = TSTAT; delta_t = 0.01; 
    flag_missing_data = false; min_valid_obs = 0.0; flag_write_binary = false; flag_freedman_lane = false; flag_triangulate = false;
    max_threads = 0;
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

void read_matrix(const std::string &fn, const std::string &desc, vnl_matrix<double> &mat)
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

void read_vector(const std::string &fn, const std::string &desc, vnl_vector<double> &vec)
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

void copy_array_replace_nans(vtkDataArray *src, vtkDataArray *trg, float nan_value = 0.0)
{
  for(int j = 0; j < src->GetNumberOfTuples(); j++)
    {
    for(int k = 0; k < src->GetNumberOfComponents(); k++)
      {
      float c = src->GetComponent(j,k);
      trg->SetComponent(j,k,std::isnan(c) ? nan_value : c);
      }
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

  if(row_mask.size() != mat.rows())
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

  // Standardize the design matrix if requested
  if(p.flag_z_transform)
    {
    for(size_t i = 0; i < mat.columns(); i++)
      {
      int n = 0;
      double sum = 0, sum_sq = 0;
      std::set<double> uniq_val;
      for(size_t j = 0; j < mat.rows(); j++)
        {
        double x = mat(j,i);
        if(!std::isnan(x))
          {
          sum += x;
          sum_sq += x * x;
          n++;
          uniq_val.insert(x);
          }
        }
      if(uniq_val.size() > 2)
        {
        // This is not a categorical variable
        double x_bar = sum / n;
        double var = sum_sq / n - x_bar * x_bar;
        double sq = sqrt(var);
        printf("Standardized predictor %02d, µ = %6.4f  σ = %6.4f\n", (int) i, x_bar, sq);
        for(size_t j = 0; j < mat.rows(); j++)
          {
          double x = mat(j,i);
          if(!std::isnan(x))
            mat(j,i) = (x - x_bar) / sq;
          }
        }
      }
    }

  // Create the shuffle arrays
  vector<int> true_order;
  for(size_t i = 0; i < mat.rows(); i++)
    true_order.push_back(i);

  // If the threshold is negative, throw up
  if(p.threshold < 0)
    throw MCException("Negative threshold provided. This is not supported. You should invert the contrast instead");

  // Names for the statistics arrays
  string an_contrast = "Contrast", an_tstat = "T-statistic";
  string an_pval = "P-value (uncorr)", an_pcorr = "P-value (FWER corrected)", an_fdr = "P-value (FDR corrected)";
  string an_beta = "Beta", an_res = "Residuals", an_betan = "Beta_Nuissance", an_dfarray="DOF";
  string an_cluster_array = "ClusterVariable";
  string an_tfce = "tfce", an_tfce_pcorr = "tfce P-value (FWER corrected)";

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
    vtkFloatArray *data = vtkArrayDownCast<vtkFloatArray>(GetArrayFromMesh(mesh[i], p.dom, p.array_name));
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
            data->SetComponent(j, k, glm_NAN);
            ++n_excl;
            }

      cout << "  Converted " << n_excl << " values to NaN based on exclusion array" << endl;
      }

    // Clean the data if necessary
    if(rows_kept.size() < n_rows_before_cleanup)
      {
      vtkDataArray *rawData = data;
      char newName[1024];
      snprintf(newName,1024,"MissingRows_%s", p.array_name.c_str());
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

    // If missing data is specified, check for number of nans at each vertex, and if the number
    // is too low, replace all values in that column with NaNs
    if(p.min_valid_obs > 0.0)
      {
      int min_not_nan = p.min_valid_obs < 1.0
                        ? int(data->GetNumberOfComponents() * p.min_valid_obs)
                        : int(p.min_valid_obs);
      int n_excluded = 0;
      for(unsigned int j = 0; j < data->GetNumberOfTuples(); j++)
        {
        int n_nans = 0;
        for(int k = 0; k < data->GetNumberOfComponents(); k++)
          if(std::isnan(data->GetComponent(j, k)))
            n_nans++;
        if(data->GetNumberOfComponents() - n_nans < min_not_nan)
          {
          set_tuple_to_nan(data, j);
          n_excluded++;
          }
        }
      cout << "Excluded " << n_excluded << " vertices with fewer than " << min_not_nan << " valid observations" << endl;
      }

    // Perform z-transformation
    if(p.flag_z_transform)
      {
      // Standardize the dependent variable
      int n_standardized = 0;
      for(unsigned int j = 0; j < data->GetNumberOfTuples(); j++)
        {
        int n = 0;
        double sum = 0, sum_sq = 0;
        for(int k = 0; k < data->GetNumberOfComponents(); k++)
          {
          double x = data->GetComponent(j, k);
          if(!std::isnan(x))
            {
            sum += x;
            sum_sq += x * x;
            n++;
            }
          }
        if(n > 1)
          {
          n_standardized++;
          double x_bar = sum / n;
          double var = sum_sq / n - x_bar * x_bar;
          double sq = sqrt(var);
          for(int k = 0; k < data->GetNumberOfComponents(); k++)
            {
            double x = data->GetComponent(j, k);
            if(!std::isnan(x))
              data->SetComponent(j, k, (x - x_bar) / sq);
            }
          }
        }
      printf("Standardized %d of %d dependent variables\n", n_standardized, (int) data->GetNumberOfTuples());
      }

    // Add the statistics array for output. The actual permutation testing always uses
    // CONTRAST or TSTAT, never the PVALUE (as that would waste time computing tcdf)
    vtkFloatArray * contrast = AddArrayToMesh(mesh[i], p.dom, an_contrast, 1, 0, false);
    vtkFloatArray * tstat = AddArrayToMesh(mesh[i], p.dom, an_tstat, 1, 0, false);
    vtkFloatArray * dfarray = AddArrayToMesh(mesh[i], p.dom, an_dfarray, 1, 0, false);
    vtkFloatArray * pval = AddArrayToMesh(mesh[i], p.dom, an_pval, 1, 0, false);
    vtkFloatArray * fdr = AddArrayToMesh(mesh[i], p.dom, an_fdr, 1, 0, false);
    vtkFloatArray * beta = AddArrayToMesh(mesh[i], p.dom, an_beta, mat.cols(), 0, false);
    vtkFloatArray * tfce_array = AddArrayToMesh(mesh[i], CELL, an_tfce, 1, 0, false);
    vtkFloatArray * cluster_array = AddArrayToMesh(mesh[i], p.dom, an_cluster_array, 1, 0, true);
    vtkFloatArray * residual {0} ;
    vtkFloatArray * beta_nuiss {0};

    // Create new GLM
    glm.push_back(new GeneralLinearModel(data, contrast, tstat, pval, beta, residual, dfarray, mat, con));

    // For FL, create a reduced GLM
    if(p.flag_freedman_lane)
      {
      residual = AddArrayToMesh(mesh[i], p.dom, an_res, mat.rows(), NAN, false);
      beta_nuiss = AddArrayToMesh(mesh[i], p.dom, an_betan, nuiss_mat.cols(), 0, false);
      glm_reduced.push_back(new GeneralLinearModel(data, contrast, tstat, pval, beta_nuiss, residual, dfarray, nuiss_mat, nuiss_con));
      }
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
  vector<double> hArea(p.np), hPower(p.np), hStat(p.np), hTFCE(p.np);

  // Name of the array that is involved in thresholding for permutation testing
  string an_ttype;
  switch(p.ttype)
    {
    case CONTRAST: an_ttype = an_contrast; break;
    case TSTAT: an_ttype = an_tstat; break;
    case PVALUE: an_ttype = an_pval; break;
    }

  // Perform permutation test
  if(p.np > 0)
    {
    // Parallelize over threads
    int nthreads = std::thread::hardware_concurrency();
    if(p.max_threads > 0)
      nthreads = std::min(nthreads, p.max_threads);

    std::vector<std::thread> threads;
    printf("Executing GLM on %d random permutations using %d threads\n", (int) p.np, (int) nthreads);

    // Mutex for printing dots
    std::mutex critical;
    int n_done = 0;

    // For each thread, give it a set of permutations to run
    for(int it = 0; it < nthreads; it++)
      {
      int ip_start = ceil(p.np / nthreads) * it;
      int ip_end = std::min((int) p.np, (int) ceil(p.np / nthreads) * (it + 1));
      threads.push_back(std::thread(std::bind([&](const int ip_start, const int ip_end)
        {
        // Create thread's copy of permutations
        vector<int> permutation = true_order;

        // We need to create a copy of the meshes for this thread because the GLM
        // class updates arrays in the mesh
        vector<MeshPointer> mesh_t;
        vector<GeneralLinearModel *> glm_t;
        vector<ClusterComputer *> clustcomp_t;
        vector<TFCEComputer<TMeshType> *> tfcecomp_t;
        for(auto &m : mesh)
          {
          // Create a copy of the mesh
          vtkNew<TMeshType> m_copy;
          m_copy->DeepCopy(m);
          mesh_t.push_back(m_copy);

          // Create a copy of the GLM
          glm_t.push_back(new GeneralLinearModel(
                            vtkArrayDownCast<vtkFloatArray>(GetArrayFromMesh(m_copy, p.dom, p.array_name)),
                            vtkArrayDownCast<vtkFloatArray>(GetArrayFromMesh(m_copy, p.dom, an_contrast)),
                            vtkArrayDownCast<vtkFloatArray>(GetArrayFromMesh(m_copy, p.dom, an_tstat)),
                            vtkArrayDownCast<vtkFloatArray>(GetArrayFromMesh(m_copy, p.dom, an_pval)),
                            vtkArrayDownCast<vtkFloatArray>(GetArrayFromMesh(m_copy, p.dom, an_beta)),
                            nullptr,
                            vtkArrayDownCast<vtkFloatArray>(GetArrayFromMesh(m_copy, p.dom, an_dfarray)),
                            mat, con));

          // Create a cluster computer
          clustcomp_t.push_back(new ClusterComputer(m_copy, p.dom, p.threshold));

          // Create a TFCE computer
          if(p.tfce_delta_h > 0.0)
            {
            tfcecomp_t.push_back(new TFCEComputer<TMeshType>(m_copy, an_ttype.c_str(), p.dom, p.tfce_delta_h, p.tfce_E, p.tfce_H));
            }
          }

        // Perform permutations for this thread
        for(unsigned int ip = ip_start; ip < ip_end; ip++)
          {
          // Shuffle the permutation
          std::mt19937 rng(std::time(nullptr));
          std::shuffle(permutation.begin(), permutation.end(), rng);

          // Initialize the histogram at zero
          hArea[ip] = 0; hPower[ip] = 0; hStat[ip] = 0; hTFCE[ip] = 0;

          // Build up the histogram of cluster areas (and powers)
          for(size_t i = 0; i < mesh_t.size(); i++)
            {
            vnl_matrix <double> Ytrue = glm_t[i]->Y;
            vnl_matrix<double> Yperm = Ytrue;

            if(p.flag_freedman_lane) // Permute Y using nuissance model
              {
              // Read residual array from mesh
              vtkDataArray * res = GetArrayFromMesh(mesh_t[i], p.dom, an_res);
              vtkDataArray * beta_n = GetArrayFromMesh(mesh_t[i], p.dom, an_betan);
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
              Yperm = res_matperm + nuiss_mat * beta_nuiss;
              }
            else
              {
              for (size_t i = 0; i < permutation.size(); i++)
                Yperm.set_row(i,Ytrue.get_row(permutation[i]));
              }

            if(p.flag_missing_data)
              glm_t[i]->ComputeWithMissingData(Yperm, p.ttype != CONTRAST, p.ttype == PVALUE, false);
            else
              glm_t[i]->Compute(Yperm, p.ttype != CONTRAST, p.ttype == PVALUE, false);

            // Generate a list of clusters (based on current scalars)
            if(p.threshold > 0)
              {
              copy_array_replace_nans(GetArrayFromMesh(mesh_t[i], p.dom, an_ttype),
                                      GetArrayFromMesh(mesh_t[i], p.dom, an_cluster_array), 0);
              ClusterArray ca = clustcomp_t[i]->ComputeClusters(false);

              // Now find the largest cluster
              for(size_t c = 0; c < ca.size(); c++)
                {
                if(ca[c].area > hArea[ip]) hArea[ip] = ca[c].area;
                if(ca[c].power > hPower[ip]) hPower[ip] = ca[c].power;
                }
              }

            // Compute TFCE if requested
            if(p.tfce_delta_h > 0.0)
              {
              // Compute TFCE
              vtkDataArray *arr_tfce = GetArrayFromMesh(mesh[i], CELL, an_tfce);
              tfcecomp_t[i]->Compute(arr_tfce);
              double tfce_stat_max = arr_tfce->GetMaxNorm();
              if(tfce_stat_max > hTFCE[ip])
                hTFCE[ip] = tfce_stat_max;
              }

            // Likewise, find the largest scalar
            vtkDataArray *stat = GetScalarsFromMesh(mesh_t[i], p.dom);
            double stat_max = stat->GetMaxNorm();
            if(stat_max > hStat[ip])
              hStat[ip] = stat_max;
            }

          std::lock_guard<std::mutex> guard(critical);
          n_done++;
          cout << "." << flush;
          if(n_done % 100 == 0 || n_done == p.np)
            cout << " " << n_done << endl;
          }
        }, ip_start, ip_end))); // End of in-thread permutation loop
      } // Threads have been created

    // Run the threads
    std::for_each(threads.begin(),threads.end(),[](std::thread& x){x.join();});

    // Sort the histograms
    sort(hArea.begin(), hArea.end());
    sort(hPower.begin(), hPower.end());
    sort(hStat.begin(), hStat.end());
    sort(hTFCE.begin(), hTFCE.end());

    } // End permutation conditional


  /*
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
    */

  // --------------------------------
  // Perform GLM on unpermuted meshes
  // --------------------------------

  // Going back to the original meshes, assign a cluster p-value to each mesh
  for(size_t i = 0; i < mesh.size(); i++)
    {
    vtkSmartPointer<TMeshType> mout;

    // Compute GLM (get all arrays, this is for keeps)
    if(p.flag_missing_data)
      glm[i]->ComputeWithMissingData(glm[i]->Y, true, true, false);
    else
      glm[i]->Compute(glm[i]->Y, true, true, false);
    }


  // ------------------------------------
  // Compute TFCE for the original meshes
  // ------------------------------------
  if(p.tfce_delta_h > 0.0)
    {
    for(size_t i = 0; i < mesh.size(); i++)
      {
      // Create a TFCE computer for this mesh
      TFCEComputer<TMeshType> tfce(mesh[i], an_ttype.c_str(), p.dom, p.tfce_delta_h, p.tfce_E, p.tfce_H);

      // Compute TFCE
      vtkDataArray *arr_tfce = GetArrayFromMesh(mesh[i], CELL, an_tfce);
      tfce.Compute(arr_tfce);
      }
    }


  // ------------------------------
  // Compute FDR corrected p-values
  // ------------------------------

  // Create a sorted array of all the p-values
  using PValueWithIndex = std::tuple<double, unsigned int, unsigned int, double>;
  std::vector<PValueWithIndex> p_sorted;
  std::vector<vtkDataArray *> fdr_arr;
  for(unsigned int i = 0; i < mesh.size(); i++)
    {
    vtkDataArray *p_val = GetArrayFromMesh(mesh[i], p.dom, an_pval);
    fdr_arr.push_back(GetArrayFromMesh(mesh[i], p.dom, an_fdr));
    for(unsigned int j = 0; j < p_val->GetNumberOfTuples(); j++)
      {
      double p_uncorr = p_val->GetTuple1(j);
      if(!std::isnan(p_uncorr))
        p_sorted.push_back(std::make_tuple(p_uncorr, i, j, 1.0));
      else
        set_tuple_to_nan(fdr_arr.back(), j);
      }
    }

  // Sort the p-values
  std::sort(p_sorted.begin(), p_sorted.end());

  // Perform adjustment using the BH procedure
  for(unsigned int k = 0; k < p_sorted.size(); k++)
    {
    double p_uncorr = std::get<0>(p_sorted[k]);
    double p_scaled = p_uncorr / (k * 1.0 / p_sorted.size());
    std::get<3>(p_sorted[k]) = p_scaled;
    }

  // Ensure that the p-values are in non-descending order
  for(unsigned int k = p_sorted.size()-1; k > 0; k--)
    {
    if(std::get<3>(p_sorted[k-1]) > std::get<3>(p_sorted[k]))
      std::get<3>(p_sorted[k-1]) = std::get<3>(p_sorted[k]);
    }

  // Ensure that the p-values are less than one
  for(unsigned int k = 0; k < p_sorted.size(); k++)
    std::get<3>(p_sorted[k]) = std::min(1.0, std::get<3>(p_sorted[k]));

  for(unsigned int k = 0; k < p_sorted.size(); k++)
    {
    double p_uncorr = std::get<0>(p_sorted[k]);
    double p_fdr = std::get<3>(p_sorted[k]);
    unsigned int i_mesh = std::get<1>(p_sorted[k]), i_elt = std::get<2>(p_sorted[k]);
    // printf("Tuple %f, %f, %d, %d   thresh %f \n", p_uncorr, p_fdr, i_mesh, i_elt, k * 0.05 / p_sorted.size());
    fdr_arr[i_mesh]->SetTuple1(i_elt, p_fdr);
    }

  // -------------------------------
  // Generate and save output meshes
  // -------------------------------

  // Going back to the original meshes, assign a cluster p-value to each mesh
  for(size_t i = 0; i < mesh.size(); i++)
    {
    vtkSmartPointer<TMeshType> mout;

    // The rest is done only if permutations > 0
    if(p.np > 0)
      {
      if(p.threshold > 0)
        {
        // Generate true clusters (with full output)
        ClusterComputer clustcomp(mesh[i], p.dom, p.threshold);
        copy_array_replace_nans(GetArrayFromMesh(mesh[i], p.dom, an_ttype),
                                GetArrayFromMesh(mesh[i], p.dom, an_cluster_array), 0);
        ClusterArray ca = clustcomp.ComputeClusters(true);
        printf("MESH %s HAS %d CLUSTERS \n", p.fn_mesh_input[i].c_str(), (int) ca.size());

        if(ca.size() > 0)
          {
          ConvertToMeshType<TMeshType>(clustcomp.GetFullMesh(), mout);

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

      // Compute the corrected TFCE p-value
      if(p.tfce_delta_h > 0.0)
        {
        vtkDataArray *tfce = GetArrayFromMesh(mout, CELL, an_tfce);
        vtkFloatArray *tfce_pcorr = AddArrayToMesh(mout, CELL, an_tfce_pcorr, 1, 0);
        for(int j = 0; j < tfce->GetNumberOfTuples(); j++)
          {
          double su = tfce->GetTuple1(j);
          double pc = 1.0 - (lower_bound(hTFCE.begin(), hTFCE.end(), su) - hTFCE.begin()) / ((double) p.np);
          tfce_pcorr->SetTuple1(j, pc);
          }
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
      else if(arg == "-R")
        {
        p.row_mask = argv[++i];
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
        p.min_valid_obs = atof(argv[++i]);
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
      else if(arg == "--tfce")
        {
        p.tfce_delta_h = atof(argv[++i]);
        }
      else if(arg == "--tfce-h")
        {
        p.tfce_H = atof(argv[++i]);
        }
      else if(arg == "--tfce-e")
        {
        p.tfce_E = atof(argv[++i]);
        }
      else if(arg == "-z")
        {
        p.flag_z_transform = true;
        }
      else if(arg == "--threads")
        {
        p.max_threads = atof(argv[++i]);
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
