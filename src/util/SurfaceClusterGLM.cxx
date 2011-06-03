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
#include "vtkUnstructuredGridToPolyDataFilter.h"
#include "vtkMergeCells.h"
#include "vtkGeometryFilter.h"
#include "vtkTriangle.h"
#include "vnl/vnl_file_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_rank.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

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

int usage()
{
  cout << 
    "meshglm - a statistical tool for VTK meshes\n"
    "usage:\n"
    "  meshglm [options]\n"
    "options:\n"
    "  -m / --mesh mesh_in.vtk mesh_out.vtk      Input/output mesh (may be repeated)\n"
    "  -a / --array NAME                         Name of array\n"
    "  -p / --perm N                             Number of random permutations\n"
    "  -g / --glm design.txt contrast.txt        GLM specification (like SPM)\n"
    "  -c / --cell                               [optional] : run on cell data\n"
    "  -t / --threshold T                        The threshold to define clusters\n"
    "  -d / --diffusion T                        Amount of diffusion to apply\n"
    "                                            before statistical analysis\n"
    "  -e / --edges                              Generate separate outputs containing edges of\n"
    "                                            clusters.\n"
    "mesh specification:\n"
    "  - The input mesh should contain a k-variate array NAME, where k is the number\n"
    "    of entries in the GLM\n"
    "  - The output mesh will contain various statistical maps derived from the GLM\n"
    "  as well as a cluster membership function (in/out) for each cluster     \n";
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
  Cluster() : n(0), area(0.0), power(0.0), tvalue(0.0) {};
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
void WriteMesh(TMeshType *mesh, const char *fname)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname)
{
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fname);
  writer->SetInput(mesh);
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fname);
  writer->SetInput(mesh);
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
    for(size_t j = 0; j < cell->GetNumberOfPoints(); j++)
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

    vtkFloatArray *arr = vtkFloatArray::New();
    arr->SetName(src1->GetName());
    arr->SetNumberOfComponents(src1->GetNumberOfComponents());
    arr->SetNumberOfTuples(out->GetNumberOfCells());
  
    for(vtkIdType i = 0; i < p1->GetNumberOfCells(); i++)
      for(int m = 0; m < src1->GetNumberOfComponents(); m++)
        arr->SetComponent(i, m, src1->GetComponent(i, m));

    for(vtkIdType i = 0; i < p2->GetNumberOfCells(); i++)
      for(int m = 0; m < src1->GetNumberOfComponents(); m++)
        arr->SetComponent(i + p1->GetNumberOfCells(), m, src2->GetComponent(i, m));

    out->GetCellData()->AddArray(arr);
    }

  // Merge all the point arrays
  for(int k = 0; k < p1->GetPointData()->GetNumberOfArrays(); k++)
    {
    int index;
    vtkDataArray *src1 = p1->GetPointData()->GetArray(k);
    vtkDataArray *src2 = p2->GetPointData()->GetArray(src1->GetName(), index);

    if(index < 0 || src1->GetNumberOfComponents() != src2->GetNumberOfComponents())
      continue;

    vtkFloatArray *arr = vtkFloatArray::New();
    arr->SetName(src1->GetName());
    arr->SetNumberOfComponents(src1->GetNumberOfComponents());
    arr->SetNumberOfTuples(out->GetNumberOfPoints());

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

    out->GetPointData()->AddArray(arr);
    }
}

template <class TMeshType>
ClusterArray ComputeClusters(
  TMeshType *mesh,
  const char *data, 
  Domain dom,
  double thresh, 
  TMeshType **mout = NULL)
{ ClusterArray ca(1) ; return ca; }

template <>
ClusterArray ComputeClusters(
  vtkPolyData *mesh,
  const char *data, 
  Domain dom,
  double thresh, 
  vtkPolyData **mout )
{
  vtkClipPolyData *fContour;
  vtkThreshold *fThresh;
  vtkDataSet *f; 

  if (dom == POINT)
    {
    // Initialize mesh
    mesh->GetPointData()->SetActiveScalars(data);
    fContour = vtkClipPolyData::New();

    // Clip the data field at the threshold
    fContour->SetInput(mesh);
    fContour->SetValue(thresh);
    fContour->Update();
    fContour->GenerateClippedOutputOn();
    f = fContour->GetOutput();
    } 
  else
    {  
    // Initialize mesh
    mesh->GetCellData()->SetActiveScalars(data);
    fThresh = vtkThreshold::New();

    // Threshold the cell data field
    fThresh->SetInput(mesh);
    fThresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
    fThresh->ThresholdByUpper(thresh);
    fThresh->Update();
    f = fThresh->GetOutput();
    } 

  vtkDataSetAttributes *fdata = (dom == POINT) ? 
    (vtkDataSetAttributes *) f->GetPointData() : 
    (vtkDataSetAttributes *) f->GetCellData();

  // Get the connected components
  vtkConnectivityFilter * fConnect = vtkConnectivityFilter::New();
  fConnect->SetInput(f);
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  fConnect->Update();

  vtkGeometryFilter *geo = vtkGeometryFilter::New();
  geo->SetInput(fConnect->GetOutput());
  geo->Update();
  vtkPolyData *p = geo->GetOutput();

  // Get the attribute data from the connected component filter
  vtkDataSetAttributes *pdata = (dom == POINT) ? 
    (vtkDataSetAttributes *) p->GetPointData() : 
    (vtkDataSetAttributes *) p->GetCellData();

  // Number of elements to process
  int nelt = (dom == POINT) ? 
    p->GetNumberOfPoints() : p->GetNumberOfCells();

  // Create output data arrays for computing area element
  vtkFloatArray *daArea = vtkFloatArray::New();
  daArea->SetName("area_element");
  daArea->SetNumberOfComponents(1);
  daArea->SetNumberOfTuples(nelt);
  daArea->FillComponent(0, 0.0);

  // Compute the area of each triangle in the cluster set
  for(int k = 0; k < p->GetNumberOfCells(); k++)
    {
    vtkCell *cell = p->GetCell(k);
    if(cell->GetCellType() != VTK_TRIANGLE)
      throw MCException("Wrong cell type, should be VTK_TRIANGLE");

    // Compute the area of the triangle
    vtkIdType a0 = cell->GetPointId(0);
    vtkIdType a1 = cell->GetPointId(1);
    vtkIdType a2 = cell->GetPointId(2);
    double p0[3], p1[3], p2[3];
    p->GetPoint(a0, p0);
    p->GetPoint(a1, p1);
    p->GetPoint(a2, p2);

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

  // The the important arrays in the resulting meshes
  vtkDataArray *daRegion = pdata->GetScalars();
  vtkDataArray *daData = fdata->GetArray(data);

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
    if(dom == POINT)
      {
      // Get the cut away data
      vtkPolyData *cut = fContour->GetClippedOutput();

      // Assign each cell in cut away piece a ClusterId of -1
      vtkFloatArray* cutrid = vtkFloatArray::New();
      cutrid->SetNumberOfComponents(1);
      cutrid->SetNumberOfTuples(cut->GetNumberOfCells());
      cutrid->SetName("RegionId");
      cutrid->FillComponent(0, -1);
      cut->GetCellData()->AddArray(cutrid);

      // Assign a RegionId to each cell in the thresholded part
      vtkFloatArray* prid = vtkFloatArray::New();
      prid->SetNumberOfComponents(1);
      prid->SetNumberOfTuples(p->GetNumberOfCells());
      prid->SetName("RegionId");

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

      p->GetCellData()->AddArray(prid);

      // Merge the two arrays
      vtkPolyData *merge = vtkPolyData::New();
      MergeMeshPieces(p, cut, merge);

      *mout = merge;
      }
    else
      {
      // Apply an inverted threshold
      vtkThreshold *fThreshInv = vtkThreshold::New();

      // Threshold the cell data field
      fThresh->SetInput(mesh);
      fThresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
      fThresh->ThresholdByLower(thresh - 1.e-6);
      fThresh->Update();
      vtkUnstructuredGrid *cut = fThresh->GetOutput();

      // Each cell in p already has a region id. Just need to add regionids to the cut
      vtkFloatArray* cutrid = vtkFloatArray::New();
      cutrid->SetNumberOfComponents(1);
      cutrid->SetNumberOfTuples(cut->GetNumberOfCells());
      cutrid->SetName("RegionId");
      cutrid->FillComponent(0, -1);
      cut->GetCellData()->AddArray(cutrid);

      // Merge the pieces
      vtkPolyData *merge = vtkPolyData::New();
      MergeMeshPieces(p, cut, merge);

      *mout = merge;
      }
    }
  else
    {
    // Delete the intermediates
    daArea->Delete();
    fConnect->Delete();
    if(dom == POINT)
       fContour->Delete();
    else
       fThresh->Delete();
    }

  // Return the cluster array
  return ca;
}

template <>
ClusterArray ComputeClusters(
  vtkUnstructuredGrid *mesh,
  const char *data,
  Domain dom,
  double thresh, 
  vtkUnstructuredGrid **mout )
{
  vtkClipDataSet *fContour;
  vtkThreshold *fThresh;

  vtkUnstructuredGrid *f;
  if(dom == POINT)
    {  
    // Initialize mesh
    mesh->GetPointData()->SetActiveScalars(data);
    fContour = vtkClipDataSet::New();

    // Clip the data field at the threshold
    fContour->SetInput(mesh);
    fContour->SetValue(thresh);
    fContour->Update();
    f = fContour->GetOutput();
    }
  else
    {  
     // Initialize mesh
     mesh->GetCellData()->SetActiveScalars(data);
     fThresh = vtkThreshold::New();

     // Threshold the cell data field
     fThresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
     fThresh->SetInput(mesh);
     fThresh->ThresholdByUpper(thresh);
     fThresh->Update();
     f = fThresh->GetOutput();
  }

  vtkDataSetAttributes *fdata = (dom == POINT) ? 
    (vtkDataSetAttributes *) f->GetPointData() : 
    (vtkDataSetAttributes *) f->GetCellData();

  // Get the connected components
  vtkConnectivityFilter * fConnect = vtkConnectivityFilter::New();
  fConnect->SetInput(f);
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  fConnect->Update();
  vtkUnstructuredGrid *p = fConnect->GetOutput();

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
  for(size_t k = 0; k < p->GetNumberOfCells(); k++)
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
  vtkDataArray *daData = fdata->GetArray(data);

  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
  for(size_t i = 0; i < nelt; i++)
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

struct Parameters
{
  // Input and output meshes
  vector<string> fn_mesh_input, fn_mesh_output;

  // Number of permutations
  size_t np;
  
  // Array
  string array_name;

  // Domain
  Domain dom;

  // GLM parameters
  string glm_design, glm_contrast;

  // T statistic threshold
  double threshold;

  // Diffusion amount
  double diffusion;

  // Whether tubes are being used
  bool flag_edges;

  Parameters()
    {
    np = 0; dom = POINT; threshold = 0; diffusion = 0; flag_edges = false;
    }
};

template <unsigned int VDim>
class IdHash : public std::vector<vtkIdType>
{
public:
  IdHash(vtkIdType *src) : std::vector<vtkIdType>(VDim, 0) 
    {
    for(int i = 0; i < VDim; i++)
      (*this)[i] = src[i];
    std::sort(this->begin(), this->end());
    }

  IdHash() : std::vector<vtkIdType>(VDim, -1) {}

  bool operator < (const IdHash<VDim> &b)
    {
    for(int k = 0; k < VDim; k++)
      {
      if(*this[k] < b[k]) return true;
      if(*this[k] > b[k]) return false;
      }
    return false;
    }
};

template <class TMeshType, unsigned int VDim>
int AppendOpenEdgesFromMesh(
  TMeshType *src, vtkDataArray *regarr, int region, TMeshType *trg)
{
  // Hashtable
  typedef IdHash<VDim> HashEntry;
  typedef std::map<HashEntry, int> Hash;
  Hash hash;
  int nadded = 0;

  // Traverse each cell
  for(int i = 0; i < src->GetNumberOfCells(); i++)
    {
    if(regarr->GetTuple1(i) == region)
      {
      vtkCell *cell = src->GetCell(i);

      // Get all edges/faces in the cell
      int ne = (VDim == 2) ? cell->GetNumberOfEdges() : cell->GetNumberOfFaces();
      for(size_t k = 0; k < ne; k++)
        {
        vtkCell *edge = (VDim == 2) ? cell->GetEdge(k) : cell->GetFace(k);
        if(edge->GetNumberOfPoints() != VDim)
          throw MCException("Wrong number of points in edge (%d)", edge->GetNumberOfPoints());
        HashEntry he(edge->GetPointIds()->GetPointer(0));
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

      // Get all edges/faces in the cell
      int ne = (VDim == 2) ? cell->GetNumberOfEdges() : cell->GetNumberOfFaces();
      for(size_t k = 0; k < ne; k++)
        {
        vtkCell *edge = (VDim == 2) ? cell->GetEdge(k) : cell->GetFace(k);
        HashEntry he(edge->GetPointIds()->GetPointer(0));

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

template <class TMeshType>
int meshcluster(Parameters &p, bool isPolyData)
{
  // Load design matrix and contrast vector files
  vnl_matrix<double> mat, con;

  try 
    {
    mat = vnl_file_matrix<double>(p.glm_design.c_str());
    }
  catch(...)
    {
    throw MCException("Error reading design matrix from %s", p.glm_design.c_str());
    }

  try
    { 
    con = vnl_file_matrix<double>(p.glm_contrast.c_str());
    }
  catch(...)
    {
    throw MCException("Error reading contrast vector from %s", p.glm_contrast.c_str());
    }

  // Make sure matrices are compatible
  if(con.columns() != mat.columns())
    throw MCException("Matrix and contrast vector must have same number of columns");

  if(con.rows() != 1)
    throw MCException("Contrast vector must have one row");

  // Create the shuffle arrays
  vector<int> true_indiv_labels, indiv_labels;
  for(int i = 0; i < mat.rows(); i++)
    indiv_labels.push_back(i);
  true_indiv_labels = indiv_labels;

  // If the threshold is negative, simply change the direction of the test
  string posneg("pos");
  if(p.threshold < 0)
    { 
    p.threshold = -p.threshold; 
    posneg = "neg";
    }
  string VOIttest = p.array_name + "_tstat_" + posneg;

  // Load each of the meshes and create a cluster analyzer
  vector<ClusterArray> claTrue;

  // Read the meshes
  vector<TMeshType *> mesh;
  for(size_t i = 0; i < p.fn_mesh_input.size(); i++)
    {
    // Read mesh
    cout << "Reading mesh " << p.fn_mesh_input[i] << endl;
    mesh.push_back(ReadMesh<TMeshType>(p.fn_mesh_input[i].c_str()));

    // Get the appropriate data type
    vtkDataSetAttributes *dsa = (p.dom == POINT) ? 
      (vtkDataSetAttributes *) mesh.back()->GetPointData() : 
      (vtkDataSetAttributes *) mesh.back()->GetCellData();

    // Create a t statistic array
    vtkFloatArray *array = vtkFloatArray::New();
    array->SetName(VOIttest.c_str());
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(dsa->GetNumberOfTuples());

    // Add array to the input mesh
    dsa->AddArray(array);
    dsa->SetActiveScalars(array->GetName());
    }
    
  // Run permutation analysis
  vector<double> hArea(p.np), hPower(p.np);
  for(size_t ip = 0; ip < p.np; ip++)
    {
    // Initialize the histogram at zero
    hArea[ip] = 0; hPower[ip] = 0;

    // Apply a random permutation to the labels array
    random_shuffle(indiv_labels.begin(), indiv_labels.end());

    // Build up the histogram of cluster areas (and powers)
    for(size_t i = 0; i < mesh.size(); i++)
      {
      GeneralLinearModel<TMeshType>(mat, con, indiv_labels, mesh[i], p.array_name.c_str(), VOIttest, p.dom);
      ClusterArray ca = ComputeClusters<TMeshType>(mesh[i], VOIttest.c_str(), p.dom, p.threshold);

      // Now find the largest cluster
      for(size_t c = 0; c < ca.size(); c++)
        {
        if(ca[c].area > hArea[ip]) hArea[ip] = ca[c].area;
        if(ca[c].power > hPower[ip]) hPower[ip] = ca[c].power;
        }
      }
    cout << "." << flush;
    if((ip+1) % 80 == 0 || (ip+1) == p.np) 
      cout << " " << ip+1 << endl;
    }

  // Sort the histograms
  sort(hArea.begin(), hArea.end());
  sort(hPower.begin(), hPower.end());

  // Going back to the original meshes, assign a cluster p-value to each mesh
  for(size_t i = 0; i < mesh.size(); i++)
    {
    TMeshType *mout;

    GeneralLinearModel<TMeshType>(mat, con, true_indiv_labels, mesh[i], p.array_name.c_str(), VOIttest, p.dom);
    ClusterArray ca = ComputeClusters<TMeshType>(mesh[i], VOIttest.c_str(), p.dom, p.threshold, &mout);

    printf("MESH %s HAS %d CLUSTERS: \n", p.fn_mesh_input[i].c_str(), (int) ca.size());

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
    string snArea = p.array_name + string(" p-cluster-area");
    vtkFloatArray *aArea = vtkFloatArray::New();
    aArea->SetName(snArea.c_str());
    aArea->SetNumberOfComponents(1);
    aArea->SetNumberOfTuples(mout->GetNumberOfCells());
    mout->GetCellData()->AddArray(aArea);

    string snPower = p.array_name + string(" p-cluster-power");
    vtkFloatArray *aPower = vtkFloatArray::New();
    aPower->SetName(snPower.c_str());
    aPower->SetNumberOfComponents(1);
    aPower->SetNumberOfTuples(mout->GetNumberOfCells());
    mout->GetCellData()->AddArray(aPower);

    // Get the cluster ID array
    vtkDataArray *aClusterId = mout->GetCellData()->GetArray("RegionId");

    // Set the mesh arrays' p-values
    for(size_t ic = 0; ic < mout->GetNumberOfCells(); ic++)
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
      vtkFloatArray *daArea = vtkFloatArray::New();
      daArea->SetNumberOfComponents(1);
      daArea->Allocate(0);
      daArea->SetName(snArea.c_str());
      emesh->GetCellData()->AddArray(daArea);

      vtkFloatArray *daPower = vtkFloatArray::New();
      daPower->SetNumberOfComponents(1);
      daPower->Allocate(0);
      daPower->SetName(snPower.c_str());
      emesh->GetCellData()->AddArray(daPower);

      vtkFloatArray *daClusterId = vtkFloatArray::New();
      daClusterId->SetNumberOfComponents(1);
      daClusterId->Allocate(0);
      daClusterId->SetName("RegionId");
      emesh->GetCellData()->AddArray(daClusterId);

      // Repeat for each cluster
      for(size_t c = 0; c < ca.size(); c++)
        {
        int added = AppendOpenEdgesFromMesh<TMeshType,2> (mout, aClusterId, c, emesh);
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
      WriteMesh<TMeshType>(emesh, fnedge.c_str());
      }

    // Threshold out the non-significant clusters
    /*
    mout->GetPointData()->SetActiveScalars(snPower.c_str());
    vtkClipPolyData *fClip = vtkClipPolyData::New();
    fClip->SetInput(mout);
    fClip->SetValue(0.05);
    fClip->InsideOutOn();
    fClip->Update();
    */

    // Save the output mesh 
    cout << "Saving output mesh.." << endl;
    WriteMesh<TMeshType>(mout, p.fn_mesh_output[i].c_str());
    // cout << "Saving input mesh.." << endl;
    // WriteMesh<TMeshType>(mesh[i], fnMeshes[i].c_str());
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
      else if(arg == "-p" || arg == "--perm")
        {
        p.np = atoi(argv[++i]);
        }
      else if(arg == "-g" || arg == "--glm")
        {
        p.glm_design = argv[++i];
        p.glm_contrast = argv[++i];
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
        p.diffusion = atoi(argv[++i]);
        }
      else if(arg == "-e" || arg == "--edges")
        {
        p.flag_edges = true;
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
    cerr << "Excepetion caught: " << exc.what() << endl; 
    return -1;
    }
}
