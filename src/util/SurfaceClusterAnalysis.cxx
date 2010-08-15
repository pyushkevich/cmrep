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
#include "vtkTriangle.h"
#include "vnl/vnl_file_matrix.h"
#include "vnl/vnl_rank.h"
#include "vnl/algo/vnl_matrix_inverse.h"


#include "Registry.h"

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int usage()
{
  cout << "This program performs cluster analysis on a VTK mesh (both surface and volume meshes)" << endl;
  cout << "Usage: " << endl;
  cout << "    meshcluster configfile.txt DataArray Threshold [Suffix] [TypeOfTest] [design_mat.txt] [contrast_vec.txt]" << endl;
  cout << "    TypeOfTest = 0 or none means two sample t-test" << endl;
  cout << "    TypeOfTest = 1 means correlation t-test" << endl;
  cout << "    TypeOfTest = 2 means paired t-test" << endl;
  cout << "    TypeOfTest = 3 means correlation with a summary measure" << endl;
  cout << "    TypeOfTest = 4 means GLM: design_mat.txt and contrast_vec.txt are required" << endl;
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

/*
class ClusterGenerator
{
public:
  // Cluster definition
  typedef vector<Cluster> ClusterArray;

  // Constructor
  ClusterGenerator(vtkPolyData *mesh);

  // Compute clusters
  ClusterArray ComputeClusters(const char *data, const char *aelt, double thresh);

  // Get the output mesh
  vtkPolyData *GetOutputMesh()
    { return fConnect->GetOutput(); }

private:
  // Pipeline elements
  vtkClipPolyData *fContour;
  vtkThreshold *fThresh;
  vtkPolyDataConnectivityFilter *fConnect;

  // Meshes
  vtkPolyData *mesh;
};

ClusterGenerator::ClusterGenerator(vtkPolyData *inmesh)
{
  // Store the mesh
  mesh = inmesh;

  // Generate the pipeline
  fContour = vtkClipPolyData::New();
  fContour->SetInput(mesh);

  // fThresh = vtkThreshold::New();
  // fThresh->SetInput(fContour->GetOutput());
  // fThresh->ThresholdByLower(0);
  // fThresh->SetAttributeModeToUsePointData();

  static double srange[] = {-0.5, 0.5};
  
  fConnect = vtkPolyDataConnectivityFilter::New();
  fConnect->SetInput(fContour->GetOutput());
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  // fConnect->ScalarConnectivityOn();
  fConnect->SetScalarRange(srange);
}

ClusterGenerator::ClusterArray
ClusterGenerator::ComputeClusters(const char *data, const char *aelt, double thresh)
{
  // Update the mesh
  // mesh->GetPointData()->SetActiveScalars(data);
  // mesh->GetPointData()->CopyAllOn();

  // Compute the clusters
  fContour->SetValue(thresh);
  fConnect->Update();
  vtkPolyData *f = fContour->GetOutput();
  vtkPolyData *p = fConnect->GetOutput();

  vtkDataArray *daRegion = p->GetPointData()->GetScalars();
  vtkDataArray *daData = f->GetPointData()->GetArray(data);
  vtkDataArray *daArea = f->GetPointData()->GetArray(aelt);

  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
  for(size_t i = 0; i < p->GetNumberOfPoints(); i++)
    {
    size_t region = (size_t) (daRegion->GetTuple1(i));
    double x = daData->GetTuple1(i);
    double a = daArea->GetTuple1(i);
    ca[region].n++;
    ca[region].area += a;
    ca[region].power += a * x;
    }

  
  //for(size_t c = 0; c < ca.size(); c++)
  //  {
  //  printf("Cluster %d: n = %d, area = %f, power = %f, mean_t = %f\n",
  //    c, ca[c].n, ca[c].area, ca[c].power, ca[c].power / ca[c].area);
  //  }
    


  // Return the cluster array
  return ca;
}*/

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


template <class TMeshType>
ClusterArray ComputeClusters(
  TMeshType *mesh,
  const char *data, 
  const char *domain,
  double thresh, 
  TMeshType **mout = NULL)
{ ClusterArray ca(1) ; return ca; }


template <>
ClusterArray ComputeClusters(
  vtkPolyData *mesh,
  const char *data, 
  const char *domain,
  double thresh, 
  vtkPolyData **mout )
{
  vtkClipPolyData *fContour;
  vtkThresholdPoints *fThresh;
//    cout << "creating " << fThresh ; 
  vtkPolyData *f; 
  int NumberOfPoints;
  if (!strcmp(domain, "Point"))
  {
     // Initialize mesh
     mesh->GetPointData()->SetActiveScalars(data);
     fContour = vtkClipPolyData::New();

     // Clip the data field at the threshold
     fContour->SetInput(mesh);
     fContour->SetValue(thresh);
     fContour->Update();
     f = fContour->GetOutput();
  } 
  else
  {  
     cerr << "*********polydata clustering based on cell attribute arrays not supported yet**********" << endl;
     throw("polydata with cells not supported");
     // Initialize mesh
     mesh->GetCellData()->SetActiveScalars(data);
     fThresh = vtkThresholdPoints::New();
    
     // Threshold the cell data field
     //fThresh->SetAttributeModeToUseCellData();
     fThresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, data);
     fThresh->SetInput(mesh);
     fThresh->ThresholdByUpper(thresh);
     fThresh->Update();
     f = fThresh->GetOutput();
  } 


  // Get the connected components
  vtkPolyDataConnectivityFilter * fConnect = vtkPolyDataConnectivityFilter::New();
  if (!strcmp(domain, "Point"))
     fConnect->SetInput(fContour->GetOutput());
  else
     fConnect->SetInput(fThresh->GetOutput());
  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  fConnect->Update();
  vtkPolyData *p = fConnect->GetOutput();
//WriteMesh<vtkPolyData>(mesh, "input.vtk");
//WriteMesh<vtkPolyData>(f, "output.vtk");

  if (!strcmp(domain, "Point"))
     NumberOfPoints = p->GetNumberOfPoints();
  else
     NumberOfPoints = p->GetNumberOfCells();


  // Create output data arrays for computing area element
  vtkFloatArray *daArea = vtkFloatArray::New();
  daArea->SetName("area_element");
  daArea->SetNumberOfComponents(1);
  daArea->SetNumberOfTuples(NumberOfPoints);
  daArea->FillComponent(0, 0.0);

  // Compute the area of each triangle in the cluster set
  for(int k = 0; k < p->GetNumberOfCells(); k++)
    {
    vtkCell *cell = p->GetCell(k);
    if(cell->GetCellType() != VTK_TRIANGLE)
{cout << "Warning: cell type " << cell->GetCellType() ;
      throw("Wrong cell type");
}
    // Compute the area of the triangle
    vtkIdType a0 = cell->GetPointId(0);
    vtkIdType a1 = cell->GetPointId(1);
    vtkIdType a2 = cell->GetPointId(2);
    double p0[3], p1[3], p2[3];
    p->GetPoint(a0, p0);
    p->GetPoint(a1, p1);
    p->GetPoint(a2, p2);

    double area = vtkTriangle::TriangleArea(p0, p1, p2);

    if (!strcmp(domain, "Point"))
    {  
       // Split the volume between neighbors
       daArea->SetTuple1(a0, area / 3.0 + daArea->GetTuple1(a0));
       daArea->SetTuple1(a1, area / 3.0 + daArea->GetTuple1(a1));
       daArea->SetTuple1(a2, area / 3.0 + daArea->GetTuple1(a2));
    }
    else // No need to split, working with cells
       daArea->SetTuple1(k, area);

    }

  // The the important arrays in the resulting meshes
  vtkDataArray *daRegion;
  vtkDataArray *daData;
  if (!strcmp(domain, "Point"))
  {  
     // The important arrays in the resulting meshes
     daRegion = p->GetPointData()->GetScalars();
     daData = f->GetPointData()->GetArray(data);
  }
  else
  {
     daRegion = p->GetCellData()->GetScalars();
     daData = f->GetCellData()->GetArray(data);
  }


  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
  for(int i = 0; i < NumberOfPoints; i++)
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
    if (!strcmp(domain, "Point"))
       p->GetPointData()->AddArray(daArea);
    else
       p->GetCellData()->AddArray(daArea);

    *mout = p;
    }
  else
    {
    // Delete the intermediates
    daArea->Delete();
    fConnect->Delete();
    if (!strcmp(domain, "Point"))
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
  const char *domain, 
  double thresh, 
  vtkUnstructuredGrid **mout )
{

/*** Segment 1 old code begin 
  // Initialize mesh
  mesh->GetPointData()->SetActiveScalars(data);

  // Clip the data field at the threshold
  vtkClipDataSet *fContour = vtkClipDataSet::New();
  fContour->SetInput(mesh);
  fContour->SetValue(thresh);
  vtkUnstructuredGrid *f = fContour->GetOutput();
 Segment 1 old code end ****/
  
/*** Segment 1 new code begin ****/
  vtkClipDataSet *fContour;
  vtkThreshold *fThresh;
//    cout << "creating " << fThresh ;
  vtkUnstructuredGrid *f;
  int NumberOfPoints;
  if (!strcmp(domain, "Point"))
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
     //fThresh->SetAttributeModeToUseCellData();
     fThresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS); 
     fThresh->SetInput(mesh);
     fThresh->ThresholdByUpper(thresh);
     fThresh->Update();
     f = fThresh->GetOutput();
  }
/*** Segment 1 new code end ****/


  // Get the connected components
  vtkConnectivityFilter * fConnect = vtkConnectivityFilter::New();
/*** Segment 2 old code begin 
  fConnect->SetInput(fContour->GetOutput());
 Segment 2 old code end ****/
/*** Segment 2 new code begin ****/
  if (!strcmp(domain, "Point"))
     fConnect->SetInput(fContour->GetOutput());
  else
     fConnect->SetInput(fThresh->GetOutput());
/*** Segment 2 new code end ****/

  fConnect->SetExtractionModeToAllRegions();
  fConnect->ColorRegionsOn();
  fConnect->Update();
  vtkUnstructuredGrid *p = fConnect->GetOutput();

/*** Segment 3 new code begin ****/
  if (!strcmp(domain, "Point"))
     NumberOfPoints = p->GetNumberOfPoints();
  else
     NumberOfPoints = p->GetNumberOfCells();
/*** Segment 3 new code end ****/


  // Create output data arrays for computing volume element
  vtkFloatArray *daArea = vtkFloatArray::New();
  daArea->SetName("volume_element");
  daArea->SetNumberOfComponents(1);
// Segment 4 old code
//  daArea->SetNumberOfTuples(p->GetNumberOfPoints());
// Segment 4 new code
  daArea->SetNumberOfTuples(NumberOfPoints);
  daArea->FillComponent(0, 0.0);

  // Compute the volume of each tetra in the cluster set
  for(size_t k = 0; k < p->GetNumberOfCells(); k++)
    {
    vtkCell *cell = p->GetCell(k);
/* Some cells are converted to VTK_WEDGE -- TODO
    if(cell->GetCellType() != VTK_TETRA)
      throw("Wrong cell type");
*/
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

/**** Segment 5 old code begin
    // Split the volume between neighbors
    daArea->SetTuple1(a0, abs(area) / 4.0 + daArea->GetTuple1(a0));
    daArea->SetTuple1(a1, abs(area) / 4.0 + daArea->GetTuple1(a1));
    daArea->SetTuple1(a2, abs(area) / 4.0 + daArea->GetTuple1(a2));
    daArea->SetTuple1(a3, abs(area) / 4.0 + daArea->GetTuple1(a3));
 Segment 5 old code end ****/
/**** Segment 5 new code begin ***/
    if (!strcmp(domain, "Point"))
    {  
       // Split the volume between neighbors
       daArea->SetTuple1(a0, abs(area) / 4.0 + daArea->GetTuple1(a0));
       daArea->SetTuple1(a1, abs(area) / 4.0 + daArea->GetTuple1(a1));
       daArea->SetTuple1(a2, abs(area) / 4.0 + daArea->GetTuple1(a2));
       daArea->SetTuple1(a3, abs(area) / 4.0 + daArea->GetTuple1(a3));
    }
    else // No need to split, working with cells
       daArea->SetTuple1(k, abs(area));

/**** Segment 5 new code end ***/
    }

/**** Segment 6 old code begin 
  // The the important arrays in the resulting meshes
  vtkDataArray *daRegion = p->GetPointData()->GetScalars();
  vtkDataArray *daData = f->GetPointData()->GetArray(data);
 Segment 6 old code end ***/
/**** Segment 6 new code begin ***/
  vtkDataArray *daRegion;
  vtkDataArray *daData;
  if (!strcmp(domain, "Point"))
  {  
     // The important arrays in the resulting meshes
     daRegion = p->GetPointData()->GetScalars();
     daData = f->GetPointData()->GetArray(data);
  }
  else
  {  
     daRegion = p->GetCellData()->GetScalars();
     daData = f->GetCellData()->GetArray(data);
  }
/**** Segment 6 new code end ***/


  // Build up the cluster array
  ClusterArray ca(fConnect->GetNumberOfExtractedRegions());
// Segment 7 old code
//  for(size_t i = 0; i < p->GetNumberOfPoints(); i++)
// Segment 7 new code
  for(size_t i = 0; i < NumberOfPoints; i++)
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
/**** Segment 8 old code begin 
    p->GetPointData()->AddArray(daArea);
 Segment 8 old code end ***/
/**** Segment 8 new code begin ***/
    if (!strcmp(domain, "Point"))
       p->GetPointData()->AddArray(daArea);
    else
       p->GetCellData()->AddArray(daArea);
/**** Segment 8 new code end ***/

    *mout = p;
    }
  else
    {
    // Delete the intermediates
    daArea->Delete();
    fConnect->Delete();
/**** Segment 9 old code begin 
    fContour->Delete();
 Segment 9 old code end ****/
/**** Segment 9 new code begin ****/
    if (!strcmp(domain, "Point"))
       fContour->Delete();
    else
       fThresh->Delete();
/**** Segment 9 new code end ****/ 

    }

  // Return the cluster array
  return ca;
}


template <class TMeshType>
void ComputeTTest(
  TMeshType *pd, 
  const char *var, 
  const char *domain,
  const vector<int> &labels,
  int l1, int l2, std::string VOIttest)
{
  // Get the pointdata
  vtkDataArray *data;
  vtkDataArray *ttest;
  // Get the pointdata
  int NumberOfPoints;
  if (!strcmp(domain, "Point"))
  {  
     data = pd->GetPointData()->GetArray(var);
     ttest = pd->GetPointData()->GetArray(VOIttest.c_str());
     NumberOfPoints = pd->GetNumberOfPoints();
  }
  else // Get the cell data
  {  
     data = pd->GetCellData()->GetArray(var);
     ttest = pd->GetCellData()->GetArray(VOIttest.c_str());
     NumberOfPoints = pd->GetNumberOfCells();
  }


  // Get the sizes of the cohorts
  int n1 = 0, n2 = 0;
  for(size_t j = 0; j < labels.size(); j++)
    {
    if(labels[j] == l1) 
      n1++;
    else    
      if(labels[j] == l2) 
        n2++;
    }

  // Loop over all the points
  for(size_t i = 0; i < NumberOfPoints; i++)
    {
    // Compute class-wise sums and sums of squares
    double s1 = 0, s2 = 0, ss1 = 0, ss2 = 0;
    for(size_t j = 0; j < labels.size(); j++)
      {
      double x = data->GetComponent(i, j);
      if(labels[j] == l1)
        { s1 += x; ss1 += x * x; }
      else if(labels[j] == l2)
        { s2 += x; ss2 += x * x; } 
      }

    // Compute t2 directly from sums and sums of squares
    double m1 = s1 / n1, m2 = s2 / n2;
    double v1 = ss1 - s1 * m1, v2 = ss2 - s2 * m2;
    double den = sqrt((v1+v2) * (1./n1 + 1./n2) / (n1 + n2 - 2.));
    double t = den > 0 ? (m1 - m2) / den : 0;

    // Add the t2 value to the array
    ttest->SetTuple1(i, t);
    }
}

template <class TMeshType>
void ComputeCorrTest(
  TMeshType *pd, 
  const char *var, 
  const char *domain,
  const vector<int> &labels,
  const vector<int> &indiv_labels,
  int l1, int l2, int ispaired, const vector<float> &corrVar, std::string VOIttest)
{
  // Get the pointdata
  vtkDataArray *data;
  vtkDataArray *ttest;
  // Get the pointdata
  int NumberOfPoints;
  if (!strcmp(domain, "Point"))
  {  
     data = pd->GetPointData()->GetArray(var);
     ttest = pd->GetPointData()->GetArray(VOIttest.c_str());
     NumberOfPoints = pd->GetNumberOfPoints();
  }
  else // Get the cell data
  {  
     data = pd->GetCellData()->GetArray(var);
     ttest = pd->GetCellData()->GetArray(VOIttest.c_str());
     NumberOfPoints = pd->GetNumberOfCells();
  }


  // Create an r (correlation coefficient) array
  vtkFloatArray *array = vtkFloatArray::New();
  vtkDataArray *corrcoef ;
  if (ispaired == 1 || ispaired == 3)
  {
  	array->SetName("R");
  	array->SetNumberOfComponents(1);
  	array->SetNumberOfTuples(NumberOfPoints);
        if (!strcmp(domain, "Point"))
        {  
           pd->GetPointData()->AddArray(array);
           corrcoef = pd->GetPointData()->GetArray("R");
        }
        else
        {  
           pd->GetCellData()->AddArray(array);
           corrcoef = pd->GetCellData()->GetArray("R");
        }

  }

  int n;
  if (ispaired == 1 || ispaired == 2)
     n = indiv_labels.size()/2;
  else if (ispaired == 3)
     n = indiv_labels.size();
  // Loop over all the points
  for(size_t i = 0; i < NumberOfPoints; i++)
    {
    // Compute class-wise sums and sums of squares
    double s1 = 0, s2 = 0, ss1 = 0, ss2 = 0, s1s2 = 0, s12 = 0, ss12 = 0;
    for(size_t j = 0; j < n; j++)
      {
      double x = data->GetComponent(i, j);
      double y;
      if (ispaired == 1 || ispaired == 2)
         y = data->GetComponent(i, indiv_labels[j] + n);
      else if (ispaired == 3)
         y = corrVar[ indiv_labels[j] ];
        { s1 += x; ss1 += x * x; }
        { s2 += y; ss2 += y * y; } 
        s1s2 += x * y;
        if(labels[j] == l1)
	  { s12 += x - y; ss12 += (x - y)*(x - y); }
        else if(labels[j] == l2)
	  { s12 += y - x; ss12 += (x - y)*(x - y); }
	
      }

    double r = 0,t = 0;
    // Compute t2 directly from sums and sums of squares
    if (ispaired == 2) // paired t-test
    {
	double numerator = sqrt(n )*s12;
	double denominator = sqrt(n*ss12 - s12*s12);
	t = numerator/denominator;
    }
    else // paired correlation
    {
    	double numerator = n*s1s2 - s1 * s2;
    	double denominator = sqrt(n*ss1 - s1 * s1) * sqrt(n*ss2 - s2 * s2);
    	r = numerator/denominator;
    	t = r * sqrt((n-2) / (1 - r*r)); 
    	// If testing for negative correlation flip sign of t
    	if (l1 == 1) t = -t;
    }

    // Add the t2 value to the array
    ttest->SetTuple1(i, t);
    if (ispaired == 1 || ispaired == 3)
    	corrcoef->SetTuple1(i, r);
    }
}

template <class TMeshType>
void GeneralLinearModel(const vnl_matrix<double> &origmat, const vnl_matrix<double> &origcon, const vector<int> &indiv_labels, TMeshType *pd, const char *var, std::string VOIttest, const char *domain)
{

  vtkDataArray *data;
  vtkDataArray *ttest;
  // Get the pointdata
  int NumberOfPoints;
  if (!strcmp(domain, "Point"))
  {  
     data = pd->GetPointData()->GetArray(var);
     ttest = pd->GetPointData()->GetArray(VOIttest.c_str());
     NumberOfPoints = pd->GetNumberOfPoints();
  }
  else // Get the cell data
  {  
     data = pd->GetCellData()->GetArray(var);
     ttest = pd->GetCellData()->GetArray(VOIttest.c_str());
     NumberOfPoints = pd->GetNumberOfCells();
  }

  vnl_matrix<double> mat;
  mat.set_size(origmat.rows(), origmat.cols());
  // Shuffle the rows according to current permutation
  for (size_t i = 0; i < indiv_labels.size(); i++)
      for (size_t j = 0; j < origmat.cols(); j++)
          mat[i][j] = origmat[indiv_labels[i]][j];

  // copy contrast matrix
  vnl_matrix<double> con;
  con.set_size(origcon.rows(), origcon.cols());
  for (size_t i = 0; i < origcon.rows(); i++)
      for (size_t j = 0; j < origcon.cols(); j++)
          con[i][j] = origcon[i][j];



  // Load all images into a Y matrix (can we do this)
  size_t n = pd->GetNumberOfPoints();
  vnl_matrix<double> Y(mat.rows(), n);
  for(size_t i = 0; i < mat.rows(); i++)
    {
    for(size_t j = 0; j < n; j++)
      { Y(i,j) = data->GetComponent(j,i); }
    }

  //cout << "Running GLM on " << mat.rows() << " images" << endl;
  //cout << "  design matrix: "  << endl << mat << endl;
  //cout << "  contrast vector: " << con << endl;

  // Compute degrees of freedom
  unsigned int rank = vnl_rank( mat, vnl_rank_row);
  unsigned int df = Y.rows() - rank;

  // Some matrices
  vnl_matrix<double> A = vnl_matrix_inverse<double>(mat.transpose() * mat).pinverse(rank);
  // Compute bhat 
  vnl_matrix<double> bhat = (A * mat.transpose()) * Y;

  //cout << "df " << df << " A " << A << endl;

  // Compute the contrast
  vnl_matrix<double> res = con * bhat;

  // error
  vnl_matrix<double> errmat = Y - mat * bhat;


  // Residual variance
  vnl_matrix<double> resvar(1, n);
  resvar.fill(0.0);
  for(size_t j = 0; j < n; j++)
     {
     for(size_t i = 0; i < mat.rows(); i++)
        resvar(0, j) += errmat(i,j)*errmat(i,j);
     resvar(0,j) = resvar(0,j)/(double)df;
     }

  // t-value
  vnl_matrix<double> tmat(1, n);
  vnl_matrix<double> den(1, n);

  den = (con * (A * con.transpose())) * resvar;
  //den = den.apply(sqrt);
  tmat.fill(0.0);
  for (size_t j = 0; j < n; j++)
      if ( den(0,j) > 0 )
         { 
         tmat(0,j) = res(0,j) / sqrt(den(0,j)); 
    //     cout << tmat(0,j) << " " <<  res(0,j) << " " <<  den(0,j) << endl;
         } 

  // Write the output 
  for(size_t j = 0; j < n; j++)
    { ttest->SetTuple1(j, tmat(0,j)); }
}


template <class TMeshType>
int meshcluster(int argc, char *argv[], Registry registry, bool isPolyData)
{

  // design matrix and contrast vector files
  string fn_design("design_mat.txt");
  string fn_con("contrast_vec.txt");  
  vnl_matrix<double> mat;
  vnl_matrix<double> con;

  // Get the number of permutations
  size_t np = registry["Analysis.NumberOfPermutations"][1000];

  // Get the cohort labels
  vector<int> labels = registry.Folder("Cohort").GetArray(-1);
  vector<int> true_labels = labels;

  // Get the variables of interest
  string sVOI = registry["Analysis.TestVariable"][""];
  if(argc > 2) sVOI = argv[2];
  cout << "Variable of Interest: " << sVOI << endl;

  // Get the clustering domain
  string domain = registry["Analysis.ClusterDomain"][""];
  cout << "Cluster Domain: " << domain << endl;

  // Get the clustering threshold
  double thresh = registry["Analysis.ClusteringThreshold"][0.0];
  if(argc > 3) thresh = atof(argv[3]);
  cout << "Threshold: " << thresh << endl;

  // Paired test desired. TODO registry ?
  int  ispaired = 0;
  if (argc > 5) ispaired = atoi(argv[5]);


  // SD If a correlation is desired, assume that the true labels are in the form of 000.. 111..
  // and that in each group the corresponding positions are paired, as in 123.. 123..
  // In this case, the permutation needs to be done only on one of the group's individual labels 123.. 
  // first create the individual labels with the above assumption
  
  
  int groupSize, Nlabels;
  vector<int> indiv_labels, true_indiv_labels ;
  vector<float> corrVar;
  if (ispaired == 1 ) // correlation between two variables 
  {
     Nlabels = (int)labels.size();
     if (Nlabels % 2 !=0)
        { cerr << "Total array size is odd, must be even for paired tests" << endl; return -1; }
     else
        groupSize = Nlabels/2;
     cout << "Generating individual labels for correlation, group size is " << groupSize << endl;
     for (int cohort = 0; cohort < 2; cohort++)
         for (int i=0; i< groupSize; i++) 
             indiv_labels.push_back( i ); 
     true_indiv_labels = indiv_labels;
  }
  else if (ispaired == 2) // paired t-test
  {
     Nlabels = (int)labels.size();
     if (Nlabels % 2 !=0)
        { cerr << "Total array size is odd, must be even for paired tests" << endl; return -1; }
     else
        groupSize = Nlabels/2;
     cout << "Generating individual labels for paired test, group size is " << groupSize << endl;
     for (int cohort = 0; cohort < 2; cohort++)
         for (int i=0; i< groupSize; i++)
             indiv_labels.push_back( i );
     true_indiv_labels = indiv_labels;
  }
  else if (ispaired == 3) // correlation with summary measurement
  {
     FILE *fd;
     fd = fopen("variable.txt","r");
     Nlabels = (int)labels.size();
     groupSize = Nlabels;
     float val;
     cout << "Variable values: " ;
     for (int i=0; i< groupSize; i++) 
     {
         indiv_labels.push_back( i );
         fscanf(fd, "%f\n", &val);
         corrVar.push_back( val ); 
         cout << val << " " ; 
     }
     true_indiv_labels = indiv_labels;
     cout << endl;
  }
  else if (ispaired == 4) // GLM
  {
     // design matrix and contrast vector
     if (argc == 8 ) 
     { 
       fn_design = string(argv[6]); 
       fn_con = string(argv[7]);
     } 
     // Read the matrix from file
     try{
     mat = vnl_file_matrix<double>(fn_design.c_str());
     }
     catch (char * str)
       { cout << "Unable to read matrix from file given" << endl; return -1;}

     // Read the contrast from file
     try{
     con = vnl_file_matrix<double>(fn_con.c_str());
     }
     catch (char * str)
       { cout << "Unable to read contrast from file given" << endl; return -1;}

     // Check that the number of images matches
     if(labels.size() != mat.rows())
       { throw string("Matrix number of rows does not match number of observations"); return -1;}

     // Check that the columns in matrix match contrast vector
     if(con.columns() != mat.columns())
       { throw string("Matrix and contrast vector must have same number of columns"); return -1;}

     Nlabels = (int)labels.size();
     groupSize = Nlabels;
     cout << "GLM: group size is " << groupSize << endl;
     for (int i=0; i< groupSize; i++) 
         indiv_labels.push_back( i );
     true_indiv_labels = indiv_labels;
  }

  // If the threshold is negative, simply change the direction of the test
  string posneg("_pos");
  int l1 = 0, l2 = 1;
  if(thresh < 0)
    { l1 = 1; l2 = 0; thresh = -thresh; posneg = "_neg";}
  string VOIttest = sVOI + "ttest" + posneg;

  // Load each of the meshes and create a cluster analyzer
  vector<string> fnMeshes, fnOutMeshes;
  vector<ClusterArray> claTrue;
  fnMeshes = registry.Folder("Mesh").GetArray(string(""));
  fnOutMeshes = registry.Folder("OutputMesh").GetArray(string(""));
  if(fnMeshes.size() == 0)
    { cerr << "Missing mesh specification" << endl; return -1; }


  // If there is a suffix, add it to the output meshes
  if(argc > 4) 
    {
    string suffix = argv[4];
    string teststring;
    std::stringstream ss;
    ss << thresh;
    if (ispaired ==0)
    	teststring = "ttst";
    else if (ispaired == 1)
        teststring = "corr";
    else if (ispaired == 2)
    	teststring = "ptst";
    else if (ispaired == 3)
    	teststring = "vcorr";
    else if (ispaired == 4)
    	teststring = "glm";
    else
	teststring = "unkn"; 
    for(size_t i = 0; i < fnOutMeshes.size(); i++)
      fnOutMeshes[i] = fnOutMeshes[i] + teststring + ss.str() + posneg + suffix;
    }
  
  // Read the meshes
  vector<TMeshType *> mesh;
  for(size_t i = 0; i < fnMeshes.size(); i++)
    {
    // Read mesh
    cout << "Reading mesh " << fnMeshes[i] << endl;
    mesh.push_back(ReadMesh<TMeshType>(fnMeshes[i].c_str()));

    // Create a t-test array
    vtkFloatArray *array = vtkFloatArray::New();
    array->SetName(VOIttest.c_str());
    array->SetNumberOfComponents(1);
    // cell or point based analysis ?
    if (!strcmp(domain.c_str(), "Point"))
    {
       array->SetNumberOfTuples(mesh[i]->GetNumberOfPoints());
       mesh[i]->GetPointData()->AddArray(array);
       mesh[i]->GetPointData()->SetActiveScalars(VOIttest.c_str());
    }
    else if (!strcmp(domain.c_str(), "Cell"))
    {
       array->SetNumberOfTuples(mesh[i]->GetNumberOfCells());
       mesh[i]->GetCellData()->AddArray(array);
       mesh[i]->GetCellData()->SetActiveScalars(VOIttest.c_str());
    }
    else
    { cerr << "Unknown clustering domain: must be Point or Cell" << endl; return -1; }
    }
    
  // Run permutation analysis
  vector<double> hArea(np), hPower(np);
  for(size_t ip = 0; ip < np; ip++)
    {
    // Apply a random permutation to the labels array
    random_shuffle(labels.begin(), labels.end());
    hArea[ip] = 0; hPower[ip] = 0;

    // SD shuffle individual member labels within group for correlations
    if (ispaired == 1 || ispaired == 3 || ispaired == 4)
    {
       vector<int>::iterator it;
       random_shuffle(indiv_labels.begin(), indiv_labels.begin()+groupSize );
       /*
       for (it=indiv_labels.begin(); it!=indiv_labels.end(); ++it)
           cout << *it << " " << flush ;
       cout << endl;
       return 0;
       */
    }    
    // Build up the histogram of cluster areas (and powers)
    for(size_t i = 0; i < fnMeshes.size(); i++)
      {
      // SD paired tests
      if (ispaired > 0 && ispaired < 4)
         ComputeCorrTest<TMeshType>(mesh[i], sVOI.c_str(), domain.c_str(), labels, indiv_labels, l1, l2, ispaired, corrVar, VOIttest);
      else if (ispaired == 4)
         {
      // GLM
         GeneralLinearModel<TMeshType>(mat, con, indiv_labels, mesh[i], sVOI.c_str(), VOIttest, domain.c_str());
         }
      // For each mesh, compute the t-test and the clusters
      else
         ComputeTTest<TMeshType>(mesh[i], sVOI.c_str(), domain.c_str(), labels, l1, l2, VOIttest);
      ClusterArray ca = ComputeClusters<TMeshType>(mesh[i], VOIttest.c_str(), domain.c_str(), thresh);
    //WriteMesh<TMeshType>(mesh[i], fnMeshes[i].c_str());
    //exit(0);

      // Now find the largest cluster
      for(size_t c = 0; c < ca.size(); c++)
        {
        if(ca[c].area > hArea[ip]) hArea[ip] = ca[c].area;
        if(ca[c].power > hPower[ip]) hPower[ip] = ca[c].power;
        }
      }
    cout << "." << flush;
    if((ip+1) % 80 == 0 || (ip+1) == np) 
      cout << " " << ip+1 << endl;
    }

  // Sort the histograms
  sort(hArea.begin(), hArea.end());
  sort(hPower.begin(), hPower.end());

  // Going back to the original meshes, assign a cluster p-value to each mesh
  for(size_t i = 0; i < fnMeshes.size(); i++)
    {
    TMeshType *mout;

    // SD paired tests
    if (ispaired > 0 && ispaired < 4)
       ComputeCorrTest<TMeshType>(mesh[i], sVOI.c_str(), domain.c_str(), true_labels, true_indiv_labels, l1, l2, ispaired, corrVar, VOIttest);
    else if (ispaired == 4)
    // GLM
       GeneralLinearModel<TMeshType>(mat, con, true_indiv_labels, mesh[i], sVOI.c_str(), VOIttest, domain.c_str());
    else
    // Compute the t-test for the mesh with the correct labels
       ComputeTTest<TMeshType>(mesh[i], sVOI.c_str(), domain.c_str(), true_labels, l1, l2, VOIttest);

    // Compute the clusters for the mesh with the correct labels
    ClusterArray ca = ComputeClusters<TMeshType>(mesh[i], VOIttest.c_str(), domain.c_str(), thresh, &mout);
    printf("MESH %s HAS %d CLUSTERS: \n", fnMeshes[i].c_str(), ca.size());

    // Assign a p-value to each cluster
    for(size_t c = 0; c < ca.size(); c++)
      {
      // Brute force search in the histogram :(
      size_t zArea = 0, zPower = 0;
      while(zArea < np && hArea[zArea] < ca[c].area) zArea++;
      while(zPower < np && hPower[zPower] < ca[c].power) zPower++;
      ca[c].pArea = 1.0 - zArea * 1.0 / np;
      ca[c].pPower = 1.0 - zPower * 1.0 / np;
      bool sig = (ca[c].pArea <= 0.05 || ca[c].pPower <= 0.05);
      printf("Cluster %03d:  AvgT = %6f; Area = %6f (p = %6f);  Power = %6f (p = %6f); %s\n",
        c, ca[c].tvalue/ca[c].n, ca[c].area, ca[c].pArea, ca[c].power, ca[c].pPower, 
        sig ? "***" : "");
      }

    // Create output mesh arrays for p-values
    string snArea = sVOI + string(" p-cluster-area");
    vtkFloatArray *aArea = vtkFloatArray::New();
    aArea->SetName(snArea.c_str());
    aArea->SetNumberOfComponents(1);
    aArea->SetNumberOfTuples(mout->GetNumberOfPoints());
    mout->GetPointData()->AddArray(aArea);

    string snPower = sVOI + string(" p-cluster-power");
    vtkFloatArray *aPower = vtkFloatArray::New();
    aPower->SetName(snPower.c_str());
    aPower->SetNumberOfComponents(1);
    aPower->SetNumberOfTuples(mout->GetNumberOfPoints());
    mout->GetPointData()->AddArray(aPower);

    // Set the mesh arrays' p-values
    for(size_t p = 0; p < mout->GetNumberOfPoints(); p++)
      {
      size_t r = (size_t) mout->GetPointData()->GetScalars()->GetTuple1(p);
      aArea->SetTuple1(p, ca[r].pArea);
      aPower->SetTuple1(p, ca[r].pPower);
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
    WriteMesh<TMeshType>(mout, fnOutMeshes[i].c_str());
    cout << "Saving input mesh.." << endl;
    WriteMesh<TMeshType>(mesh[i], fnMeshes[i].c_str());
    }
    return EXIT_SUCCESS;

}

int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();
  // Read the registry file
  Registry registry(argv[1]);

  vector<string> fnMeshes;
  fnMeshes = registry.Folder("Mesh").GetArray(string(""));
  if(fnMeshes.size() == 0)
    { cerr << "Missing mesh specification" << endl; return -1; }


  
  // Read the meshes
  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(fnMeshes[0].c_str());
  reader->OpenVTKFile();
  reader->ReadHeader();

  bool isPolyData = true;
  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    isPolyData = false;
    return meshcluster<vtkUnstructuredGrid>( argc, argv, registry, isPolyData);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return meshcluster<vtkPolyData>( argc, argv, registry, isPolyData);

    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}
