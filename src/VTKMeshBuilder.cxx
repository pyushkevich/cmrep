#include "VTKMeshBuilder.h"
#include "MeshTraversal.h"

#include "vtkCellData.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkPolyDataWriter.h"


template <class TDataSet>
VTKMeshBuilder<TDataSet>::VTKMeshBuilder()
{
  pd = TDataSet::New();
}

template <class TDataSet>
VTKMeshBuilder<TDataSet>::VTKMeshBuilder(TDataSet* other)
{
  pd = TDataSet::New();
  pd->DeepCopy(other);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::SetPoints(const vnl_matrix<double> &x)
{
  assert(x.columns() == 3);
  vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
  pts->SetNumberOfPoints(x.rows());
  for(int i = 0; i < x.rows(); i++)
    pts->SetPoint(i, x(i,0), x(i,1), x(i,2));
  pd->SetPoints(pts);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::SetTriangles(const TriangleMesh &mesh)
{
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();
  for(auto T : mesh.triangles)
    {
    cells->InsertNextCell(3);
    for(unsigned int a = 0; a < 3; a++)
      cells->InsertCellPoint(T.vertices[a]);
    }
  pd->SetPolys(cells);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::SetTriangles(const vnl_matrix<unsigned int> &tri)
{
  vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();
  for(unsigned int i = 0; i < tri.rows(); i++)
    {
    cells->InsertNextCell(3);
    for(unsigned int a = 0; a < 3; a++)
      cells->InsertCellPoint(tri(i,a));
    }
  pd->SetPolys(cells);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::SetNormals(const vnl_matrix<double> &x)
{
  assert(x.columns() == 3);
  assert(pd->GetNumberOfPoints() == x.rows());

  vtkSmartPointer<vtkDoubleArray> arr_nrm = vtkDoubleArray::New();
  arr_nrm->SetNumberOfComponents(3);
  arr_nrm->SetNumberOfTuples(x.rows());

  // Update the points
  for(int i = 0; i < x.rows(); i++)
    arr_nrm->SetTuple3(i, x(i,0), x(i,1), x(i,2));

  pd->GetPointData()->SetNormals(arr_nrm);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::AddArrayInternal(const vnl_matrix<double> &x, const char *name, bool cell)
{
  assert(pd->GetNumberOfPoints() == x.rows());

  vtkSmartPointer<vtkDoubleArray> arr= vtkDoubleArray::New();
  arr->SetNumberOfComponents(x.columns());
  arr->SetNumberOfTuples(x.rows());
  arr->SetName(name);

  // Update the points
  for(int i = 0; i < x.rows(); i++)
    for(int a = 0; a < x.columns(); a++)
      arr->SetComponent(i, a, x(i,a));

  if(cell)
    pd->GetCellData()->AddArray(arr);
  else
    pd->GetPointData()->AddArray(arr);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::AddArrayInternal(const vnl_vector<double> &x, const char *name, bool cell)
{
  assert(pd->GetNumberOfPoints() == x.size());

  vtkSmartPointer<vtkDoubleArray> arr= vtkDoubleArray::New();
  arr->SetNumberOfComponents(1);
  arr->SetNumberOfTuples(x.size());
  arr->SetName(name);

  // Update the points
  for(int i = 0; i < x.size(); i++)
    arr->SetTuple1(i, x[i]);

  if(cell)
    pd->GetCellData()->AddArray(arr);
  else
    pd->GetPointData()->AddArray(arr);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::AddArrayInternal(const vnl_matrix<int> &x, const char *name, bool cell)
{
  assert(pd->GetNumberOfPoints() == x.rows());

  vtkSmartPointer<vtkIntArray> arr= vtkIntArray::New();
  arr->SetNumberOfComponents(x.columns());
  arr->SetNumberOfTuples(x.rows());
  arr->SetName(name);

  // Update the points
  for(int i = 0; i < x.rows(); i++)
    for(int a = 0; a < x.columns(); a++)
      arr->SetComponent(i, a, x(i,a));

  if(cell)
    pd->GetCellData()->AddArray(arr);
  else
    pd->GetPointData()->AddArray(arr);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::AddArrayInternal(const vnl_vector<int> &x, const char *name, bool cell)
{
  assert(pd->GetNumberOfPoints() == x.size());

  vtkSmartPointer<vtkIntArray> arr= vtkIntArray::New();
  arr->SetNumberOfComponents(1);
  arr->SetNumberOfTuples(x.size());
  arr->SetName(name);

  // Update the points
  for(int i = 0; i < x.size(); i++)
    arr->SetTuple1(i, x[i]);

  if(cell)
    pd->GetCellData()->AddArray(arr);
  else
    pd->GetPointData()->AddArray(arr);
}

template <class TDataSet>
void VTKMeshBuilder<TDataSet>::Save(const char *fn)
{
  vtkSmartPointer<vtkPolyDataWriter> w = vtkPolyDataWriter::New();
  w->SetInputData(pd);
  w->SetFileName(fn);
  w->Update();
}

template class VTKMeshBuilder<vtkPolyData>;
