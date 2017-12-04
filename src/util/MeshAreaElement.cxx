#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkSmartPointer.h"
#include <vtkPointData.h>
#include "vtkTriangle.h"

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
  writer->SetInputData(mesh);
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  writer->Update();
}

int usage()
{
  cout << "Usage: mesh_area_element input.vtk output.vtk" << endl;
  return -1;
}


template <class TMeshType>
int
MeshAreaElement(const char *fin, const char *fout)
{
  // Read the mesh
  vtkSmartPointer<TMeshType> mesh = ReadMesh<TMeshType>(fin);
  int n = mesh->GetNumberOfPoints();

  // Create area element array
  vtkSmartPointer<vtkFloatArray> daArea = vtkFloatArray::New();
  daArea->SetNumberOfComponents(1);
  daArea->SetNumberOfTuples(n);
  daArea->SetName("AreaElement");
  daArea->FillComponent(0, 0.0);
  mesh->GetPointData()->AddArray(daArea);

  // Compute area element
  for(int k = 0; k < mesh->GetNumberOfCells(); k++)
    {
    vtkCell *cell = mesh->GetCell(k);
    if(cell->GetCellType() != VTK_TRIANGLE)
      {
      std::cerr << "Wrong cell type, should be VTK_TRIANGLE" << std::endl;
      return -1;
      }

    // Compute the area of the triangle
    double p0[3], p1[3], p2[3];
    vtkIdType a0 = cell->GetPointId(0), a1 = cell->GetPointId(1), a2 = cell->GetPointId(2);
    mesh->GetPoint(a0, p0); mesh->GetPoint(a1, p1); mesh->GetPoint(a2, p2); 
    double area = vtkTriangle::TriangleArea(p0, p1, p2);

    // If area is negative complain
    if(area < 0)
      {
      std::cerr << "Negative cell area, suggests flipped triangles" << std::endl;
      return -1;
      }

    // Split the volume between neighbors
    daArea->SetTuple1(a0, area / 3.0 + daArea->GetTuple1(a0));
    daArea->SetTuple1(a1, area / 3.0 + daArea->GetTuple1(a1));
    daArea->SetTuple1(a2, area / 3.0 + daArea->GetTuple1(a2));
    }

  // Write the mesh
  WriteMesh<TMeshType>(mesh, fout);

  return 0;
}

int main(int argc, char **argv)
{
  if(argc < 3)
    return usage();

  return MeshAreaElement<vtkPolyData>(argv[1], argv[2]);
}
