#include "ReadWriteVTK.h"
#include <vtkBYUReader.h>
#include <vtkBYUWriter.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>

using namespace std;

vtkPolyData *ReadVTKData(string fn)
{
  vtkPolyData *p1 = NULL;

  // Choose the reader based on extension

  if(fn.rfind(".byu") == fn.length() - 4)
    {
    vtkBYUReader *reader = vtkBYUReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".stl") == fn.length() - 4)
    {
    vtkSTLReader *reader = vtkSTLReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".vtk") == fn.length() - 4)
    {
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".obj") == fn.length() - 4)
    {
    vtkOBJReader *reader = vtkOBJReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".ply") == fn.length() - 4)
    {
    vtkPLYReader *reader = vtkPLYReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else
    {
    cout << "Could not find a reader for " << fn << endl;
    return NULL;
    }

  return p1;
}

void WriteVTKData(vtkPolyData *data, string fn, bool force_binary)
{
  if(fn.rfind(".byu") == fn.length() - 4)
    {
    vtkBYUWriter *writer = vtkBYUWriter::New();
    writer->SetGeometryFileName(fn.c_str());
    writer->SetInputData(data);
    writer->Update();
    }
  else if(fn.rfind(".stl") == fn.length() - 4)
    {
    vtkSTLWriter *writer = vtkSTLWriter::New();
    writer->SetFileName(fn.c_str());
    writer->SetInputData(data);
    writer->Update();
    }
  else if(fn.rfind(".vtk") == fn.length() - 4)
    {
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(fn.c_str());
    writer->SetInputData(data);
    if(force_binary)
      writer->SetFileTypeToBinary();
    writer->Update();
    }
  else
    {
    cerr << "Could not find a writer for " << fn << endl;
    return;
    }
}

// Templated mesh IO methods
template <class TMeshType>
vtkSmartPointer<TMeshType> ReadMesh(const char *fname)
{ return NULL; }

template <>
vtkSmartPointer<vtkUnstructuredGrid> ReadMesh<>(const char *fname)
{
  vtkNew<vtkUnstructuredGridReader> reader;
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}

template <>
vtkSmartPointer<vtkPolyData> ReadMesh<>(const char *fname)
{
  vtkNew<vtkPolyDataReader> reader;
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}

template <class TMeshType>
void WriteMesh(TMeshType *mesh, const char *fname, bool vtk_binary)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname, bool vtk_binary)
{
  vtkNew<vtkUnstructuredGridWriter> writer;
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(vtk_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname, bool vtk_binary)
{
  vtkNew<vtkPolyDataWriter> writer;
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(vtk_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}
