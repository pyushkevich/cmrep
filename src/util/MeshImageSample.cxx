#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkOrientedRASImage.h"
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "ReadWriteVTK.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "Usage: mesh_image_sample [options] mesh.vtk sample.img output.vtk array_name" << endl;
  cout << "Options: " << endl;
  cout << "   -rms i n       : Root mean square mode: recursively computes RMS of n images, e.g. " << endl;
  cout << "                    mesh_image_sample -rms 0 5 mesh.vtk sample.img mesh.vtk arr " << endl;
  cout << "                    mesh_image_sample -rms 1 5 mesh.vtk sample.img mesh.vtk arr " << endl;
  return -1;
}

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
int MeshImageSample(int argc, char *argv[], size_t irms, size_t nrms)
{
  // Read the input mesh
  TMeshType *mesh = ReadMesh<TMeshType>(argv[argc-4]);

  // Read the input image
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[argc-3]);
  fltReader->Update();
  ImageType::Pointer imgInput = fltReader->GetOutput();

  // Create an interpolation function
  //typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpType;
  InterpType::Pointer interp = InterpType::New();
  interp->SetInputImage(imgInput);

  // Create an array to hold the output
  vtkFloatArray *array = vtkFloatArray::New();
  array->Allocate(mesh->GetNumberOfPoints());

  // Sample the data
  size_t n = mesh->GetNumberOfPoints();
  for(size_t i = 0; i < n; i++)
    {
    double *x = mesh->GetPoint(i);
    itk::Point<double,3> pt;
    itk::ContinuousIndex<double,3> idx;
    pt[0] = x[0]; pt[1] = x[1]; pt[2] = x[2];
    imgInput->TransformRASPhysicalPointToContinuousIndex(pt, idx);
    float val = interp->EvaluateAtContinuousIndex(idx);
    array->InsertTuple1(i, val);
    }

  // Get the point data
  vtkPointData *pdat = mesh->GetPointData();
  char *arrname = argv[argc-1];

  // If we are in RMS mode 
  if(irms > 0)
    {
    // Get the array from the polydata
    vtkDataArray *arrold = pdat->GetArray(arrname);
    if(arrold == NULL)
      {
      cerr << "RMS with non-zero first parameter requires array " << arrname << " in the mesh" << endl;
      return usage();
      }

    // Add the new data to the old data
    for(size_t i = 0; i < n; i++)
      {
      double dsum = array->GetTuple1(i) * array->GetTuple1(i) + arrold->GetTuple1(i);
      if(irms == nrms - 1)
        dsum = sqrt(dsum / nrms);
      array->SetTuple1(i, dsum);
      }

    // Remove the array
    pdat->RemoveArray(arrname);
    }
  else if(nrms > 0)
    {
    for(size_t i = 0; i < n; i++)
      array->SetTuple1(i, array->GetTuple1(i) * array->GetTuple1(i));
    }

  // Stick the array in the output
  array->SetName(arrname);
  pdat->AddArray(array);

  // Write the output
  WriteMesh<TMeshType>(mesh, argv[argc - 2]);
  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  if(argc < 5)
    return usage();

  // Read optional parameters
  size_t irms = 0, nrms = 0;
  for(size_t ip = 1; ip < argc-4; ip++)
    {
    if(strcmp(argv[ip], "-rms") == 0)
      {
      if(ip + 2 < argc-4)
        {
        irms = atoi(argv[++ip]);
        nrms = atoi(argv[++ip]);
        }
      else
        {
        cerr << "error: -rms flag needs two parameters" << endl;
        return usage();
        }
      }
    else
      {
      cerr << "error: unrecognized parameter " << argv[ip] << endl;
      return usage();
      }
    }
  
  // Read the meshes
  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(argv[argc-4]);
  reader->OpenVTKFile();
  reader->ReadHeader();

  bool isPolyData = true;
  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    isPolyData = false;
    return MeshImageSample<vtkUnstructuredGrid>( argc, argv, irms, nrms);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return MeshImageSample<vtkPolyData>( argc, argv, irms, nrms);

    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}

