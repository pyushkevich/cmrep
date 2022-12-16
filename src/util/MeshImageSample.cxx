#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkImageRegionIterator.h"
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "itkVectorImage.h"
#include "vtkThreshold.h"
#include "vtkConnectivityFilter.h"
#include "vtkSmartPointer.h"
#include "vtkGeometryFilter.h"
#include "ReadWriteVTK.h"
#include <iostream>
#include <limits>

using namespace std;

int usage()
{
  cout << "Usage: mesh_image_sample [options] mesh.vtk sample.img output.vtk array_name" << endl;
  cout << "Options: " << endl;
  cout << "   -rms i n       : Root mean square mode: recursively computes RMS of n images, e.g. " << endl;
  cout << "                    mesh_image_sample -rms 0 5 mesh.vtk sample.img mesh.vtk arr " << endl;
  cout << "                    mesh_image_sample -rms 1 5 mesh.vtk sample.img mesh.vtk arr " << endl;
  cout << "   -i n           : Interpolation, 0 for nearest neighbor, 1 for linear (*default*)" << endl;
  cout << "   -t min max     : Trim away anything in the mesh that falls outside of range" << endl;
  cout << "   -C n           : After trimming, retain n largest connected components" << endl;
  cout << "   -V             : Voting mode - each vertex gets the label of the closest " << endl;
  cout << "                    voxel with a non-zero label " << endl;
  cout << "   -B             : Write VTK files as binary" << endl;
  cout << "   -b <value>     : Background value (when vertex falls outside of the image)" << endl;
  cout << "                    defaults to NaN" << endl;
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
void WriteMesh(TMeshType *mesh, const char *fname, bool vtk_binary)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname, bool vtk_binary)
{
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(vtk_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname, bool vtk_binary)
{
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(vtk_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

vtkSmartPointer<vtkUnstructuredGrid> ThresholdMesh(
  vtkUnstructuredGrid *mesh, float trim_min, float trim_max, const char *arrname, bool trim_comp)
{
  vtkSmartPointer<vtkThreshold> thresh = vtkThreshold::New();
  thresh->SetInputData(mesh);
  thresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, arrname);
  thresh->SetLowerThreshold(trim_min);
  thresh->SetUpperThreshold(trim_max);
  thresh->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
  thresh->Update();
  vtkSmartPointer<vtkUnstructuredGrid> result = thresh->GetOutput();

  if(trim_comp)
    {
    vtkSmartPointer<vtkConnectivityFilter> conn = vtkConnectivityFilter::New();
    conn->SetInputData(result);
    conn->SetExtractionModeToLargestRegion();
    conn->Update();
    result = conn->GetUnstructuredGridOutput();
    }

  return result;
}

vtkSmartPointer<vtkPolyData> ThresholdMesh(
  vtkPolyData *mesh, float trim_min, float trim_max, const char *arrname, bool trim_comp)
{
  vtkSmartPointer<vtkThreshold> thresh = vtkThreshold::New();
  thresh->SetInputData(mesh);
  thresh->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, arrname);
  thresh->SetLowerThreshold(trim_min);
  thresh->SetUpperThreshold(trim_max);
  thresh->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
  thresh->Update();
  vtkSmartPointer<vtkUnstructuredGrid> grid = thresh->GetOutput();
  std::cout << "After thresholding " 
    << grid->GetNumberOfPoints() << " of " << mesh->GetNumberOfPoints() 
    << " remain." << std::endl;

  // Apply connected components
  if(trim_comp)
    {
    vtkSmartPointer<vtkConnectivityFilter> conn = vtkConnectivityFilter::New();
    conn->SetInputData(grid);
    conn->SetExtractionModeToLargestRegion();
    conn->Update();
    vtkSmartPointer<vtkPolyData> result = conn->GetPolyDataOutput();

    std::cout << "After connectivity filter " 
      << grid->GetNumberOfPoints() << " of " << mesh->GetNumberOfPoints() 
      << " remain." << std::endl;

    return result;
    }
  else
    {
    vtkSmartPointer<vtkGeometryFilter> geom = vtkGeometryFilter::New();
    geom->SetInputData(grid);
    geom->Update();
    vtkSmartPointer<vtkPolyData> result = geom->GetOutput();
    return result;
    }
}

template <class TMeshType>
int MeshImageSample(int argc, char *argv[], size_t irms, size_t nrms, int interp_mode,
  bool flag_trim, float trim_min, float trim_max, bool trim_comp, bool voting_mode, bool vtk_binary, float background_value)
{
  // Read the input mesh
  TMeshType *mesh = ReadMesh<TMeshType>(argv[argc-4]);
  size_t n = mesh->GetNumberOfPoints();

  // Create an array to hold the output
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetNumberOfTuples(mesh->GetNumberOfPoints());

  // Special case - voting mode
  if(voting_mode)
    {
    // Load the image as a scalar image
    typedef itk::Image<float, 3> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer fltReader = ReaderType::New();
    fltReader->SetFileName(argv[argc-3]);
    fltReader->Update();
    ImageType::Pointer imgInput = fltReader->GetOutput();

    // Find all the unique labels in the image
    int last_label = 0;
    std::set<int> unique_labels;
    for(itk::ImageRegionIterator<ImageType> it(imgInput, imgInput->GetBufferedRegion());
      !it.IsAtEnd(); ++it)
      {
      int label = (int) it.Get() + 0.5;
      if(label != last_label && label != 0)
        {
        unique_labels.insert(label);
        last_label = label;
        }
      }

    // If too many unique labels, crash
    if(unique_labels.size() > 128) 
      {
      std::cerr << "Too many unique intensity values in input image, max is 128" << std::endl;
      return -1;
      }
    else
      {
      std::cout << "Voting among " << unique_labels.size() << " unique labels " << std::endl;
      }

    // Create an array to hold the highest vote
    array->SetNumberOfComponents(1);
    vtkFloatArray *vote = vtkFloatArray::New();
    vote->SetNumberOfComponents(1);
    vote->SetNumberOfTuples(mesh->GetNumberOfPoints());
    vote->FillComponent(0, 0.0);

    bool first_label = true;
    for(auto label : unique_labels) 
      {
      // Threshhold on the label
      typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdFilter;
      typename ThresholdFilter::Pointer f_thresh = ThresholdFilter::New();
      f_thresh->SetInput(imgInput);
      f_thresh->SetLowerThreshold(label);
      f_thresh->SetUpperThreshold(label);
      f_thresh->SetInsideValue(1.0);
      f_thresh->SetOutsideValue(0.0);

      // Smooth the thresholded image
      typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> DistanceFilter;
      typename DistanceFilter::Pointer f_dist = DistanceFilter::New();
      f_dist->SetInput(f_thresh->GetOutput());
      f_dist->Update();

      // Create an interpolation function
      typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpType;
      typename InterpType::Pointer interp = InterpType::New();
      interp->SetInputImage(f_dist->GetOutput());

      // Sample the vertices
      for(unsigned int i = 0; i < n; i++)
        {
        // Get the physical position of the point
        double *x = mesh->GetPoint(i);
        itk::Point<double,3> pt;
        itk::ContinuousIndex<double,3> idx;

        // The first two coordinates are flipped for RAS/LPS conversion
        pt[0] = -x[0]; pt[1] = -x[1]; pt[2] = x[2];

        imgInput->TransformPhysicalPointToContinuousIndex(pt, idx);
        double pix = background_value;
        if(interp->IsInsideBuffer(idx)) 
          pix = interp->EvaluateAtContinuousIndex(idx);

        if(first_label || vote->GetComponent(i, 0) > pix)
          {
          vote->SetComponent(i, 0, pix);
          array->SetComponent(i, 0, label);
          }
        }
      first_label = false;
      }
    }
  else
    {
    // Read the input image
    typedef itk::VectorImage<float, 3> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer fltReader = ReaderType::New();
    fltReader->SetFileName(argv[argc-3]);
    fltReader->Update();
    ImageType::Pointer imgInput = fltReader->GetOutput();

    // Create an interpolation function
    typedef itk::InterpolateImageFunction<ImageType, double> InterpType;
    InterpType::Pointer interp;
    if(interp_mode == 0)
      interp = itk::NearestNeighborInterpolateImageFunction<ImageType, double>::New();
    else
      interp = itk::LinearInterpolateImageFunction<ImageType, double>::New();

    interp->SetInputImage(imgInput);

    // Set the output number of components
    size_t c = imgInput->GetNumberOfComponentsPerPixel();
    array->SetNumberOfComponents(c);

    // Sample the data
    for(size_t i = 0; i < n; i++)
      {
      // Get the physical position of the point
      double *x = mesh->GetPoint(i);
      itk::Point<double,3> pt;
      itk::ContinuousIndex<double,3> idx;

      // The first two coordinates are flipped for RAS/LPS conversion
      pt[0] = -x[0]; pt[1] = -x[1]; pt[2] = x[2];

      imgInput->TransformPhysicalPointToContinuousIndex(pt, idx);
      if(interp->IsInsideBuffer(idx))
        {
        typename ImageType::PixelType pix = interp->EvaluateAtContinuousIndex(idx);
        for(int j = 0; j < c; j++)
          array->SetComponent(i, j, pix[j]);
        }
      else
        {
        for(int j = 0; j < c; j++)
          array->SetComponent(i, j, background_value);
        }
      }
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

  // If trimming
  if(flag_trim) 
    {
    // Apply thresholding
    mesh = ThresholdMesh(mesh, trim_min, trim_max, arrname, trim_comp);
    }

  // Write the output
  WriteMesh<TMeshType>(mesh, argv[argc - 2], vtk_binary);
  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  if(argc < 5)
    return usage();

  // Read optional parameters
  size_t irms = 0, nrms = 0;
  int interp_mode = 1;

  // Optional trimming parameters
  float trim_min, trim_max;
  bool flag_trim = false, trim_comp = false, voting_mode = false, vtk_binary = false;
  double background_value = std::numeric_limits<double>::quiet_NaN();

  for(int ip = 1; ip < argc-4; ip++)
    {
    if(strcmp(argv[ip], "-i") == 0)
      {
      if(ip + 1 < argc-4)
        {
        interp_mode = atoi(argv[++ip]);
        cout << "Interpolation mode: " << interp_mode << endl;
        }
      else
        {
        cerr << "error: -i flag needs one parameter" << endl;
        return usage();
        }
      }
    else if(strcmp(argv[ip], "-rms") == 0)
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
    else if(strcmp(argv[ip], "-t") == 0)
      {
      if(ip + 2 < argc-4)
        {
        flag_trim = true;
        trim_min = atof(argv[++ip]);
        trim_max = atof(argv[++ip]);
        }
      else
        {
        cerr << "error: -t flag needs two parameters" << endl;
        return usage();
        }
      }
    else if(strcmp(argv[ip], "-C") == 0)
      {
      trim_comp = true;
      }
    else if(strcmp(argv[ip], "-V") == 0)
      {
      voting_mode = true;
      }
    else if(strcmp(argv[ip], "-B") == 0)
      {
      vtk_binary = true;
      }
    else if(strcmp(argv[ip], "-b") == 0)
      {
      background_value = atof(argv[++ip]);
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
    return MeshImageSample<vtkUnstructuredGrid>( argc, argv, irms, nrms, interp_mode,
      flag_trim, trim_min, trim_max, trim_comp, voting_mode, vtk_binary, background_value);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return MeshImageSample<vtkPolyData>( argc, argv, irms, nrms, interp_mode,
      flag_trim, trim_min, trim_max, trim_comp, voting_mode, vtk_binary, background_value);

    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}

