#include "ScriptInterface.h"
#include "ITKImageWrapper.h"

#include <itkOrientedRASImage.h>
#include <itkCastImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkDerivativeImageFilter.h>
#include <vnl/vnl_erf.h>

#include <string>

// ITK includes 
#include "itk_to_nifti_xform.h"
#include <itkOrientedRASImage.h>
#include <itkVTKImageExport.h>

// VTK includes
#include <vtkImageData.h>
#include <vtkImageImport.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>

using namespace std;

FloatImage::FloatImage()
{
  xImage = ITKImageWrapperFactory<float>::NewImageWrapper();
  xGradient[0] = ITKImageWrapperFactory<float>::NewImageWrapper();
  xGradient[1] = ITKImageWrapperFactory<float>::NewImageWrapper();
  xGradient[2] = ITKImageWrapperFactory<float>::NewImageWrapper();
}
  
FloatImage::~FloatImage()
{
  delete xImage;
  delete xGradient[0];
  delete xGradient[1];
  delete xGradient[2];
}

void FloatImage::GetVoxelIndex(const SMLVec3d &x, int &vx, int &vy, int &vz)
{
  xImage->GetVoxel(x[0], x[1], x[2], vx, vy, vz);

}

// Interpolate the image at a given position
float FloatImage::Interpolate(const SMLVec3d &x)
{ 
  return xImage->Interpolate(x[0], x[1], x[2], xOutsideValue); 
}

// Interpolate the image at a given position
float FloatImage::InterpolateNearestNeighbor(const SMLVec3d &x)
{ 
  return xImage->InterpolateNearestNeighbor(x[0], x[1], x[2], xOutsideValue); 
}

// Interpolate the image gradient
void FloatImage::InterpolateImageGradient(const SMLVec3d &x, SMLVec3d &g)
{
  xImage->InterpolateGradient(x, g);
  // g[0] = xGradient[0]->Interpolate(x[0], x[1], x[2], 0.0f);
  // g[1] = xGradient[1]->Interpolate(x[0], x[1], x[2], 0.0f);
  // g[2] = xGradient[2]->Interpolate(x[0], x[1], x[2], 0.0f);
}




void FloatImage::LoadFromFile(const char *file)
{
  // Load an image from ITK
  xImage->LoadFromFile(file);

  // Clear the gradient images
  for(unsigned int i = 0; i < 3; i++)
    xGradient[i]->SetInternalImage(NULL);
}

void FloatImage::LoadGradientFromFile(unsigned int iComponent, const char *file)
{
  // Load an image from ITK
  xGradient[iComponent]->LoadFromFile(file);
}

// Load the image and the gradient from a path
void FloatImage::LoadFromPath(const char *filebase, const char *ext)
{
  string fMain = string(filebase) + "." + string(ext);
  string fGradX = string(filebase) + ".dx." + string(ext);
  string fGradY = string(filebase) + ".dy." + string(ext);
  string fGradZ = string(filebase) + ".dz." + string(ext);

  LoadFromFile(fMain.c_str());
  LoadGradientFromFile(0, fGradX.c_str());
  LoadGradientFromFile(1, fGradY.c_str());
  LoadGradientFromFile(2, fGradZ.c_str());
}

void FloatImage::SaveToFile(const char *file)
{
  // Save the ITK image
  xImage->SaveToFile(file);
}

void FloatImage::SaveGradientToFile(unsigned int iComponent, const char *file)
{
  // Save the ITK image
  xGradient[iComponent]->SaveToFile(file);
}

  // Load the image and the gradient from a path
void FloatImage::SaveToPath(const char *filebase, const char *ext)
{
  string fMain = string(filebase) + "." + string(ext);
  string fGradX = string(filebase) + ".dx." + string(ext);
  string fGradY = string(filebase) + ".dy." + string(ext);
  string fGradZ = string(filebase) + ".dz." + string(ext);

  SaveToFile(fMain.c_str());
  SaveGradientToFile(0, fGradX.c_str());
  SaveGradientToFile(1, fGradY.c_str());
  SaveGradientToFile(2, fGradZ.c_str());
}

void FloatImage::SetToGradientMagnitude(BinaryImage *imgSource, double xSigma)
{
  // Blur the image with a Gaussian
  typedef itk::OrientedRASImage<unsigned char, 3> BinaryImageType;
  typedef itk::OrientedRASImage<float, 3> FloatImageType;
  typedef itk::CastImageFilter<BinaryImageType, FloatImageType> CasterType;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>
    GaussianType;
  typedef itk::GradientMagnitudeImageFilter<FloatImageType, FloatImageType>
    GradMagType;

  // Create a pipeline of filters
  CasterType::Pointer fltCaster = CasterType::New();
  fltCaster->SetInput(imgSource->xImage->GetInternalImage());

  GaussianType::Pointer fltGaussian = GaussianType::New();
  fltGaussian->SetInput(fltCaster->GetOutput());
  fltGaussian->SetVariance(xSigma * xSigma);

  GradMagType::Pointer fltGrad = GradMagType::New();
  fltGrad->SetInput(fltGaussian->GetOutput());
  fltGrad->Update();

  // Set the internal image
  xImage->SetInternalImage(fltGrad->GetOutput());
}

bool FloatImage::IsGradientAvailable()
{
  return xGradient[0]->IsImageLoaded();
}

void FloatImage::SetOutsideValue(float xOutsideValue)
{
  this->xOutsideValue = xOutsideValue;
}

void FloatImage::SetImageOrigin(double ox, double oy, double oz)
{
  double origin[] = {ox, oy, oz};
  if(xImage->GetInternalImage())
    xImage->GetInternalImage()->SetOrigin(origin);
  for(size_t i = 0; i < 3; i++)
    if(xGradient[i]->GetInternalImage())
      xGradient[i]->GetInternalImage()->SetOrigin(origin);
}

class BinaryFloatFunctor 
{
public:
  float operator() (unsigned char x)
    { return (x == 0) ? -1.0f : 1.0f; } 
};

class ErrorFunctor
{
public:
  ErrorFunctor(float alpha = 1.0f)
    { this->alpha = alpha; }

  float operator() (float x)
    { return vnl_erf( alpha * x ); }

  bool operator == (const ErrorFunctor &x) const
    { return x.alpha == alpha; }

  bool operator != (const ErrorFunctor &x) const
    { return !(*this == x); }

private:
  float alpha;
};

void FloatImage
::SetToBlurredBinary(BinaryImage *imgSource, double xSigma, double xErrorScale)
{
  // Blur the image with a Gaussian
  typedef itk::OrientedRASImage<unsigned char, 3> BinaryImageType;
  typedef itk::OrientedRASImage<float, 3> FloatImageType;
  typedef itk::UnaryFunctorImageFilter<
    BinaryImageType, FloatImageType, BinaryFloatFunctor> MapperType;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>
    GaussianType;
  typedef itk::DerivativeImageFilter<FloatImageType, FloatImageType>
    DerivativeType;

  // Create a pipeline of filters to blur the image
  MapperType::Pointer fltMapper = MapperType::New();
  fltMapper->SetInput(imgSource->xImage->GetInternalImage());

  GaussianType::Pointer fltGaussian = GaussianType::New();
  fltGaussian->SetInput(fltMapper->GetOutput());
  fltGaussian->SetVariance(xSigma * xSigma);
  fltGaussian->SetMaximumError(0.001);
  fltGaussian->SetUseImageSpacing(true);
  fltGaussian->Update();

  // Optionally apply the error functor
  if(xErrorScale > 0.0)
    {
    typedef itk::UnaryFunctorImageFilter<
      FloatImageType, FloatImageType, ErrorFunctor> ErrorMapper;
    ErrorMapper::Pointer fltError = ErrorMapper::New();
    fltError->SetInput(fltGaussian->GetOutput());
    fltError->SetFunctor( ErrorFunctor( xErrorScale ) );
    fltError->Update();
    xImage->SetInternalImage(fltError->GetOutput());
    }
  else
    xImage->SetInternalImage(fltGaussian->GetOutput());
  
  // Now compute the gradient magnitude by convolution with a Gaussian directional
  // derivative filter in each cardinal direction
  for(unsigned int d = 0; d < 3; d++)
    {
    DerivativeType::Pointer fltDerivative = DerivativeType::New();
    fltDerivative->SetInput(fltGaussian->GetOutput());
    fltDerivative->SetOrder( 1 );
    fltDerivative->SetDirection( d );
    fltDerivative->Update();
    xGradient[d]->SetInternalImage(fltDerivative->GetOutput());
    }
}

double FloatImage::IntegratePositiveVoxels()
{
  // Create an iterator
  typedef itk::OrientedRASImage<float, 3> FloatImageType;
  typedef itk::ImageRegionConstIterator<FloatImageType> IteratorType;

  // Take the sum of all positive voxels in the image
  double xPositiveSum = 0.0;

  FloatImageType::Pointer img = xImage->GetInternalImage();
  IteratorType it(img, img->GetBufferedRegion());
  while(!it.IsAtEnd())
    {
    if(it.Get() > 0)
      xPositiveSum += it.Get();
    ++it;
    }

  // Scale the positive sum by the volume of the voxel
  double xVoxelVolume = 
    img->GetSpacing()[0] * img->GetSpacing()[1] * img->GetSpacing()[2];

  // Return the total probability of the positive region
  return xPositiveSum * xVoxelVolume;
}

double FloatImage::ComputeImageVolume()
{
  // Scale the positive sum by the volume of the voxel
  typedef itk::OrientedRASImage<float, 3> FloatImageType;
  FloatImageType::Pointer img = xImage->GetInternalImage();
  double xVoxelVolume = 
    img->GetSpacing()[0] * img->GetSpacing()[1] * img->GetSpacing()[2];

  return xVoxelVolume * img->GetBufferedRegion().GetNumberOfPixels();
}

double FloatImage::ComputeObjectVolume()
{
  // Create an iterator
  typedef itk::OrientedRASImage<float, 3> FloatImageType;
  typedef itk::ImageRegionConstIterator<FloatImageType> IteratorType;

  // Create the volume accumulator
  long nPositiveVoxels = 0;
  double xVolume = 0.0;

  // Add each voxel to the accumulator
  FloatImageType::Pointer img = xImage->GetInternalImage();
  IteratorType it(img, img->GetBufferedRegion());
  while(!it.IsAtEnd())
    {
    if(it.Get() > 0)
      nPositiveVoxels++;
    xVolume += 0.5 * (it.Get() + 1.0);
    ++it;
    }

  // Scale by the voxel's spacing
  double xVoxelVolume = 
    img->GetSpacing()[0] * img->GetSpacing()[1] * img->GetSpacing()[2];
  cout << "NPositive = " << nPositiveVoxels << endl;
  // return nPositiveVoxels * xVoxelVolume;
  return xVolume * xVoxelVolume;
}

BinaryImage::BinaryImage()
{
  xImage = ITKImageWrapperFactory<unsigned char>::NewImageWrapper();
}
  
BinaryImage::~BinaryImage()
{
  delete xImage;
}

void BinaryImage::LoadFromFile(const char *file)
{
  // Load an image from ITK
  xImage->LoadFromFile(file);
}


vtkSmartPointer<vtkPolyData> 
GenerateContour(FloatImage *image)
{
  typedef itk::OrientedRASImage<float,3> ImageType;
  typedef itk::VTKImageExport<ImageType> VTKExportType;
  typedef itk::SmartPointer<VTKExportType> VTKExportPointer;
  typedef ImageType::TransformMatrixType TransformMatrixType; 
  VTKExportPointer m_VTKExporter;

  vtkSmartPointer<vtkImageImport> m_VTKImporter;
  vtkSmartPointer<vtkMarchingCubes> m_MarchingCubesFilter;
  vtkSmartPointer<vtkTransformPolyDataFilter> m_TransformFilter;
  vtkSmartPointer<vtkTransform> m_Transform;

  // Initialize the VTK Exporter
  m_VTKExporter = VTKExportType::New();
  m_VTKExporter->ReleaseDataFlagOn();
  m_VTKExporter->SetInput(image->GetInternalImage()->GetInternalImage());
  
  // Initialize the VTK Importer
  m_VTKImporter = vtkSmartPointer<vtkImageImport>::New();
  m_VTKImporter->ReleaseDataFlagOn();

  // Pipe the importer into the exporter (that's a lot of code)
  m_VTKImporter->SetUpdateInformationCallback(
    m_VTKExporter->GetUpdateInformationCallback());
  m_VTKImporter->SetPipelineModifiedCallback(
    m_VTKExporter->GetPipelineModifiedCallback());
  m_VTKImporter->SetWholeExtentCallback(
    m_VTKExporter->GetWholeExtentCallback());
  m_VTKImporter->SetSpacingCallback(
    m_VTKExporter->GetSpacingCallback());
  m_VTKImporter->SetOriginCallback(
    m_VTKExporter->GetOriginCallback());
  m_VTKImporter->SetScalarTypeCallback(
    m_VTKExporter->GetScalarTypeCallback());
  m_VTKImporter->SetNumberOfComponentsCallback(
    m_VTKExporter->GetNumberOfComponentsCallback());
  m_VTKImporter->SetPropagateUpdateExtentCallback(
    m_VTKExporter->GetPropagateUpdateExtentCallback());
  m_VTKImporter->SetUpdateDataCallback(
    m_VTKExporter->GetUpdateDataCallback());
  m_VTKImporter->SetDataExtentCallback(
    m_VTKExporter->GetDataExtentCallback());
  m_VTKImporter->SetBufferPointerCallback(
    m_VTKExporter->GetBufferPointerCallback());  
  m_VTKImporter->SetCallbackUserData(
    m_VTKExporter->GetCallbackUserData());

  // Create and configure the marching cubes filter
  m_MarchingCubesFilter = vtkSmartPointer<vtkMarchingCubes>::New();
  m_MarchingCubesFilter->ReleaseDataFlagOn();
  m_MarchingCubesFilter->ComputeScalarsOff();
  m_MarchingCubesFilter->ComputeGradientsOff();
  m_MarchingCubesFilter->SetNumberOfContours(1);
  m_MarchingCubesFilter->SetValue(0,0.0f);
  m_MarchingCubesFilter->SetInputConnection(m_VTKImporter->GetOutputPort());

  // Create a transform into RAS coordinates
  vnl_matrix_fixed<double, 4, 4> vtk2nii =
    ConstructVTKtoNiftiTransform(
      image->GetInternalImage()->GetInternalImage()->GetDirection().GetVnlMatrix(),
      image->GetInternalImage()->GetInternalImage()->GetOrigin().GetVnlVector(),
      image->GetInternalImage()->GetInternalImage()->GetSpacing().GetVnlVector());          

  m_Transform = vtkSmartPointer<vtkTransform>::New();
  m_Transform->SetMatrix(vtk2nii.data_block());

  m_TransformFilter = vtkTransformPolyDataFilter::New();
  m_TransformFilter->SetTransform(m_Transform);
  m_TransformFilter->SetInputConnection(m_MarchingCubesFilter->GetOutputPort());
  m_TransformFilter->Update();

  cout << "Transform: " << m_Transform << endl;
  cout << "VOX2RAS: " << image->GetInternalImage()->GetInternalImage()->
    GetVoxelSpaceToRASPhysicalSpaceMatrix() << endl;

  cout << "SOS2RAS: " << image->GetInternalImage()->GetInternalImage()->
    GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix() << endl;


  vtkSmartPointer<vtkPolyData> result = m_TransformFilter->GetOutput();
  return result;
}


