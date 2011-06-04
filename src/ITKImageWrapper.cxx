#include "ITKImageWrapper.h"

#include "itkOrientedRASImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "vnl/vnl_inverse.h"

template<typename TPixel>
class ITKImageWrapperImpl : public virtual ITKImageWrapper<TPixel>
{
public:
  // The typedefs for this class
  typedef ITKImageWrapper<TPixel> Superclass;
  typedef typename Superclass::ImageType ImageType;
  // typedef itk::LinearInterpolateImageFunction<ImageType, float> InterpolatorType;
  typedef itk::BSplineInterpolateImageFunction<ImageType> InterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NNInterpolator;
  
  // Load the image from a file
  void LoadFromFile(const char *file);

  // Save the image to a file
  void SaveToFile(const char *file);

  // Set the internal image to some ITK image
  void SetInternalImage(ImageType *xImage) 
    {
    if(xImage == NULL)
      {
      this->xImage = ImageType::New();
      }
    else
      {
      this->xImage = xImage;
      OnImageUpdate();
      }
    }

  // Check whether an image is loaded into the wrapper
  virtual bool IsImageLoaded()
    { return xImage->GetBufferedRegion().GetSize()[0] > 0; }

  // Get the basic image properties
  unsigned int GetImageSize(unsigned int d)
    { return (unsigned int) xImage->GetBufferedRegion().GetSize()[d]; }
  
  double GetImageSpacing(unsigned int d)
    { return xImage->GetSpacing()[d]; }
  
  double GetImageOrigin(unsigned int d)
    { return xImage->GetOrigin()[d]; }

  // Interpolate the image at a continuous index, or return background if out of bounds
  float Interpolate(double x, double y, double z, float xBackground);

  // Interpolate the image at a continuous index, or return background if out of bounds
  float InterpolateNearestNeighbor(double x, double y, double z, float xBackground);

  // Interpolate the gradient of the image
  void InterpolateGradient(const SMLVec3d &X, SMLVec3d &G);

  // Get a voxel for a point in space
  void GetVoxel(double px, double py, double pz, int &vx, int &vy, int &vz);

  // Get the internal ITK image pointer
  ImageType *GetInternalImage()
    { return xImage.GetPointer(); }

private:
  // The constructor
  ITKImageWrapperImpl();
  
  // The internaly stored image
  typename ImageType::Pointer xImage;

  // The linear interpolator for the currently stored image
  typename InterpolatorType::Pointer fnInterpolator;
  typename NNInterpolator::Pointer fnInterpolatorNN;

  // Function called when the image is updated externally
  void OnImageUpdate();

  // Factory is our friend
  friend class ITKImageWrapperFactory<TPixel>;
};

template<typename TPixel>
ITKImageWrapperImpl<TPixel>
::ITKImageWrapperImpl()
{
  xImage = ImageType::New();
  fnInterpolator = InterpolatorType::New();
  fnInterpolatorNN = NNInterpolator::New();

  // CHANGED BY PAUL 6/17/2010 - Why do we need cubic interpolation?
  fnInterpolator->SetSplineOrder(1);
}

template<typename TPixel>
void
ITKImageWrapperImpl<TPixel>
::GetVoxel(double px, double py, double pz, int &vx, int &vy, int &vz)
{
  // Create a point with the passed in coordinates
  itk::Point<double, 3> xPoint;
  xPoint[0] = px; xPoint[1] = py; xPoint[2] = pz;

  // Create a continuous index from the point
  itk::Index<3> xIndex;
  xImage->TransformRASPhysicalPointToIndex(xPoint, xIndex);

  // Interpolate at the index
  if(fnInterpolator->IsInsideBuffer(xIndex))
    {
    vx = (int) xIndex[0];
    vy = (int) xIndex[1];
    vz = (int) xIndex[2];
    }
  else
    {
    vx = vy = vz = -1;
    }
}


template<typename TPixel>
float
ITKImageWrapperImpl<TPixel>
::Interpolate(double x, double y, double z, float xBackground)
{
  // Create a point with the passed in coordinates
  itk::Point<double, 3> xPoint;
  xPoint[0] = x; xPoint[1] = y; xPoint[2] = z;

  // Create a continuous index from the point
  itk::ContinuousIndex<double, 3> xIndex;
  xImage->TransformRASPhysicalPointToContinuousIndex(xPoint, xIndex);

  // Interpolate at the index
  if(fnInterpolator->IsInsideBuffer(xIndex))
    return (float) fnInterpolator->EvaluateAtContinuousIndex(xIndex);
  else
    return xBackground;
}

template<typename TPixel>
void
ITKImageWrapperImpl<TPixel>
::InterpolateGradient(const SMLVec3d &X, SMLVec3d &G)
{
  // Create a point with the passed in coordinates
  itk::Point<double, 3> xPoint;
  xPoint[0] = X[0]; xPoint[1] = X[1]; xPoint[2] = X[2];

  // Create a continuous index from the point
  itk::ContinuousIndex<double, 3> xIndex;
  xImage->TransformRASPhysicalPointToContinuousIndex(xPoint, xIndex);

  // Interpolate at the index
  if(fnInterpolator->IsInsideBuffer(xIndex))
    {
    typename InterpolatorType::CovariantVectorType grad = 
      fnInterpolator->EvaluateDerivativeAtContinuousIndex(xIndex);

    // The above method applies Jacobian based on the ITK Physical-Voxel mapping
    // To convert to the nIFTI mapping we use, multiply by -1, -1, 1
    G[0] = -grad[0]; G[1] = -grad[1]; G[2] = grad[2];

    /* 

    // test derivative
    SMLVec3d Xd; double eps = 0.001;
    for(size_t d = 0; d < 3; d++)
      {
      SMLVec3d X1 = X, X2 = X;
      X1[d] -= eps; X2[d] += eps;
      Xd[d] = (this->Interpolate(X2[0],X2[1],X2[2],0.0) - this->Interpolate(X1[0],X1[1],X1[2],0.0)) / (2 * eps);
      }
    std::cout << "Analytic = " << G << "\t Numeric = " << Xd << std::endl; */
    }
  else
    { G.fill(0.0); }

  // Apply jacobian of the image coordinate mapping to the gradient
}

template<typename TPixel>
float
ITKImageWrapperImpl<TPixel>
::InterpolateNearestNeighbor(double x, double y, double z, float xBackground)
{
  // Create a point with the passed in coordinates
  itk::Point<double, 3> xPoint;
  xPoint[0] = x; xPoint[1] = y; xPoint[2] = z;

  // Create a continuous index from the point
  itk::ContinuousIndex<double, 3> xIndex;
  xImage->TransformRASPhysicalPointToContinuousIndex(xPoint, xIndex);

  // Interpolate at the index
  if(fnInterpolatorNN->IsInsideBuffer(xIndex))
    return (float) fnInterpolatorNN->EvaluateAtContinuousIndex(xIndex);
  else
    return xBackground;
}

template<typename TPixel>
void ITKImageWrapperImpl<TPixel>
::OnImageUpdate()
{
  fnInterpolator->SetInputImage(xImage);
  fnInterpolatorNN->SetInputImage(xImage);
}

template<typename TPixel>
void ITKImageWrapperImpl<TPixel>
::LoadFromFile(const char *file)
{
  // Load the image
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(file);
  fltReader->Update();
  xImage = fltReader->GetOutput();

  // Refresh other components
  OnImageUpdate();
}

template<typename TPixel>
void ITKImageWrapperImpl<TPixel>
::SaveToFile(const char *file)
{
  // Load the image
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer fltWriter = WriterType::New();
  fltWriter->SetInput(xImage);
  fltWriter->SetFileName(file);
  fltWriter->Update();
}

template <typename TPixel>
ITKImageWrapper<TPixel> *
ITKImageWrapperFactory<TPixel>
::NewImageWrapper()
{
  // Create a new instance of image wrapper impl
  return new ITKImageWrapperImpl<TPixel>;
}

// Define image wrappers for necessary types
template class ITKImageWrapperFactory<float>;
template class ITKImageWrapperFactory<unsigned char>;

