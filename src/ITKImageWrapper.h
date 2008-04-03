#ifndef __ITKImageWrapper_h_
#define __ITKImageWrapper_h_

#include "smlmath.h"

namespace itk {
  template <typename TPixel, unsigned int VDim> class Image;
};

template<typename TPixel>
class ITKImageWrapper
{
public:
  ITKImageWrapper() {};
  virtual ~ITKImageWrapper() {};
    
  // The wrapped ITK image
  typedef itk::Image<TPixel, 3> ImageType;
  
  // Load the image from a file
  virtual void LoadFromFile(const char *file) = 0;

  // Save the image to a file
  virtual void SaveToFile(const char *file) = 0;

  // Set the internal image to some ITK image, can be NULL
  virtual void SetInternalImage(ImageType *xImage) = 0;

  // Check whether an image is loaded into the wrapper
  virtual bool IsImageLoaded() = 0;

  // Get the basic image properties
  virtual unsigned int GetImageSize(unsigned int d) = 0;
  virtual double GetImageSpacing(unsigned int d) = 0;
  virtual double GetImageOrigin(unsigned int d) = 0;

  // Get a voxel corresponding to point in space
  virtual void GetVoxel(double px, double py, double pz, int &vx, int &vy, int &vz) = 0;

  // Interpolate the image at a continuous index, or return background if out of bounds
  virtual float Interpolate(double x, double y, double z, float xBackground) = 0;

  // Interpolate the image at a continuous index, or return background if out of bounds
  virtual float InterpolateNearestNeighbor(double x, double y, double z, float xBackground) = 0;

  // Interpolate the image gradient
  virtual void InterpolateGradient(const SMLVec3d &X, SMLVec3d &G) = 0;

  // Get the internal ITK image pointer
  virtual ImageType *GetInternalImage() = 0;
};

/**
 * A factory class for creating new instances of image wrappers
 */
template<class TPixel>
class ITKImageWrapperFactory
{
public:
  static ITKImageWrapper<TPixel> *NewImageWrapper();
};

#endif // __ITKImageWrapper_h_
