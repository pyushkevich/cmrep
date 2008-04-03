#ifndef __CoefficientMask_h_
#define __CoefficientMask_h_

#include <vector>
#include <smlmath.h>
#include "BasisFunctions2D.h"
#include "PrincipalComponents.h"

class PrincipalComponents;

using namespace std;


/**
 * A medial coefficient mask is a mapping from one set of coefficients to
 * another set of coefficients. The second set of coefficients is what defines
 * a medial model, while the first set is what we optimize over. An example is
 * optimization over the affine transforms. The input set of coefficients
 * are the parameters of the affine transform and the output set of
 * coefficients are the parameters of the medial model.
 */
class IMedialCoefficientMask : public ICoefficientSettable
{
public:
  typedef vnl_vector<double> Vec;
  
  /**
   * This method returns a variation (direction in output coefficient space)
   * that corresponds to a single component in the input coefficient space.
   */
  virtual Vec GetVariationPerComponent(size_t iComponent) 
    { throw std::bad_exception; }
  
  /**
   * This method returns a variation (direction in output coefficient space)
   * that corresponds to a variation in the input coefficient space.
   */
  virtual Vec GetVariationPerVariation(const double *xVar) 
    { throw std::bad_exception; }

};

class SelectionMedialCoefficientMask : public IMedialCoefficientMask
{
public:
  SelectionMedialCoefficientMask(
    ICoefficientSettable *source, const vector<size_t> &mask);

  SelectionMedialCoefficientMask(
    ICoefficientSettable *source, size_t iMaskSize, const size_t *iMaskIndex);

  double *GetCoefficientArray();
  void SetCoefficientArray(const double *xData);
  size_t GetNumberOfCoefficients() const
    { return nCoeff; }
  double GetCoefficient(size_t i) const; 
  void SetCoefficient(size_t i, double x);


  // The methods for returning a variation
  typedef IMedialCoefficientMask::Vec Vec;
  virtual Vec GetVariationPerComponent(size_t iComponent) = 0;
  virtual Vec GetVariationPerVariation(const double *xVar) = 0;

  // IHyperSurface2D *GetComponentSurface(size_t iCoefficient);
  // void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  // IHyperSurface2D *GetVariationSurface(const double *xData);
  // void ReleaseVariationSurface(IHyperSurface2D *xSurface);

private:
  ICoefficientSettable *xSource;
  vnl_vector<double> xRawCoefficients, xMaskCoefficients;
  vector<size_t> xMask;
  size_t nCoeff;
};

class PassThroughCoefficientMask : public virtual IMedialCoefficientMask
{
public:
  PassThroughCoefficientMask(ICoefficientSettable *source);
  double *GetCoefficientArray();
  void SetCoefficientArray(const double *xData);
  size_t GetNumberOfCoefficients() const;
  double GetCoefficient(size_t i) const;
  void SetCoefficient(size_t i, double x);

  // Get a surface corresponding to a single component
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  // Get a surface corresponding to some variation
  IHyperSurface2D *GetVariationSurface(const double *xData);
  void ReleaseVariationSurface(IHyperSurface2D *xSurface);
  
private:
  ICoefficientSettable *xSource;
};

class AffineTransformCoefficientMask : virtual public IMedialCoefficientMask
{
public:
  AffineTransformCoefficientMask(IBasisRepresentation2D *surface);

  double *GetCoefficientArray()
    { return xData.data_block(); }
  
  void SetCoefficientArray(const double *inData);
  
  size_t GetNumberOfCoefficients() const
    { return nCoeff; }

  double GetCoefficient(size_t i) const
    { return xData[i]; }

  void SetCoefficient(size_t i, double x)
    { xData[i] = x; SetCoefficientArray(xData.data_block()); }

  // Get a surface corresponding to a single component
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  // Get a surface corresponding to some variation
  IHyperSurface2D *GetVariationSurface(const double *xData);
  void ReleaseVariationSurface(IHyperSurface2D *xSurface);

protected:

  // This virtual method is called before applying an affine transform to set 
  // the state of the coefficients to the pre-transform values.
  virtual double *GetCoefficientsBeforeTransform()
    { return xOriginalCoefficients.data_block(); }

  // Create an affine variation based on some particular variation of the
  // underlying surface
  virtual IHyperSurface2D *CreateAffineVariation(
    IHyperSurface2D *xSurfaceVariation, const double *xAffineVariation);

  
  // The underlying surface
  IBasisRepresentation2D *xSurface;
  
protected:
  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  
  size_t nDim, iIndexB, nCoeff;
  
  VectorType xOriginalCoefficients;
  MatrixType A;
  VectorType b, c;
  VectorType xData;
  double xRhoScale;
};


/** 
 * This is an affine transform in 3x1D, it uses the standard 12 rather than
 * 20 coefficients. It's basically a selection mask wrapped around an affine
 * mask
 */
class AffineTransform3DCoefficientMask :
  public SelectionMedialCoefficientMask
{
public:
  AffineTransform3DCoefficientMask(IBasisRepresentation2D *surface) : 
    SelectionMedialCoefficientMask(
      xAffineMask = new AffineTransformCoefficientMask(surface), 12, xIndexArray) 
  { }

  ~AffineTransform3DCoefficientMask()
    { delete xAffineMask; }
    
private:
  AffineTransformCoefficientMask *xAffineMask;
  static const size_t xIndexArray[];
  friend class Affine3DAndPCACoefficientMask;
};

/**
 * This is a PCA-based coefficient mask. In it, PCA-space coefficients are
 * mapped to changes in the actual coefficients
 */
class PCACoefficientMask : public IMedialCoefficientMask
{
public:
  // Constructor: takes a surface and a data matrix
  PCACoefficientMask(
    const vnl_matrix<double> &mPCA, 
    ICoefficientSettable *xSurface, size_t nModes);
  
  double *GetCoefficientArray();
  void SetCoefficientArray(const double *xData);
  size_t GetNumberOfCoefficients() const
    { return nModes; }
  double GetCoefficient(size_t i) const; 
  void SetCoefficient(size_t i, double x);

  // Get a surface corresponding to a single component
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  // Get a surface corresponding to some variation
  IHyperSurface2D *GetVariationSurface(const double *xData);
  void ReleaseVariationSurface(IHyperSurface2D *xSurface);

private:
  ICoefficientSettable *xSurface;
  PrincipalComponents pca;
  size_t nModes;
  vnl_vector<double> xCoeffArray, xOriginalOffsetFromMean;
};

class AffineAndPCACoefficientMask : public AffineTransformCoefficientMask
{
public:
  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  // Constructor: takes a surface and a data matrix
  AffineAndPCACoefficientMask( const MatrixType &mPCA, 
    IBasisRepresentation2D *xSurface, size_t nModes);
  
  // Get the array of coefficients
  double *GetCoefficientArray()
    { return xCoeffArray.data_block(); }

  // Set the array of coefficients
  void SetCoefficientArray(const double *xData);

  // This method does nothing
  virtual double *GetCoefficientsBeforeTransform()
    { return xPreAffineCoeffArray.data_block(); }

  // Get the number of coefficients
  size_t GetNumberOfCoefficients() const
    { return n; }

  // Get a single coefficient
  double GetCoefficient(size_t i) const
    { return xCoeffArray[i]; }
  
  // Set a single coefficient
  void SetCoefficient(size_t i, double x)
    { xCoeffArray[i] = x; SetCoefficientArray(xCoeffArray.data_block()); }

  // Get a surface corresponding to a single component
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  // Get a surface corresponding to some variation
  IHyperSurface2D *GetVariationSurface(const double *xData);
  void ReleaseVariationSurface(IHyperSurface2D *xSurface);

private:

  // Numbers of coefficients of each kind
  size_t n, nAffine, nPCA;
  VectorType xCoeffArray, xPreAffineCoeffArray;
  PCACoefficientMask maskPCA;
};

/** 
 * This is a wrapper around AffineAndPCACoefficientMask that takes out the 
 * unused affine coefficients by defining a selection mask.
 */
class Affine3DAndPCACoefficientMask :
  public SelectionMedialCoefficientMask
{
public:
 
  Affine3DAndPCACoefficientMask(const vnl_matrix<double> &mPCA, 
    IBasisRepresentation2D *xSurface, size_t nModes);

  ~Affine3DAndPCACoefficientMask()
    { delete xIndexArray; delete xWrappedMask; }

private:
  // The actual mask that does all the work
  AffineAndPCACoefficientMask *xWrappedMask;

  // The selection index array
  size_t *xIndexArray;
  
  // Method that generates the index array
  size_t *FillIndexArray(size_t nModes);
};


#endif //__CoefficientMask_h_

