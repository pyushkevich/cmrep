#ifndef __BasisFunctions2D_h_
#define __BasisFunctions2D_h_

#include <valarray>
#include <smlmath.h>
#include "Registry.h"
#include "AffineTransformDescriptor.h"
#include "CoarseToFineMappingDescriptor.h"

using namespace std;

class IMedialCoefficientMask;

/*
 * We need to design a basic interface that can encapsulate all types of inputs
 * to medial PDE solvers. These inputs can be subdivision surfaces, Fourier
 * surfaces or they eventually could be stratified sets. So we need to design a
 * minimal set of functionality that all of these inputs support.
 *
 * Each input maps a number of coefficients (C1...CN) to an array of points in
 * space (X1...XM). It is safe to assume that the neighborhood relationship
 * between these points will stay fixed regardless of the coefficient values.
 *
 * This abstract class should not encapsulate what its output should be. It can
 * be a grid of points with a local frame, or it can be a mesh of points without
 * any geometrical information. It can be rho-values or R-values. The design
 * should be as flexible as possible.
 *
 * What the class does need to support is the ability to apply affine
 * transformations. Basically, there needs to be some logic as to how the
 * coefficients C1...CN must change in order to apply an affine transformation
 * to the output.
 */

class IHyperSurface2D
{
public:
  virtual ~IHyperSurface2D() {}

  // Typed definitions
  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  virtual void Evaluate(double u, double v, double *x) = 0;
 
  /** Evaluate a partial derivative of the function at the u, v coordinate.
   * The first component is c0, number of components nc */
  virtual void EvaluateDerivative(
    double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x) = 0;

  /** Evaluate the function at a grid point */
  virtual void EvaluateAtGridIndex(
    size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x) = 0;

  /** Get the evaluation grid parameters */
  virtual void GetEvaluationGrid(VectorType &uu, VectorType &vv) const = 0;

  /** Get the number of dimensions of the surface */
  virtual size_t GetNumberOfDimensions() = 0;

  /** Print a report about this surface */
  virtual void PrintReport() {};
};

class IMutableHyperSurface2D : virtual public IHyperSurface2D
{
public:
  virtual ~IMutableHyperSurface2D () {}

  /** Specify a grid of values along which the function will need to be
   * evaluated repeatedly. */
  virtual void SetEvaluationGrid(const VectorType &uu, const VectorType &vv) = 0;

  /** Get the affine transform descriptor for this surface */
  virtual const AffineTransformDescriptor *GetAffineTransformDescriptor() const = 0;

  /** Get the coarse to fine descriptor for this type of surface */
  virtual const CoarseToFineMappingDescriptor * 
    GetCoarseToFineMappingDescriptor() const = 0;
};

/**
 * This is an interface that contains functions for setting and getting
 * coefficients and corresponding component surfaces
 */
class ICoefficientSettable
{
public:
  virtual ~ICoefficientSettable() {}

  virtual double *GetCoefficientArray() = 0;
  virtual void SetCoefficientArray(const double *xData) = 0;
  virtual size_t GetNumberOfCoefficients() const = 0 ;
  virtual double GetCoefficient(size_t i) const = 0;
  virtual void SetCoefficient(size_t i, double x) = 0;

  // Get a surface corresponding to a single component
  virtual IHyperSurface2D *GetComponentSurface(size_t iCoefficient) = 0;
  virtual void ReleaseComponentSurface(IHyperSurface2D *xSurface) = 0;

  // Get a variation surface corresponding to a given change in coeffs
  virtual IHyperSurface2D *GetVariationSurface(const double *xCoeff) = 0;
  virtual void ReleaseVariationSurface(IHyperSurface2D *xSurface) = 0;

  typedef vnl_vector<double> Vec;
};

/** This is an interface that should be implemented by any 2D basis function
 * implementation */
class IBasisRepresentation2D : 
  virtual public IMutableHyperSurface2D, 
  virtual public ICoefficientSettable
{
public:
  virtual ~IBasisRepresentation2D () {}


  /** Fit the i-th component of the coefficients to some data points */
  virtual void FitToData(size_t n, size_t i, double *uu, double *vv, double *xx) = 0;

  /** Save the coefficients to a registry */
  virtual void SaveToRegistry(Registry &R) = 0;

  /** Read the coefficients from a registry */
  virtual bool ReadFromRegistry(Registry &R) = 0;
};
  
/**
 * A 3D index into an array. As the array is traversed, the x increases
 * fastest, z slowest
 */
class Index3D
{
public:
  Index3D(size_t nx, size_t ny, size_t nz) 
    { 
    stride_y = nx; 
    stride_z = nx * ny; 
    sz = nx * ny * nz;
    data = new double[sz];
    }

  ~Index3D()
    {
    delete data;
    }

  void resize(size_t nx, size_t ny, size_t nz)
    {
    stride_y = nx; 
    stride_z = nx * ny;
    sz = nx * ny * nz;
    delete data;
    data = new double[sz];    
    }

  size_t size() const
    { return sz; }

  double &operator[] (size_t idx)
    { return data[idx]; }
  
  const double &operator[] (size_t idx) const
    { return data[idx]; }
  
  double &operator() (size_t ix, size_t iy, size_t iz) 
    { 
    size_t index = iz * stride_z + iy * stride_y + ix;
    return data[index];
    }

  const double &operator() (size_t ix, size_t iy, size_t iz) const
    { 
    size_t index = iz * stride_z + iy * stride_y + ix;
    return data[index];
    }

  double *GetPointer()
    { return data; }

  const double *GetPointer() const
    { return data; }

private:
  size_t stride_z, stride_y, sz;
  double *data;
};

/*
class Index3D : public valarray<double> 
{
public:
  Index3D(size_t nx, size_t ny, size_t nz) 
    : valarray<double> (nx * ny * nz)
    { stride_y = nx; stride_z = nx * ny; }

  void resize(size_t nx, size_t ny, size_t nz)
    {
    stride_y = nx; 
    stride_z = nx * ny; 
    valarray<double>::resize(nx * ny * nz);
    (*this) *= 0.0; 
    }
    
  double &operator() (size_t ix, size_t iy, size_t iz) 
    { 
    size_t index = iz * stride_z + iy * stride_y + ix;
    double &ref = (*this)[index];
    return ref;
    }

  const double &operator() (size_t ix, size_t iy, size_t iz) const
    { 
    size_t index = iz * stride_z + iy * stride_y + ix;
    double &ref = (*this)[index];
    return ref;
    }

  double *GetPointer()
    { return &((*this)[0]); }

  const double *GetPointer() const
    { return &((*this)[0]); }

private:
  size_t stride_z, stride_y;
};
*/

/** 
 * A generic implementation that uses some arbitrary basis function class.
 * This representation is modeled after Fourier basis functions, but it could
 * be arbitrary. The basic idea is that the 2D basis functions are separable
 * into 1D components, i.e., F_k(u,v) = f_i(u) f_j(v)
 */
template< size_t NComponents, size_t NOrder, 
  typename BasisFunctionU, typename BasisFunctionV >
class GenericBasisRepresentation2D : virtual public IBasisRepresentation2D
{
public:
  virtual ~GenericBasisRepresentation2D () {}

  // Some useful typedefs
  typedef GenericBasisRepresentation2D<
    NComponents,NOrder,BasisFunctionU,BasisFunctionV> Self;

  /** Class that allows you to return a single basis function as a surface;
   * used in variational calculus computations */
  class SingleFunctionAdapter : public IHyperSurface2D
  {
  public:
    void Evaluate(double u, double v, double *x);
    void EvaluateDerivative(
      double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x);

    void GetEvaluationGrid(VectorType &uu, VectorType &vv) const;
    void EvaluateAtGridIndex(
      size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x);
    size_t GetNumberOfDimensions();

  private:
    GenericBasisRepresentation2D *xParent;
    size_t icu, icv, iComp;

    friend class GenericBasisRepresentation2D;
  };
      
  /** Construct the representation with a specific number of coefficients in u
   * and v components */
  GenericBasisRepresentation2D(size_t ncu, size_t ncv);

  void SetNumberOfCoefficients(size_t ncu, size_t ncv);
  
  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  void Evaluate(double u, double v, double *x)
    { return EvaluateDerivative(u, v, 0, 0, 0, NComponents, x); } 
 
  /** Evaluate a partial derivative of the function at the u, v coordinate */
  void EvaluateDerivative(double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x); 

  /** Specify a grid of values along which the function will need to be
   * evaluated repeatedly. Call PrecomputeGrid() before evaluating at grid and
   * after changing any of the coefficients */
  void SetEvaluationGrid(const VectorType &uu, const VectorType &vv);

  /** Get the vectors that define the evaluation grid */
  void GetEvaluationGrid(VectorType &uu, VectorType &vv) const;

  /** Evaluate the function at a grid point */
  void EvaluateAtGridIndex(
    size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x); 

  /** Get the number of components/dimensions */
  size_t GetNumberOfDimensions() 
    { return NComponents; }

  /** Get number of raw coefficients */
  size_t GetNumberOfCoefficients() const 
    { return C.size(); }

  /** Get the coefficients as an array of real values */
  double *GetCoefficientArray() 
    { return C.GetPointer(); }
  
  /** Get the coefficients as an array of real values */
  const double *GetCoefficientArray() const
    { return C.GetPointer(); }

  /** Set the array of coefficients */
  void SetCoefficientArray(const double *data)
    { for(size_t i = 0; i < C.size(); i++) C[i] = data[i]; }

  // Direct coefficient access
  double GetCoefficient(size_t iCoeff) const
    { return C[iCoeff]; }
  
  void SetCoefficient(size_t iCoeff, double xValue) 
    { C[iCoeff] = xValue; }

  void SetAllCoefficients(double xValue)
    {  
    for(size_t i=0; i < C.size(); i++)
      C[i] = xValue;
    }

  // Set a coefficient by index
  void SetCoefficient(size_t iu, size_t iv, size_t iComp, double value)
    { C(iComp,iu,iv) = value; }

  double GetCoefficient(size_t iu, size_t iv, size_t iComp) const
    { return C(iComp,iu,iv); }

  /** Fit the i-th component of the coefficients to some data points */
  void FitToData(size_t n, size_t i, double *uu, double *vv, double *xx); 

  /** Save the coefficients to a registry */
  void SaveToRegistry(Registry &R);

  /** Read the coefficients from a registry */
  bool ReadFromRegistry(Registry &R);

  /** Get the number of coefficients in U and V */
  void GetNumberOfCoefficientsUV(size_t &ncu, size_t &ncv)
    { ncu = this->ncu; ncv = this->ncv; }

  /** Debug report */
  void PrintReport();

  /** Get one of the components as a surface */
  IHyperSurface2D *GetComponentSurface(size_t iCoefficient);

  /** Free the memory associated with the component surface */
  void ReleaseComponentSurface(IHyperSurface2D *xSurface);

  /** Get the component by u, v, comp index */
  IHyperSurface2D *GetComponentSurface(size_t icu, size_t icv, size_t iComp);

protected:
  // The raw array of coefficients (without pointers)
  Index3D C;

  // The evaluation grids along u and v. 
  valarray<double> uGrid, vGrid, uGridValues, vGridValues, uEval, vEval;

  // Number of coefficients, raw coefficients
  size_t ncu, ncv, nc, ncRaw;

  // The basis functions
  BasisFunctionU fu; 
  BasisFunctionV fv;

  // Get the slice corresponding to coefficient iu, iv
  slice GetCoefficientIndex(size_t iu, size_t iv)
    { return slice(((ncu * iv) + iu) * NComponents, NComponents, 1); }

  size_t GetCoefficientIndex(size_t iu, size_t iv, size_t k)
    { return ((ncu * iv) + iu) * NComponents + k; }

  void GetRawCoefficientIndices(size_t iRaw, size_t &icu, size_t &icv, size_t &iComp)
    {
    iComp = iRaw % NComponents;
    icu = (iRaw / NComponents) % ncu;
    icv = iRaw / (NComponents * ncu);
    }

  void Initialize(size_t ncu, size_t ncv);
};

/** The cosine basis function */
class CosineBasisFunction
{
public:
  double Evaluate(double u, size_t k, size_t d);
};


/** Declare the Fourier surface class */
typedef 
GenericBasisRepresentation2D <4, 3, CosineBasisFunction, CosineBasisFunction>
FourierSurfaceBase;

class FourierSurface : public FourierSurfaceBase
{
public:

  // Constructor
  FourierSurface(size_t ncu, size_t ncv);

  // Copy constructor
  FourierSurface(const FourierSurface &source); 

  virtual ~FourierSurface () {}

  /** Get a surface corresponding to a variation from this surface */
  IHyperSurface2D *GetVariationSurface(const double *xCoeff);
  void ReleaseVariationSurface(IHyperSurface2D *xSurface);

  /** Get the affine surface descriptor pointer */
  const AffineTransformDescriptor *GetAffineTransformDescriptor() const
    { return &xAffineDescriptor; }

  /** Get the C2F surface descriptor */
  const CoarseToFineMappingDescriptor *GetCoarseToFineMappingDescriptor() const
    { return &xC2FDescriptor; }

  /**
   * Get the list of raw coefficient indices that correspond to the first 
   * (ncu x ncv) components in X and Rho
   */
  vector<size_t> GetCoefficientSubset(
    size_t ncuX, size_t ncvX, size_t ncuRho, size_t ncvRho);

private:
  // The affine transform descriptor
  FourierAffineTransformDescriptor xAffineDescriptor;

  // The coarse-to-fine mapping descriptor
  FourierCoarseToFineMappingDescriptor xC2FDescriptor;
};

#endif

