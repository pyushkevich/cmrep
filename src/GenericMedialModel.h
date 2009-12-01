#ifndef __GenericMedialModel_h_
#define __GenericMedialModel_h_

#include "MedialAtom.h"
#include "MedialAtomIterators.h"
#include "MedialException.h"
#include "CoarseToFineMappingDescriptor.h"
#include "AffineTransformDescriptor.h"
#include "Registry.h"
#include <iostream>
#include <vector>

class AffineTransformDescriptor;

/**
 * This class represents a generic medial model, whether Cartesian-based or
 * subdivision surface based. The user can inquire basic information about the
 * model, such as the number of atoms, etc.
 *
 * This is how gradient computation works. A gradient, in a general sense,
 * is a set of directional (i.e. variational) derivatives with respect to 
 * a given set of directions (variations). The set of directions may simply
 * be the coordinate axes, but in this code, any set of directions is 
 * admitted. We refer to this set of directions as "Variational Basis"
 *
 * Certain terms in the gradient computation depend only on the variational 
 * basis itself (i.e. the directions themselves) and not on the state of the
 * model at the time the gradient is evaluated. These terms are linear with
 * respect to the directions. For instance, the X, Xu, Xv, etc. on the medial
 * surface are typically such terms. 
 *
 * To make the computation more efficient, the code requires the children of the 
 * GenericMedialModel class to precompute these terms once before a set of 
 * gradient computations (since the gradient is usually evaluated 100s of times
 * in an optimization problem). This is done by the children in the method
 * SetVariationalBasis(). The user specifies the set of directions with respect
 * to which derivatives will be taken, and the method computes whatever terms
 * are going to pepend on these directions and not on the state of the model. 
 *
 * The actual computation of the gradient consists of three steps. First the 
 * method BeginGradientComputation() is called. Then, for each direction, the
 * method ComputeAtomVariationalDerivative() is called. Finally, the method
 * EndGradientComputation() is executed. The Begin and End methods are intented 
 * to allow the MedialModel to precompute a set of terms that only depend on the
 * state of the model and not on the variation. Again, this is done to minimize
 * duplication of work by the program. 
 *
 * Note that the call to ComputeAtomVariationalDerivative() does not take as a 
 * parameter the actual variation. It is the responsibility of the MedialModel
 * to store the variations - or necessary terms derived from the variations - 
 * during the call to SetVariationalBasis
 */
class GenericMedialModel
{
public:
  // Vector typedef
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;
  
  size_t GetNumberOfAtoms() const
    { return xIterationContext->GetNumberOfAtoms(); }

  size_t GetNumberOfBoundaryPoints() const
    { return xIterationContext->GetNumberOfBoundaryPoints(); }

  size_t GetNumberOfInternalPoints(size_t ncuts) const
    { return xIterationContext->GetNumberOfInternalPoints(ncuts); }

  size_t GetNumberOfTriangles() const
    { return xIterationContext->GetNumberOfTriangles(); }

  size_t GetNumberOfBoundaryTriangles() const
    { return xIterationContext->GetNumberOfBoundaryTriangles(); }

  size_t GetNumberOfCells(size_t ncuts) const
    { return xIterationContext->GetNumberOfCells(ncuts); }

  size_t GetNumberOfProfileIntervals(size_t ncuts) const
    { return xIterationContext->GetNumberOfProfileIntervals(ncuts); }
  
  MedialAtomIterator GetAtomIterator() const
    { return MedialAtomIterator(xIterationContext); }

  MedialTriangleIterator GetMedialTriangleIterator() const
    { return MedialTriangleIterator(xIterationContext); }

  MedialBoundaryPointIterator GetBoundaryPointIterator() const
    { return MedialBoundaryPointIterator(xIterationContext); }

  MedialBoundaryTriangleIterator GetBoundaryTriangleIterator() const
    { return MedialBoundaryTriangleIterator(xIterationContext); }

  MedialInternalPointIterator GetInternalPointIterator(size_t ncuts) const
    { return MedialInternalPointIterator(xIterationContext, ncuts); }

  MedialInternalCellIterator GetInternalCellIterator(size_t ncuts) const
    { return MedialInternalCellIterator(xIterationContext, ncuts); }

  MedialProfileIntervalIterator GetProfileIntervalIterator(size_t ncuts) const
    { return MedialProfileIntervalIterator(xIterationContext, ncuts); }

  /** Get the array of medial atoms - child implements */
  MedialAtom *GetAtomArray() const 
    { return xAtoms; }

  /** Get the iteration context */
  MedialIterationContext *GetIterationContext() const
    { return xIterationContext; }

  /** Get the boundary point corresponding to the current position of a
   * boundary point iterator */
  BoundaryAtom &GetBoundaryPoint(const MedialBoundaryPointIterator &it)
    { return xAtoms[it.GetAtomIndex()].xBnd[it.GetBoundarySide()]; }

  /** Get the number of optimizable coefficients that define this model */
  virtual size_t GetNumberOfCoefficients() const = 0;

  /** Get the array of (optimizable) coefficients that define this model */
  virtual const Vec GetCoefficientArray() const = 0;

  /** Get one of the coefficients defining this model */
  virtual double GetCoefficient(size_t i) const = 0;

  /** Set one of the coefficients defining this model */
  virtual void SetCoefficient(size_t i, double x) = 0;

  /** Set the array of coefficients defining this model */
  virtual void SetCoefficientArray(const Vec &xData) = 0;

  /** 
   * Get 'hint' for computing atoms near the current solution. The 'hint' is
   * a vector specific to the individual model, and for some models there may
   * not be any hints. The hints are relevant for PDE-based models, where the
   * hint vector is the solution of the PDE. By initializing the solver with
   * this solution, solutions at nearby points can be computed easily
   */
  virtual Vec GetHintArray() const = 0;

  /** 
   * Compute the atoms from given a set of coefficients. In PDE-based models
   * this step involves solving the PDE. The optional parameter is the hint
   * array that may increase the efficiency of the computation
   */
  virtual void ComputeAtoms(bool flagAllowErrors, const double *xHint = NULL) = 0;

  /** 
   * Specify the set of directions (variations) for repeated gradient computations
   * Each row in xBasis specifies a variation.
   */
  virtual void SetVariationalBasis(const Mat &xBasis) = 0;

  /** 
   * Method called before multiple calls to ComputeAtomVariationalDerivative()
   */
  virtual void BeginGradientComputation() = 0;

  /** 
   * This method computes the variational derivative of the model with respect to 
   * a variation. Index iVar points into the variations passed in to the method
   * SetVariationalBasis();
   */
  virtual void ComputeAtomVariationalDerivative(size_t iVar, MedialAtom *dAtoms) = 0;

  /** 
   * Method called before multiple calls to ComputeAtomVariationalDerivative()
   */
  virtual void EndGradientComputation() = 0;

  /**
   * Get a pointer to the affine transform descriptor corresponding to this class. 
   * The descriptor is a lightweight object, and its allocation should be
   * managed by the child of GenericMedialModel.
   */
  virtual const AffineTransformDescriptor *GetAffineTransformDescriptor() const = 0;

  /**
   * Get a pointer to the coarse-to-fine masking descriptor corresponding to
   * this class. This descriptor is managed internally by the child of
   * GenericMedialModel.
   */
  virtual const CoarseToFineMappingDescriptor *
    GetCoarseToFineMappingDescriptor() const = 0;

  /**
   * Get the mask that selects radial coefficients in the model. For all models
   * defined so far, this mask selects every fourth component of the coefficient
   * vector, so we define this function at the top level
   */
  virtual vnl_vector<size_t> GetRadialCoefficientMask() 
    {
    size_t n = this->GetNumberOfCoefficients();
    vnl_vector<size_t> mask(n, 0);
    for(size_t i = 3; i < n; i+=4)
      mask[i] = 1;
    return mask;
    }

  /**
   * Get the mask that selects radial coefficients in the model. For all models
   * defined so far, this mask selects every fourth component of the coefficient
   * vector, so we define this function at the top level
   */
  virtual vnl_vector<size_t> GetSpatialCoefficientMask() 
    {
    size_t n = this->GetNumberOfCoefficients();
    vnl_vector<size_t> mask(n, 1);
    for(size_t i = 3; i < n; i+=4)
      mask[i] = 0;
    return mask;
    }

  /** Get the center of rotation for the model */
  virtual SMLVec3d GetCenterOfRotation() 
    { 
    return this->GetAffineTransformDescriptor()->GetCenterOfRotation(
      this->GetCoefficientArray());
    }

  virtual ~GenericMedialModel() {}

protected:

  // Hidden constructor
  GenericMedialModel()
    { xIterationContext = NULL; xAtoms = NULL; }
  

  // The context (data structure) used to facilitate iteration over atoms
  // This data structure must be initialized before iterators can be created
  MedialIterationContext *xIterationContext;

  // Medial atoms array (managed by children)
  MedialAtom *xAtoms;
};

#endif
