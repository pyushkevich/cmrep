#ifndef _CartesianMedialModel_h_
#define _CartesianMedialModel_h_

#include <iostream>
#include <smlmath.h>
#include "CodeTimer.h"
#include "MedialPDEMasks.h"
#include "MedialPDESites.h"
#include "MedialAtom.h"
#include "GenericMedialModel.h"
#include "CartesianGridMedialIterationContext.h"
#include "BasisFunctions2D.h"
#include "PardisoInterface.h"

using namespace std;

class CartesianMedialModel : public GenericMedialModel
{
public:
  // Vector and matrix typedefs
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;
  
  /** Use the parameter grid to initialize the medial PDE solver. */
  CartesianMedialModel(const Vec &uGrid, const Vec &vGrid);

  /**
   * This helper constructor generates a medial PDE grid with given numbers
   * (nu, nv) of regularly spaced grid points. In addition, there is an option
   * to generate irregularly spaced grid points at the boundary of the domain.
   * These grid points have exponentially decreasing spacing as you approach
   * the boundary and are used to make sure that the m-rep boundary is sampled
   * densely near the endcaps. The exponent is currenly 2, but that may change
   */
  CartesianMedialModel(size_t nu, size_t nv, 
    double xScale = 0.0, size_t pu = 0, size_t pv = 0);

  /** Destructor */
  ~CartesianMedialModel();

  /** Get the number of optimizable coefficients that define this model */
  size_t GetNumberOfCoefficients() const
    { return xSurface->GetNumberOfCoefficients(); }

  /** Get the array of (optimizable) coefficients that define this model */
  const Vec GetCoefficientArray() const 
    { return Vec(xSurface->GetCoefficientArray(), xSurface->GetNumberOfCoefficients()); }

  /** Get one of the coefficients defining this model */
  double GetCoefficient(size_t i) const
    { return xSurface->GetCoefficient(i); }

  /** Set one of the coefficients defining this model */
  void SetCoefficient(size_t i, double x) 
    { xSurface->SetCoefficient(i, x); }

  /** Set the array of coefficients defining this model */
  void SetCoefficientArray(const Vec &xData)
    { xSurface->SetCoefficientArray(xData.data_block()); }

  /** Get the hint array based on the current state of the model */
  Vec GetHintArray() const;

  /** Compute the PDE solution and the atoms */
  void ComputeAtoms(const double *xHint = NULL)
    { this->Solve(xHint); }

  /** 
   * Store the current solution vector as the initialization vector
   */
  void SetSolutionAsInitialGuess();
    
  /** 
   * Generate the default initial guess
   */
  void SetDefaultInitialGuess(double xMagnitude = 0.01);

  /** 
   * Specify the surface for which the problem should be solved. The medial
   * surface will be owned by this class, and this class will be responsible
   * for deleting it
   */
  void AdoptMedialSurface(IBasisRepresentation2D *xSurface);

  /** Get the pointer to the internal surface */
  IBasisRepresentation2D *GetMedialSurface() 
    { return xSurface; }

  /**
   * Solve the PDE for a given surface and a given rho function up to the required 
   * level of accuracy.
   */
  void Solve(const double *xHint = NULL, double delta = 1e-12);

  /** 
   * Specify the set of directions (variations) for repeated gradient computations
   * Each row in xBasis specifies a variation.
   */
  void SetVariationalBasis(const Mat &xBasis);

  /** 
   * Method called before multiple calls to ComputeAtomVariationalDerivative()
   */
  void BeginGradientComputation();

  /** 
   * This method computes the variational derivative of the model with respect to 
   * a variation. Index iVar points into the variations passed in to the method
   * SetVariationalBasis();
   */
  void ComputeAtomVariationalDerivative(size_t iVar, MedialAtom *dAtoms);

  /** 
   * Method called before multiple calls to ComputeAtomVariationalDerivative()
   */
  void EndGradientComputation() {}

  /** Get the number of atoms in u dimension */
  unsigned int GetNumberOfUPoints()
    { return m; }

  /** Get the number of atoms in u dimension */
  unsigned int GetNumberOfVPoints()
    { return n; }

  /** Get the radius field */
  const Mat &GetPhiField() { return y; }

  /** Get the phi field for the last variational derivative computation */
  const Mat &GetPhiDerivativeField() { return dy; }

  /** Get the u and v finite difference grids */
  const Vec &GetGridU() const { return uGrid; }
  const Vec &GetGridV() const { return vGrid; }

  /** 
   * Get the index of the atom at position iu, iv in array returned by
   * GetAtomArray() 
   */
  size_t GetGridAtomIndex(size_t iu, size_t iv) const
    { 
    return static_cast<CartesianGridMedialIterationContext *>(this->xIterationContext)
      ->GetAtomIndex(iu,iv);
    }

  /**
   * Load a medial model from the Registry.
   */
  void ReadFromRegistry(Registry &folder);

  /**
   * Save the model to registry folder
   */
  void WriteToRegistry(Registry &folder);

  /** Get the affine transform descriptor for this model */
  const AffineTransformDescriptor *GetAffineTransformDescriptor() const
    { return xSurface->GetAffineTransformDescriptor(); }

  /** Get the C2F descriptor */
  const CoarseToFineMappingDescriptor * GetCoarseToFineMappingDescriptor() const
    { return xSurface->GetCoarseToFineMappingDescriptor(); }

  void TestFiniteDifferenceConvergence();

  CodeTimer tSolver;
private:

  // Numbers of grid points
  size_t m, n;
  
  // U and V grids
  Vec uGrid, vGrid;

  /** Number of 'sites' */
  unsigned int nSites;

  /** The hypersurface on which the PDE is solved */
  IBasisRepresentation2D *xSurface;

  // Array of mask pointers
  std::vector<FiniteDifferenceMask *> xMasks;

  /** Array of sites */
  std::vector<FDAbstractSite *> xSites;

  /** Index into the sites */
  vnl_matrix<size_t> xSiteIndex;

  /** Representation for the sparse matrix */
  double *xSparseValues;
  int *xRowIndex, *xColIndex, nSparseEntries;

  /** Three vectors used in the iteration */
  Mat eps, b, y, zTest, xInitSoln, xDefaultInitSoln, dy;

  /** Representation of the variational basis */
  struct VariationalBasisAtomData {
    SMLVec3d X, Xu, Xv, Xuu, Xuv, Xvv;
    double xLapR;
    VariationalBasisAtomData() 
      : X(0.0), Xu(0.0), Xv(0.0), Xuu(0.0), Xuv(0.0), Xvv(0.0), xLapR(0.0) {}
  };

  // Matrix storing all data associated with a variational basis
  typedef std::vector<VariationalBasisAtomData> VariationRep;
  typedef std::vector<VariationRep> VariationalBasisRep;
  VariationalBasisRep xVariationalBasis;

  // Array that holds temporary data for gradient computation
  std::vector<MedialAtom::DerivativeTerms> xTempDerivativeTerms;

  // Sparse linear matrix solver
  UnsymmetricRealPARDISO xPardiso;

  // Common initialization
  void Initialize(const Vec &uGrid, const Vec &vGrid);

  /** Routine to compute medial atoms */
  void TestJacobi();
  void ReconstructAtoms(const Mat &ySolution);
  void InitializeSiteGeometry();
  void SolveOnce(const double *xHint, double delta);
  double EstimateLBOperator(const Mat &F, size_t i, size_t j);
  double ComputeNewtonRHS(const Mat& x, Mat &b);

  bool flagReuseLastSolution;

};

/* *********************** Template Code ************************** */


#endif // _CartesianMedialModel_h_
