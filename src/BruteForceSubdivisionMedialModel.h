#ifndef __BruteForceSubdivisionMedialModel_h_
#define __BruteForceSubdivisionMedialModel_h_

#include "SubdivisionMedialModel.h"
#include "MedialAtomGrid.h"

/**
 * This is a type of subdivision medial model that does not solve any PDE.
 * Instead, in this model the subdivision surface represents the function
 * (xyzr) and no equality contraints are enforced. Who knows if that is going to
 * work, but at least we can try!
 */
class BruteForceSubdivisionMedialModel : public SubdivisionMedialModel
{
public:
  // Vector typedef
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  // Mesh level definition from parent
  typedef SubdivisionMedialModel::MeshLevel MeshLevel;

  BruteForceSubdivisionMedialModel();
  ~BruteForceSubdivisionMedialModel();

  /** There are 4 components (x,y,z,r) */
  size_t GetNumberOfComponents() const { return 4; }

  /** The fifth component is scaling */
  std::list<size_t> GetScalingComponents() const
    { std::list<size_t> lst; lst.push_back(3); return lst; }

  void SetMesh(const MeshLevel &mesh, 
    const Vec &C, const Vec &u, const Vec &v,
    size_t nAtomSubs, size_t nCoeffSubs);

  /** Get the hint array. This returns a single double, it's a dummy method */
  Vec GetHintArray() const;

  /** Compute the atoms from given a set of coefficients.  */
  void ComputeAtoms(bool flagAllowErrors, const double *xHint = NULL);


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
  void EndGradientComputation();

  /**
   * This method is called before multiple calls to the ComputeJet routine. Given a
   * of variation (direction in which directional derivatives will be computed), and
   * an array of atoms corresponding to the variation, this method will precompute the
   * parts of the atoms that do not depend on the position in coefficient space where
   * the gradient is taken. This method must be called for each variation in the Jet.
   */
  // void PrepareAtomsForVariationalDerivative(
  //  const Vec &xVariation, MedialAtom *dAtoms) const;

  /**
   * Gradient computation
   */
  // void ComputeAtomVariationalDerivative(const Vec &xVariation, void *xNonVarying);

  /**
   * Load a medial model from the Registry.
   */
  void ReadFromRegistry(Registry &folder);

  /**
   * Save the model to registry folder
   */
  void WriteToRegistry(Registry &folder);

private:

  // A gradient calculator
  MeshGradientComputer xGradComp;

  // Loop scheme for computing tangents
  LoopTangentScheme xLoopScheme;

  // Sparse matrices containing the contribuion of each mesh node
  // to partial derivatives in U and V
  typedef std::pair<double, double> UVPair;
  typedef ImmutableSparseArray<UVPair> SparseMat;
  SparseMat Wuv;

  // Struct representing nonvarying terms in variational derivatives
  struct NonvaryingAtomTerms {
    SMLVec3d X, Xu, Xv;
    double R, Ru, Rv;
    int order;
    NonvaryingAtomTerms() : 
      X(0.0), Xu(0.0), Xv(0.0), R(0.0), Ru(0.0), Rv(0.0), order(3) {};
  };


  typedef ImmutableSparseArray<NonvaryingAtomTerms> NonvaryingTermsMatrix;

  // The sparse matrix representing the basis for gradient computation
  NonvaryingTermsMatrix xBasis;

  // Array of common derivative terms
  std::vector<MedialAtom::DerivativeTerms> dt;

  // Coefficients used in computing the derivatives of Ru
  struct LocalDerivativeTerms 
    { 
    double w_Rv, w_g12, w_g22, w_g, Ru_orig, Ru_fixed; 
    double sqrt_gz;
    };
    
  std::vector<LocalDerivativeTerms> dt_local;
};


#endif // __BruteForceSubdivisionMedialModel_h_
