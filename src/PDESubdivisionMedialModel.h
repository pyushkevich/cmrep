#ifndef __PDESubdivisionMedialModel_h_
#define __PDESubdivisionMedialModel_h_

#include "SubdivisionMedialModel.h"

class PDESubdivisionMedialModel : public SubdivisionMedialModel
{
public:
  // Mesh level definition from parent
  typedef SubdivisionMedialModel::MeshLevel MeshLevel;

  PDESubdivisionMedialModel();

  /** Set the mesh and initial coefficient array */
  void SetMesh(
    const MeshLevel &mesh, 
    const Vec &C, const Vec &u, const Vec &v,
    size_t nAtomSubs, size_t nCoeffSubs);
  
  /** Get the solver */
  MeshMedialPDESolver *GetSolver() 
    { return &xSolver; }

  /** Get the hint array (solution phi) */
  virtual Vec GetHintArray() const;

  /** Compute the atoms from given a set of coefficients.  */
  void ComputeAtoms(const double *xHint = NULL);

  Vec ComputeLBO(const double *phi);

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

  /**
   * Load a medial model from the Registry.
   */
  void ReadFromRegistry(Registry &folder);

  /**
   * Save the model to registry folder
   */
  void WriteToRegistry(Registry &folder);

private:
  // Vector typedef
  typedef vnl_vector<double> Vec;

  // PDE solver which does most of the work in this method
  MeshMedialPDESolver xSolver;

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

};


#endif // __PDESubdivisionMedialModel_h_
