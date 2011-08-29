#ifndef __SubdivisionMedialModel_h_
#define __SubdivisionMedialModel_h_

#include "GenericMedialModel.h"
#include "MeshMedialPDESolver.h"
#include "SubdivisionSurface.h"

/**
 * This is a generic interface for medial models based on subdivision surfaces
 * It supports models that use a PDE and non-PDE models
 */
class SubdivisionMedialModel : public GenericMedialModel
{
public:
  typedef SubdivisionSurface::MeshLevel MeshLevel;

  SubdivisionMedialModel();
  virtual ~SubdivisionMedialModel();

  /**
   * This method associates the medial model with a mesh. The two parameters
   * that follow are levels of subdivision: the number of levels at which the
   * atoms are interpolated and a smaller number, at which the coefficients
   * are generated. By default, the nodes in the mesh itself serve as the
   * coefficients, but it is possible to request that the input mesh first be
   * subdivided and the nodes of that mesh be used as coefficients
   */
  virtual void SetMesh(
    const MeshLevel &mesh, 
    const Vec &C, const Vec &u, const Vec &v,
    size_t nAtomSubs, size_t nCoeffSubs);

  /** 
   * Returns the number of coefficients per mesh control point. It's usually
   * 4 (x,y,z,r) or 5 (x,y,z,rho,tau)
   */
  virtual size_t GetNumberOfComponents() const = 0;

  // Get a list of scaling coefficients. These are the coefficients that should
  // be scaled by similarity transforms
  virtual std::list<size_t> GetScalingComponents() const = 0;

  /** Get the number of optimizable coefficients that define this model */
  size_t GetNumberOfCoefficients() const
    { return xCoefficients.size(); }

  /** Get the array of (optimizable) coefficients that define this model */
  const Vec GetCoefficientArray() const
    { return xCoefficients; }

  /** Get the u array of the coefficients */
  const Vec GetCoefficientU() const { return uCoeff; }
  const Vec GetCoefficientV() const { return vCoeff; }

  /** Get one of the coefficients defining this model */
  double GetCoefficient(size_t i) const
    { return xCoefficients[i]; }

  /** Set one of the coefficients defining this model */
  void SetCoefficient(size_t i, double x)
    { xCoefficients[i] = x; }

  /** Set the array of coefficients defining this model */
  void SetCoefficientArray(const Vec &xData)
    { xCoefficients = xData; }

  /**
   * Load a medial model from the Registry.
   */
  virtual void ReadFromRegistry(Registry &folder);

  /**
   * Save the model to registry folder
   */
  virtual void WriteToRegistry(Registry &folder);

  /**
   * Get the subdivision level
   */
  size_t GetSubdivisionLevel() const { return xSubdivisionLevel; }

  /**
   * Get a pointer to the affine transform descriptor corresponding to this class.
   * The descriptor is a lightweight object, and its allocation should be
   * managed by the child of GenericMedialModel.
   */
  const AffineTransformDescriptor *GetAffineTransformDescriptor() const
    { return xAffineDescriptor; }

  /**
   * Get a pointer to the coarse-to-fine masking descriptor corresponding to
   * this class. This descriptor is managed internally by the child of
   * GenericMedialModel.
   */
  const CoarseToFineMappingDescriptor *
    GetCoarseToFineMappingDescriptor() const;

  /** Get a pointer to the coefficient mesh stored in this model */
  const MeshLevel *GetCoefficientMesh() const
    { return &mlCoefficient; }

  /** Get a pointer to the atom mesh stored in this model */
  const MeshLevel *GetAtomMesh() const
    { return &mlAtom; }

  /** Get the vector of phi-values for this model (equal to xAtoms.F) */
  Vec GetPhi() const 
    { 
    Vec phi(this->GetNumberOfAtoms());
    for(size_t i = 0; i < phi.size(); i++)
      phi[i] = xAtoms[i].F;
    return phi;
    }

  /** Set the vector of phi-values for the model */
  void SetPhi(const Vec &phi)
    {
    for(size_t i = 0; i < phi.size(); i++)
      xAtoms[i].F = phi[i];
    }

protected:
  // Vector typedef
  typedef vnl_vector<double> Vec;

  // Mesh levels corresponding to the coefficients and the atoms
  MeshLevel mlCoefficient, mlAtom;

  // The array of coefficients (4-tuples at each vertex, it is assumed at this
  // point that the coefficients are always xyz + something, whether rho or rad
  // but that may change at some point in time
  Vec xCoefficients;

  // U and V values at the coefficients and atoms
  Vec uCoeff, vCoeff, uAtom, vAtom;

  // Coarse-to-fine mapping descriptor (dummy)
  SubdivisionSurfaceCoarseToFineMappingDescriptor *xCTFDescriptor;

  // Affine transform mapping descriptor
  PointArrayAffineTransformDescriptor *xAffineDescriptor;

  // The subdivision level (number of subdivisions from mlCoefficient to
  // mlAtom)
  size_t xSubdivisionLevel;

  friend class SubdivisionMedialModelIO;
};


#endif // __SubdivisionMedialModel_h_
