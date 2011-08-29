#include "SubdivisionMedialModel.h"
#include "SubdivisionSurfaceMedialIterationContext.h"

SubdivisionMedialModel::SubdivisionMedialModel() 
{
  xCTFDescriptor = NULL;
  xAffineDescriptor = NULL;

  xAtoms = NULL;
  xIterationContext = NULL;
  xSubdivisionLevel = 0;
}

SubdivisionMedialModel::~SubdivisionMedialModel()
{
  if(xCTFDescriptor)
    delete xCTFDescriptor;
  if(xAffineDescriptor)
    delete xAffineDescriptor;
}

void
SubdivisionMedialModel
::SetMesh(
  const MeshLevel &mesh, 
  const Vec &C, const Vec &u, const Vec &v,
  size_t nAtomSubs, size_t nCoeffSubs)
{
  // Validity check
  assert(nAtomSubs >= nCoeffSubs);

  // Set the subdivision level
  xSubdivisionLevel = nAtomSubs - nCoeffSubs;

  // Create a vector of mesh levels
  vector<const MeshLevel *> xTempLevels; xTempLevels.push_back(&mesh);

  // Subdivide the input mesh into the coefficient-level mesh
  SubdivisionSurface::RecursiveSubdivide(&mesh, &mlCoefficient, nCoeffSubs);

  // Get number of components per control point
  size_t nc = GetNumberOfComponents();

  // Compute the coefficients and u/v arrays from the input data
  if(nCoeffSubs > 0)
    {
    // Initialize the arrays
    xCoefficients.set_size(mlCoefficient.nVertices * nc);
    uCoeff.set_size(mlCoefficient.nVertices);
    vCoeff.set_size(mlCoefficient.nVertices);

    // Apply the subdivision to the coefficients
    SubdivisionSurface::ApplySubdivision(
      C.data_block(), xCoefficients.data_block(), nc, mlCoefficient);

    // Apply to the u and v arrays
    SubdivisionSurface::ApplySubdivision(
      u.data_block(), uCoeff.data_block(), 1, mlCoefficient);
    SubdivisionSurface::ApplySubdivision(
      v.data_block(), vCoeff.data_block(), 1, mlCoefficient);
    }
  else
    {
    xCoefficients = C;
    uCoeff = u; vCoeff = v;
    }

  // Set the coefficient-level mesh as the new root (forget the input)
  mlCoefficient.SetAsRoot();

  // Subdivide the coefficient-level mesh up to the atom-level mesh
  SubdivisionSurface::RecursiveSubdivide(
    &mlCoefficient, &mlAtom, nAtomSubs - nCoeffSubs);

  // Apply the subdivision to the u and v coordinates
  uAtom.set_size(mlAtom.nVertices); vAtom.set_size(mlAtom.nVertices);
  SubdivisionSurface::ApplySubdivision(
    uCoeff.data_block(), uAtom.data_block(), 1, mlAtom);
  SubdivisionSurface::ApplySubdivision(
    vCoeff.data_block(), vAtom.data_block(), 1, mlAtom);

  // Create the atoms array
  if(xAtoms) delete xAtoms;
  xAtoms = new MedialAtom[mlAtom.nVertices];

  // Copy the u, v values into the atoms
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    xAtoms[i].u = uAtom[i];
    xAtoms[i].v = vAtom[i];
    }

  // Set up the iteration context
  if(this->xIterationContext != NULL) delete this->xIterationContext;
  this->xIterationContext = new SubdivisionSurfaceMedialIterationContext(&mlAtom);

  // Set the coarse-to-fine descriptor
  if(xCTFDescriptor)
    delete xCTFDescriptor;
  xCTFDescriptor = new SubdivisionSurfaceCoarseToFineMappingDescriptor(mlCoefficient.nVertices);


  if(xAffineDescriptor)
    delete xAffineDescriptor;
  xAffineDescriptor = 
    new PointArrayAffineTransformDescriptor(
      GetNumberOfComponents(), GetScalingComponents());

}

void
SubdivisionMedialModel::
WriteToRegistry(Registry &R)
{
  // Mesh type
  R["Grid.Type"] << "LoopSubdivision";

  // Save the subdivision level info
  R["Grid.Model.Atom.SubdivisionLevel"] << GetSubdivisionLevel();

}

void
SubdivisionMedialModel::
ReadFromRegistry(Registry &R)
{
}

const CoarseToFineMappingDescriptor *
SubdivisionMedialModel
::GetCoarseToFineMappingDescriptor() const
{
  return xCTFDescriptor;

}
