#include "PDESubdivisionMedialModel.h"
#include "SubdivisionSurfaceMedialIterationContext.h"

/**
 * Constructor does nothing, just calls parent constructor
 */
PDESubdivisionMedialModel::PDESubdivisionMedialModel() :
  SubdivisionMedialModel()
{
}

void
PDESubdivisionMedialModel
::SetMesh(const MeshLevel &mesh, 
  const Vec &C, const Vec &u, const Vec &v,
  size_t nAtomSubs, size_t nCoeffSubs)
{
  // Call the parent method
  SubdivisionMedialModel::SetMesh(mesh, C, u, v, nAtomSubs, nCoeffSubs);

  // Pass the mesh information to the solver
  xSolver.SetMeshTopology(&mlAtom, xAtoms);
}



void 
PDESubdivisionMedialModel
::ComputeAtoms(bool flagAllowErrors, const double *xHint)
{
  size_t i;

  // The first step is to compute X, tau and rho of the atoms based on the
  // coefficients. This step is performed in this class, not in the solver
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i]; 
    a.X.fill(0.0); 
    a.xLapR = 0.0;
    a.R = 0.0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::RowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() * 5; 
      double w = it.Value();
      a.X += w * xCoefficients.extract(3, j);
      a.xLapR += w * xCoefficients[j + 3];
      a.R += w * xCoefficients[j + 4];
      }

    // If R is negative at an edge vertex, we have a problem, can not recover
    if(a.flagCrest && a.R <= 0.0)
      throw MedialModelException("Non-positive tau passed to PDESubdivisionMedialModel");
    }

  // Now have the solver solve the equation
  xSolver.SolveEquation(true, flagAllowErrors);
}

PDESubdivisionMedialModel::Vec
PDESubdivisionMedialModel::GetHintArray() const
{
  Vec xHint(mlAtom.nVertices, 0.0);
  for(size_t i = 0; i < xHint.size(); i++)
    xHint[i] = xAtoms[i].F;
  return xHint;
}

void
PDESubdivisionMedialModel
::SetVariationalBasis(const Mat &xBasis)
{
  // Allocate the array of terms linearly dependent on the variation
  xVariationalBasis = 
    VariationalBasisRep(xBasis.rows(), VariationRep(mlAtom.nVertices));

  // Loop over the variations
  for(size_t var = 0; var < xBasis.rows(); var++)
    {
    // The current variation
    Vec xVariation = xBasis.get_row(var);

    for(size_t i = 0; i < mlAtom.nVertices; i++)
      {
      // Set up i-th atom
      VariationalBasisAtomData &vbad = xVariationalBasis[var][i];

      // Compute the weighted sum of the coefficients
      ImmutableSparseMatrix<double>::RowIterator it = mlAtom.weights.Row(i);
      for( ; !it.IsAtEnd(); ++it)
        {
        size_t j = it.Column() * 5; 
        double w = it.Value();
        vbad.X += w * xVariation.extract(3, j);
        vbad.xLapR += w * xVariation[j+3];
        vbad.R += w * xVariation[j+4];
        }
      }
    }
}

void
PDESubdivisionMedialModel
::BeginGradientComputation()
{
  xSolver.BeginGradientComputation();
}

void
PDESubdivisionMedialModel
::ComputeAtomVariationalDerivative(size_t ivar, MedialAtom *dAtoms)
{
  // Set whatever we can in the dAtoms array
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    VariationalBasisAtomData &vbad = xVariationalBasis[ivar][i];
    dAtoms[i].X = vbad.X;
    dAtoms[i].xLapR = vbad.xLapR;
    dAtoms[i].R = vbad.R;
    }

  xSolver.ComputeAtomVariationalDerivative(dAtoms);
}


void
PDESubdivisionMedialModel::
WriteToRegistry(Registry &R)
{
  SubdivisionMedialModel::WriteToRegistry(R);

  // Set the model subtype
  R["Grid.Model.SolverType"] << "PDE";
}

void
PDESubdivisionMedialModel::
ReadFromRegistry(Registry &R)
{
}
