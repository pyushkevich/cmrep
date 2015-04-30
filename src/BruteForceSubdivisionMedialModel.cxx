#include "BruteForceSubdivisionMedialModel.h"

BruteForceSubdivisionMedialModel::BruteForceSubdivisionMedialModel() :
  SubdivisionMedialModel()
{
  xAtoms = NULL;
}

BruteForceSubdivisionMedialModel::~BruteForceSubdivisionMedialModel()
{
}

void
BruteForceSubdivisionMedialModel
::SetMesh(const MeshLevel &mesh, 
  const Vec &C, const Vec &u, const Vec &v,
  size_t nAtomSubs, size_t nCoeffSubs)
{
  // Call the parent method
  SubdivisionMedialModel::SetMesh(mesh, C, u, v, nAtomSubs, nCoeffSubs);

  // Initialize the gradient computer
  xLoopScheme.SetMesh(&mlAtom);

  // Initialize the derivative term arrays
  dt.resize(mlAtom.nVertices);
  dt_local.resize(mlAtom.nVertices);
}

void
BruteForceSubdivisionMedialModel
::ComputeAtoms(bool flagAllowErrors, const double *xHint)
{
  size_t i;

  // The first step is to compute the X and rho of the atoms based on the
  // coefficients. This step is performed in this class, not in the solver
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i]; a.X.fill(0.0); a.R = 0.0;

    // Compute the weighted sum of the coefficients
    ImmutableSparseMatrix<double>::RowIterator it = mlAtom.weights.Row(i);
    for( ; !it.IsAtEnd(); ++it)
      {
      size_t j = it.Column() << 2; double w = it.Value();
      a.X[0] += w * xCoefficients[j];
      a.X[1] += w * xCoefficients[j+1];
      a.X[2] += w * xCoefficients[j+2];
      a.R    += w * xCoefficients[j+3];
      }

    // Negative R should be disallowed
    if(!flagAllowErrors && a.R < 0.0)
      throw MedialModelException("Negative F in BruteForceModel");
    }

  // Get a pointer to the Loop scheme weights
  const LoopTangentScheme::WeightMatrix &W = xLoopScheme.GetWeightMatrix();

  // We must now compute the derivatives of X and R to make real atoms
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];

    // Compute the partial derivatives of F and X
    a.Xu.fill(0.0); a.Xv.fill(0.0); a.Ru = 0; a.Rv = 0;
    for(LoopTangentScheme::WeightMatrix::ConstRowIterator it = W.Row(i); 
      !it.IsAtEnd(); ++it)
      {
      MedialAtom &anbr = xAtoms[it.Column()];
      double wu = it.Value().w[0];
      double wv = it.Value().w[1];

      a.Xu += wu * anbr.X;
      a.Xv += wv * anbr.X;
      a.Ru += wu * anbr.R;
      a.Rv += wv * anbr.R;
      }

    // Compute the atom's differential geometry, normal and boundary
    a.G.SetOneJet(a.X.data_block(), a.Xu.data_block(), a.Xv.data_block());
    a.ComputeNormalVector();
    a.ComputeBoundaryAtomsUsingR(!mlAtom.IsVertexInternal(i));

    // For atoms that are on the edge, we correct Ru, so that |gradR|=1 and
    // Rv stays constant. This involves solving the system
    // g11 Ru Ru + 2 g12 Ru Rv + g22 Rv Rv = 1
    if(!mlAtom.IsVertexInternal(i))
      {
      double g12 = a.G.xCovariantTensor[0][1];
      double g22 = a.G.xCovariantTensor[1][1];
      double z = g22 - a.Rv * a.Rv;
      if(!flagAllowErrors && z < 0)
        throw MedialModelException("Excessive Rv in BruteForceModel");

      // Store the unfixed gradR magnitude in the atom
      a.xGradRMagSqrOrig = a.xGradRMagSqr;
      a.xNormalFactorOrig = a.xNormalFactor;

      // we store the square root of z for further derivative compn
      dt_local[i].sqrt_gz = sqrt(a.G.g * z);

      // Back up original Ru b/c we need it for derivative computation
      dt_local[i].Ru_orig = a.Ru;

      // here is the corrected Ru
      dt_local[i].Ru_fixed = a.Ru = (g12 * a.Rv - dt_local[i].sqrt_gz) / g22;
      
      // Compute atom using corrected Ru
      a.ComputeBoundaryAtomsUsingR(!mlAtom.IsVertexInternal(i));
      }

    if(!flagAllowErrors && !a.flagValid)
      {
      std::cerr << "Atom " << i << " is invalid in ComputeAtoms" << std::endl;
      std::cerr << "  Edge atom? " << a.flagCrest << std::endl;
      std::cerr << "  GradMagRSq " << a.xGradRMagSqr << std::endl;
      throw MedialModelException("Invalid Atom in BruteForceModel");
      }
    }
}

BruteForceSubdivisionMedialModel::Vec
BruteForceSubdivisionMedialModel::GetHintArray() const
{
  Vec xHint(1, 0.0);
  return xHint;
}

void 
BruteForceSubdivisionMedialModel
::SetVariationalBasis(const Mat &xBasis) 
{
  size_t i;

  // The non-varying terms are represented as a sparse array. Each entry in the
  // sparse array holds X and its partial derivatives, as well as F. 
  size_t nvar = xBasis.rows();

  // Allocate structure for holding data derived from the basis
  NonvaryingTermsMatrix::STLSourceType nvSource;

  // Get a pointer to the Loop scheme weights
  const LoopTangentScheme::WeightMatrix &W = xLoopScheme.GetWeightMatrix();

  // Iterate over all basis components
  for(size_t var = 0; var < nvar; var++)
    {
    // Get the current variation
    Vec xVariation = xBasis.get_row(var);

    // For each atom, we must decide if any of it's properties are affected by the
    // variation or not. The atom's properties depend on up to the second derivative
    // of X and R. The following array is used to label atoms as dependent or not

    // Also allocate a non-sparse array for computing the non-varying terms
    std::vector<NonvaryingAtomTerms> vTerms(mlAtom.nVertices);

    // Go through and mark the atoms whose X and R are affected
    size_t nc = 0;
    for(i = 0; i < mlAtom.nVertices; i++)
      {
      // Loop over 1-ring of neighbors
      for(ImmutableSparseMatrix<double>::RowIterator it = mlAtom.weights.Row(i);
        !it.IsAtEnd(); ++it)
        {
        // Get the index of the neighbor in the variation
        size_t j = it.Column() << 2; double w = it.Value();

        // Get the corresponding coefficients
        double x0 = xVariation[j];
        double x1 = xVariation[j+1];
        double x2 = xVariation[j+2];
        double r  = xVariation[j+3];

        // Check if there is any contribution from the variation to atom i
        if(x0 != 0.0 || x1 != 0.0 || x2 != 0.0 || r != 0.0)
          {
          // Set the atom's order to 0 (it's set to 3 by default)
          vTerms[i].order = 0;
          nc++;
          }

        vTerms[i].X[0] += w * x0;
        vTerms[i].X[1] += w * x1;
        vTerms[i].X[2] += w * x2;
        vTerms[i].R    += w * r;
        }
      }

    // Propagate the dependency flag to the 2-nd ring of neighbors. This
    // is because the neighbors are used in the computation of Xu and Xv
    for(int level = 0; level < 2; level++)
      {
      size_t nc = 0;
      std::vector<bool> vMark(mlAtom.nVertices, false);
      for(i = 0; i < mlAtom.nVertices; i++)
        {
        for(EdgeWalkAroundVertex it(&mlAtom, i); !it.IsAtEnd(); ++it)
          {
          if(vTerms[it.MovingVertexId()].order < 3)
            vMark[i] = true;
          }
        if(vMark[i]) nc++;
        }
      for(i = 0; i < mlAtom.nVertices; i++)
        if(vMark[i])
          vTerms[i].order = std::min(vTerms[i].order, level+1);
      }

    for(int ord = 0; ord < 3; ord++)
      {
      size_t nc = 0;
      for(i = 0; i < mlAtom.nVertices; i++) 
        if(vTerms[i].order == ord)
          nc++;
      }                                   

    // Next, precompute the first partial derivatives of X and R
    for(i = 0; i < mlAtom.nVertices; i++) 
  
      for(LoopTangentScheme::WeightMatrix::ConstRowIterator it = W.Row(i); 
        !it.IsAtEnd(); ++it)
        {
        double wu = it.Value().w[0];
        double wv = it.Value().w[1];
        SMLVec3d Xnbr = vTerms[it.Column()].X;
        double Rnbr = vTerms[it.Column()].R;

        vTerms[i].Xu += wu * Xnbr;
        vTerms[i].Xv += wv * Xnbr;
        vTerms[i].Ru += wu * Rnbr;
        vTerms[i].Rv += wv * Rnbr;
        }

    // Finally, place the marked nodes into a sparse structure for later use
    NonvaryingTermsMatrix::STLRowType nvRow;
    for(i = 0; i < mlAtom.nVertices; i++)
      if(vTerms[i].order < 3)
        nvRow.push_back(make_pair(i, vTerms[i]));
    nvSource.push_back(nvRow);
    }

  // The very last step is to initialize the sparse matrix
  this->xBasis.SetFromSTL(nvSource, mlAtom.nVertices);
}

void
BruteForceSubdivisionMedialModel
::BeginGradientComputation()
{
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    // Precompute the atom's terms
    xAtoms[i].ComputeCommonDerivativeTerms(dt[i]);
    
    // Now, our own thing
    if(!mlAtom.IsVertexInternal(i))
      {
      MedialAtom &a = xAtoms[i];
      double g12 = a.G.xCovariantTensor[0][1];
      double g22 = a.G.xCovariantTensor[1][1];
      double z = g22 - a.Rv * a.Rv;
      LocalDerivativeTerms &ld = dt_local[i];

      ld.w_Rv = (a.Rv*ld.sqrt_gz + g12*z)/(g22*z);
      ld.w_g12 = a.Rv / g22;
      ld.w_g22 = 
        ((g22 - 2 * a.Rv * a.Rv)*ld.sqrt_gz - 2*g12*a.Rv*z)/(2. * g22 * g22 * z);
      ld.w_g = - ld.sqrt_gz / (2 * a.G.g * g22);
      }
    }  
}

void
BruteForceSubdivisionMedialModel
::EndGradientComputation()
{
}

void
BruteForceSubdivisionMedialModel
::ComputeAtomVariationalDerivative(size_t iBasis, MedialAtom *dAtoms)
{
  // Iterator for selecting the atoms affected by the current variation
  NonvaryingTermsMatrix::RowIterator it;

  // Clear the derivative information in all atoms
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    dAtoms[i].SetAllDerivativeTermsToZero();
    dAtoms[i].order = 3;
    }

  // Iterate over the corresponding sparse matrix row to compute the 
  // relevant atoms
  size_t nc = 0;
  for(it = xBasis.Row(iBasis); !it.IsAtEnd(); ++it)
    {
    // Get the atom's index
    size_t j = it.Column();

    // Get the current atom and the derivative (which we are computing)
    MedialAtom &a = xAtoms[j];
    MedialAtom &da = dAtoms[j];

    // Label the atom as dependent
    da.order = it.Value().order;
    nc++;

    // Copy the constant terms. This is a little wasteful, but pretty much 
    // necessary in order to only hold one array of dAtoms in memory
    da.X   = it.Value().X;
    da.Xu  = it.Value().Xu;
    da.Xv  = it.Value().Xv;
    da.R   = it.Value().R;
    da.Ru  = it.Value().Ru;
    da.Rv  = it.Value().Rv;

    // Compute the metric tensor derivatives of the atom
    a.ComputeMetricTensorDerivatives(da);

    // For atoms that are on the edge, we correct Ru, so that |gradR|=1 and
    // Rv stays constant. This involves solving the system
    // g11 Ru Ru + 2 g12 Ru Rv + g22 Rv Rv = 1
    if(a.flagCrest)
      {
      LocalDerivativeTerms &ld = dt_local[j];

      // Use the original value of Ru for the first computation
      a.Ru = ld.Ru_orig;

      // Compute the derivatives of the boundary nodes using unfixed Ru
      a.ComputeBoundaryAtomDerivativesUsingR(da, dt[j]);

      // Record the values of |gradR| before applying the fix
      da.xGradRMagSqrOrig = da.xGradRMagSqr;
      da.xNormalFactorOrig = da.xNormalFactor;

      // Set Ru back to the fixed value
      a.Ru = ld.Ru_fixed;

      da.Ru = 
        ld.w_g   * da.G.g +
        ld.w_g12 * da.G.xCovariantTensor[0][1] +
        ld.w_g22 * da.G.xCovariantTensor[1][1] +
        ld.w_Rv  * da.Rv;

      // Recompute the derivatives of the boundary nodes
      a.ComputeBoundaryAtomDerivativesUsingR(da, dt[j]);
      }
    else
      {
      // Compute the derivatives of the boundary nodes
      a.ComputeBoundaryAtomDerivativesUsingR(da, dt[j]);
      }
    }

}

void
BruteForceSubdivisionMedialModel::
WriteToRegistry(Registry &R)
{
  SubdivisionMedialModel::WriteToRegistry(R);

  // Set the model subtype
  R["Grid.Model.SolverType"] << "BruteForce";
}


void
BruteForceSubdivisionMedialModel::
ReadFromRegistry(Registry &R)
{
}

