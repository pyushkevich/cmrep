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

  // Pass the mesh to the loop scheme
  xLoopScheme.SetMeshLevel(&mlAtom);

  // Compute the sparse matrix that gives u and v derivatives of 
  // functions/vectors over the mesh. These are going to be derivatives
  // with respect to the U and V arrays passed in, not to some arbitrary
  // local parameterization (which is referred to as p, q)
  SparseMat::STLSourceType srcWuv;
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    {
    // Generate the row
    SparseMat::STLRowType row;

    // Compute the weigths for this row. The weights arise in the following:
    //   dF_i/dp = Sum_{j \in N(i)} W^p_j F_j
    //   dF_i/du = Sum_{j \in N(i)} (dp/du W^p_j + dq/du W^p_j) F_j
    // So we are here computing (dp/du W^p_j + dq/du W^p_j) and storing them
    // as elements of the sparse matrix Wu. We also do the same for Wv
    
    // First, we are going to compute du/dp and dv/dp
    double dudp = xLoopScheme.GetPartialDerivative(0, i, uAtom.data_block());
    double dudq = xLoopScheme.GetPartialDerivative(1, i, uAtom.data_block());
    double dvdp = xLoopScheme.GetPartialDerivative(0, i, vAtom.data_block());
    double dvdq = xLoopScheme.GetPartialDerivative(1, i, vAtom.data_block());

    // Now compute the determinant of the Jacobian
    double detj = 1.0 / (dudp * dvdq - dudq * dvdp);

    // Now compute the inverse of the Jacobian
    double dpdu =   detj * dvdq;
    double dpdv = - detj * dudq;
    double dqdu = - detj * dvdp;
    double dqdv =   detj * dudp; 

    // Compute how much the vertex itself contributes to the weights
    row.push_back(
      make_pair(i, 
        make_pair(
          xLoopScheme.GetOwnWeight(0, i) * dpdu + 
          xLoopScheme.GetOwnWeight(1, i) * dqdu,
          xLoopScheme.GetOwnWeight(0, i) * dpdv +
          xLoopScheme.GetOwnWeight(1, i) * dqdv)));

    // Now, compute the contributing weights for each vertex
    for(EdgeWalkAroundVertex w(&mlAtom, i); !w.IsAtEnd(); ++w)
      {
      row.push_back(
        make_pair(w.MovingVertexId(), 
          make_pair(
            xLoopScheme.GetNeighborWeight(0, w) * dpdu + 
            xLoopScheme.GetNeighborWeight(1, w) * dqdu,
            xLoopScheme.GetNeighborWeight(0, w) * dpdv + 
            xLoopScheme.GetNeighborWeight(1, w) * dqdv)));
      }

    // Store the rows
    srcWuv.push_back(row);
    }

  // Set the derivative matrices
  Wuv.SetFromSTL(srcWuv, mlAtom.nVertices);
}

void
BruteForceSubdivisionMedialModel
::ComputeAtoms(const double *xHint)
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
      a.X += w * xCoefficients.extract(3, j);
      a.R += w * xCoefficients[j+3];
      }

    // Set F = R^2
    a.F = a.R * a.R;
    }

  // We must now compute the derivatives of X and R to make real atoms
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];

    // Compute the partial derivatives of F and X
    a.Xu.fill(0.0); a.Xv.fill(0.0); a.Fu = 0; a.Fv = 0;
    for(SparseMat::RowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
      {
      MedialAtom &anbr = xAtoms[it.Column()];
      double wu = it.Value().first;
      double wv = it.Value().second;

      a.Xu += wu * anbr.X;
      a.Xv += wv * anbr.X;
      a.Fu += wu * anbr.F;
      a.Fv += wv * anbr.F;
      }
    }

  // On the third pass, we compute the second order partial derivatives
  for(i = 0; i < mlAtom.nVertices; i++)
    {
    // Set up i-th atom
    MedialAtom &a = xAtoms[i];

    // Compute the partial derivatives of F and X
    a.Xuu.fill(0.0); a.Xuv.fill(0.0); a.Xvv.fill(0.0);
    a.Fuu = 0.0; a.Fuv = 0.0; a.Fvv = 0.0;
    SMLVec3d Xvu; Xvu.fill(0.0);
    for(SparseMat::RowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
      {
      MedialAtom &anbr = xAtoms[it.Column()];
      double wu = it.Value().first;
      double wv = it.Value().second;

      a.Xuu += wu * anbr.Xu;
      a.Xuv += 0.5 * (wu * anbr.Xv + wv * anbr.Xu);
      a.Xvv += wv * anbr.Xv;
      a.Fuu += wu * anbr.Fu;
      a.Fuv += 0.5 * (wu * anbr.Fv + wv * anbr.Fu);
      a.Fvv += wv * anbr.Fv;
      }

    // Compute the atom's differential geometry, normal and boundary
    a.ComputeDifferentialGeometry();
    a.ComputeNormalVector();
    a.ComputeBoundaryAtoms(!mlAtom.IsVertexInternal(i));

    // Turned off since we are not using this anywhere
    // a.ComputeBoundaryCurvature();
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
        size_t j = it.Column() << 2; 

        // Get the neighbor's X and R combined, scale by the neighbor's weight
        SMLVec4d XR = it.Value() * xVariation.extract(4,j);

        // Check if there is any contribution from the variation to atom i
        if(XR.squared_magnitude() > 0.0)
          {
          // Set the atom's order to 0 (it's set to 3 by default)
          vTerms[i].order = 0;
          nc++;
          }

        vTerms[i].X += XR.extract(3);
        vTerms[i].R += XR[3];
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
      // cout << "Order " << ord << ": " << nc << " atoms. " << endl;
      }

    // Next, precompute the first partial derivatives of X and R
    for(i = 0; i < mlAtom.nVertices; i++) 
      for(SparseMat::RowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
        {
        double wu = it.Value().first;
        double wv = it.Value().second;
        SMLVec3d Xnbr = vTerms[it.Column()].X;

        vTerms[i].Xu += wu * Xnbr;
        vTerms[i].Xv += wv * Xnbr;
        }

    // Next, precompute the second order partial derivatives of X and R
    for(i = 0; i < mlAtom.nVertices; i++) 
      for(SparseMat::RowIterator it = Wuv.Row(i); !it.IsAtEnd(); ++it)
        {
        double wu = it.Value().first;
        double wv = it.Value().second;
        SMLVec3d XUnbr = vTerms[it.Column()].Xu;
        SMLVec3d XVnbr = vTerms[it.Column()].Xv;

        vTerms[i].Xuu += wu * XUnbr;
        vTerms[i].Xuv += 0.5 * (wu * XVnbr + wv * XUnbr);
        vTerms[i].Xvv += wv * XVnbr;
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
  // Precompute common terms for atom derivatives
  dt = new MedialAtom::DerivativeTerms[mlAtom.nVertices];
  for(size_t i = 0; i < mlAtom.nVertices; i++)
    xAtoms[i].ComputeCommonDerivativeTerms(dt[i]);
}

void
BruteForceSubdivisionMedialModel
::EndGradientComputation()
{
  // Precompute common terms for atom derivatives
  delete dt;
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

    // Label the atom as dependent
    dAtoms[j].order = it.Value().order;
    nc++;

    // Copy the constant terms. This is a little wasteful, but pretty much 
    // necessary in order to only hold one array of dAtoms in memory
    dAtoms[j].R   = it.Value().R;
    dAtoms[j].X   = it.Value().X;
    dAtoms[j].Xu  = it.Value().Xu;
    dAtoms[j].Xv  = it.Value().Xv;
    dAtoms[j].Xuu = it.Value().Xuu;
    dAtoms[j].Xuv = it.Value().Xuv;
    dAtoms[j].Xvv = it.Value().Xvv;

    // Compute the derivative of F = R^2 at the atom
    dAtoms[j].F = 2.0 * dAtoms[j].R * xAtoms[j].R;
    }
  // cout << "Dep Atoms : " << nc << " of " << mlAtom.nVertices << endl;

  // Second loop to compute partials of F
  for(it = xBasis.Row(iBasis); !it.IsAtEnd(); ++it)
    {
    // Get the atom's index
    size_t j = it.Column();

    // Get the current atom and the derivative (which we are computing)
    MedialAtom &a = xAtoms[j];
    MedialAtom &da = dAtoms[j];

    // Compute the partial derivatives of Fu
    da.Fu = da.Fv = 0.0;
    for(SparseMat::RowIterator wit = Wuv.Row(j); !wit.IsAtEnd(); ++wit)
      {
      double Fnbr = dAtoms[wit.Column()].F;
      double wu = wit.Value().first;
      double wv = wit.Value().second;

      da.Fu += wu * Fnbr;
      da.Fv += wv * Fnbr;
      }

    // Compute the metric tensor derivatives of the atom
    a.ComputeMetricTensorDerivatives(da);
    a.ComputeChristoffelDerivatives(da);

    // Compute the derivatives of the boundary nodes
    a.ComputeBoundaryAtomDerivatives(da, dt[j]);
    }

  // Third loop to compute second order partials of F
  for(it = xBasis.Row(iBasis); !it.IsAtEnd(); ++it)
    {
    // Get the atom's index
    size_t j = it.Column();

    // Get the current atom and the derivative (which we are computing)
    MedialAtom &a = xAtoms[j];
    MedialAtom &da = dAtoms[j];

    // Compute the derivative of mean and gauss curvatures
    da.Fuu = da.Fuv = da.Fvv = 0.0;
    for(SparseMat::RowIterator wit = Wuv.Row(j); !wit.IsAtEnd(); ++wit)
      {
      MedialAtom &danbr = dAtoms[wit.Column()];
      double wu = wit.Value().first;
      double wv = wit.Value().second;

      da.Fuu += wu * danbr.Fu;
      da.Fvv += wv * danbr.Fv;
      da.Fuv += 0.5 * (wu * danbr.Fv + wv * danbr.Fu);
      }

    // Compute things in the atom that depend on second derivatives
    // a.ComputeBoundaryCurvatureDerivative(da);
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

