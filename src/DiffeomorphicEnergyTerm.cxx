#include "DiffeomorphicEnergyTerm.h"
#include "SparseMatrix.txx"

DiffeomorphicEnergyTerm
::DiffeomorphicEnergyTerm(GenericMedialModel *model)
{
  xMuteMat.resize(model->GetNumberOfBoundaryPoints());
}

DiffeomorphicEnergyTerm
::~DiffeomorphicEnergyTerm()
{
}

double 
DiffeomorphicEnergyTerm
::UnifiedComputeEnergy(SolutionData *S, bool flagGradient)
{
  // Reinitialize the sparse matrix if in gradient mode
  if(flagGradient)
    {
    for(MutableSparseMatrix::iterator it = xMuteMat.begin(); 
      it != xMuteMat.end(); ++it)
      {
      it->clear();
      }
    }

  // Reset the penalty
  saPenalty.Reset();

  // Do a quadratic search for intersections
  for(MedialBoundaryPointIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    size_t ib = it.GetIndex();
    BoundaryAtom bat = GetBoundaryPoint(it, S->xAtoms);
    double xb = bat.X[0], yb = bat.X[1], zb = bat.X[2];

    // Loop over all medial atoms (this is where the kd-tree could help)
    for(size_t ia = 0; ia < S->nAtoms; ia++)
      {
      MedialAtom &A = S->xAtoms[ia];
      double dx = xb - A.X[0], dy = yb - A.X[1], dz = zb - A.X[2];
      double d2 = dx * dx + dy * dy + dz * dz;
      if(d2 < A.F && ia != it.GetAtomIndex())
        {
        // Get the actual distance
        double d = sqrt(d2);
        double sqpen = d / A.R - 1.0;

        // Insert into the sparse array structure
        if(flagGradient)
          {
          // Create the structure for faster gradient computation
          Intersection p;
          p.a = (2.0 * sqpen / (d * A.R)) * Vec(dx, dy, dz);
          p.b = 2.0 * sqpen * (- d / A.F); 
          xMuteMat[ib].push_back(make_pair(ia, p));
          }

        // Associate the penalty with the term
        saPenalty.Update(sqpen * sqpen);
        }
      }
    }

  // Generate the sparse matrix
  if(flagGradient)
    {
    xDistMat.SetFromSTL(xMuteMat, S->nAtoms);
    }

  // Return the total penalty (shouldn't there be a large multiplier?)
  return saPenalty.GetSum();
}

double
DiffeomorphicEnergyTerm
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // Get the column index and row index
  size_t *rows = xDistMat.GetRowIndex();
  size_t *cols = xDistMat.GetColIndex();
  Intersection *ints = xDistMat.GetSparseData();
  size_t n = xDistMat.GetNumberOfRows();

  double xTotalDeriv = 0.0;

  // Iterate over all elements of the intersection
  for(MedialBoundaryPointIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    size_t ib = it.GetIndex();
    for(size_t j = rows[ib]; j < rows[ib+1]; j++)
      {
      // Now we have an atom-boundary intersection
      size_t ia = cols[j];

      // Compute the derivative
      Intersection &p = ints[j];
      double deriv = 
        dot_product(p.a, GetBoundaryPoint(it, dS->xAtoms).X - dS->xAtoms[ia].X)
        + p.b * dS->xAtoms[ia].R;

      // Accumulate
      xTotalDeriv += deriv;
      }
    }

  return xTotalDeriv;
}

void
DiffeomorphicEnergyTerm
::PrintReport(ostream &sout)
{
  sout << "  Diffeomorphic Penalty Term : " << endl;
  sout << "    number of encoachments   : " << saPenalty.GetCount() << endl; 
  sout << "    RMS encroachment         : " << sqrt(saPenalty.GetMean()) << endl; 
  sout << "    max encroachment         : " << sqrt(saPenalty.GetMax()) << endl; 
  sout << "    min encroachment         : " << sqrt(saPenalty.GetMin()) << endl; 
  sout << "    total penalty            : " << saPenalty.GetSum() << endl;
}

template class ImmutableSparseArray<DiffeomorphicEnergyTerm::Intersection>;
