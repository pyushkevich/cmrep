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
  saEncroach.Reset();
  saFalseEncroach.Reset();

  // Do a quadratic search for intersections
  for(MedialBoundaryPointIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {
    size_t ib = it.GetIndex();
    BoundaryAtom bat = GetBoundaryPoint(it, S->xAtoms);
    double xb = bat.X[0], yb = bat.X[1], zb = bat.X[2];

    // Set the scale for the encroachment calculation
    double scale = 10.0;

    // Set the cutoff for d^2/R^2, at which the penalty is so small that
    // it just does not matter (i.e., 10^-10). It looks like taking 
    // 4.5 / scale gives us the right cutoff
    double cutoff = (1 + 4.5 / scale) * (1 + 4.5 / scale);    

    // Loop over all medial atoms (this is where the kd-tree could help)
    for(size_t ia = 0; ia < S->nAtoms; ia++)
      {
      MedialAtom &A = S->xAtoms[ia];
      double dx = xb - A.X[0], dy = yb - A.X[1], dz = zb - A.X[2];
      double d2 = dx * dx + dy * dy + dz * dz;
      if(d2 < cutoff * A.F && ia != it.GetAtomIndex())
        {
        // The distance is a vector between -1 and 1, with 1 being total encroachment
        // and -1 being a whole 'r' distance away
        double d = sqrt(d2);
        double u = 1.0 - d / A.R;

        // Use the ERF function to compute the penalty. The slope of the ERF function
        // should be controlled by a slope variable
        double t = scale * (u - 0.2);
        double z = 0.5 * erf(t) + 0.5;

        // Insert into the sparse array structure
        if(flagGradient)
          {
          // Compute the derivative of ERF wrt u
          double dzdu = 0.5 * scale * vnl_math::two_over_sqrtpi / exp(t * t);

          // Create the structure for faster gradient computation          
          Intersection p;
          p.a = (- dzdu / (d * A.R)) * Vec(dx,dy,dz);
          p.b = dzdu * d / A.F;
          xMuteMat[ib].push_back(make_pair(ia, p));
          }

        // Associate the penalty with the term
        saPenalty.Update(z);
        if(u >= 0.0)
          saEncroach.Update(u);
        else
          saFalseEncroach.Update(u);
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
  sout << "    number of encoachments   : " << saEncroach.GetCount() << endl; 
  sout << "    mean encroachment        : " << saEncroach.GetMean() << endl; 
  sout << "    max encroachment         : " << saEncroach.GetMax() << endl; 
  sout << "    min encroachment         : " << saEncroach.GetMin() << endl; 
  sout << "    number of false encr.    : " << saFalseEncroach.GetCount() << endl; 
  sout << "    mean false encr.         : " << saFalseEncroach.GetMean() << endl; 
  sout << "    max false encr.          : " << saFalseEncroach.GetMax() << endl; 
  sout << "    min false encr.          : " << saFalseEncroach.GetMin() << endl; 
  sout << "    total penalty            : " << saPenalty.GetSum() << endl;
}

template class ImmutableSparseArray<DiffeomorphicEnergyTerm::Intersection>;
