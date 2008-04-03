#ifndef __PrinicipalComponents_h_
#define __PrinicipalComponents_h_

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

bool ReadMatrixFile(vnl_matrix<double> &mat, const char *file);

class PrincipalComponents
{
public:
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;

  PrincipalComponents(const Mat &A);

  /** Get the number of principal modes */
  size_t GetNumberOfModes() 
    { return nModes; }

  /** Get the eigenvalues in decreasing order */
  double GetEigenvalue(size_t i) const
    { return lambda[i]; }

  /** Get the eigenvectors in decreasing order of eigenvalue */
  Vec GetEigenvector(size_t i) const
    { return Q.get_row(i); }

  /** Get the shape at a certain position in PCA space. The coordinates in this
   * space are the eigenvectors of PCA and the units are standard deviations
   * (sqrt of the eigenvalues). The input vector can be shorter than the number
   * of modes; the remaining modes will be treated as zero */
  Vec MapToFeatureSpace(const Vec &xComponent) const
    {
    Vec x = mu;
    for(size_t i = 0; i < xComponent.size(); i++)
      x += xComponent[i] * sdev[i] * GetEigenvector(i);
    return x;
    }

  /**
   * Map a set of PCA coefficients to a feature vector with zero mean. This is
   * in derivative/gradient computations
   */
  Vec MapToFeatureSpaceZeroMean(const Vec &xComponent) const
    {
    Vec x(m, 0.0);
    for(size_t i = 0; i < xComponent.size(); i++)
      x += xComponent[i] * sdev[i] * GetEigenvector(i);
    return x;
    }

  const Vec &GetMean() const
    { return mu; }

  /** Project a point in feature (shape) space into the coefficient space */
  Vec MapToCoefficientSpace(const Vec &xFeature, size_t nCoeff) const
    {
    Vec z(m, 0.0);
    if(nCoeff > nModes) nCoeff = nModes;
    for(size_t i = 0; i < nCoeff; i++)
      z[i] = dot_product(xFeature - mu, GetEigenvector(i)) / sdev[i];
    return z;
    }

  /** Get a feature space vector along one of the principal modes. Useful for 
   * animating the PCA modes */
  Vec MapToFeatureSpace(size_t iComponent, double xComponentValue)
    {
    return mu + xComponentValue * sdev[iComponent] * GetEigenvector(iComponent);
    }

  /**
   * Map a set of PCA coefficients to a feature vector with zero mean. This is
   * in derivative/gradient computations
   */
  Vec MapToFeatureSpaceZeroMean(size_t iComponent, double xComponentValue)
    {
    return xComponentValue * sdev[iComponent] * GetEigenvector(iComponent);
    }

  /**
   * Perform leave-one-out analysis, i.e., how well the PCA can predict the 
   * individual members.
   */
  // static void LeaveOneOutAnalysis(const Mat &A);

private:
  // The mean vector
  Vec mu;

  // The principal component vectors (rows)
  Mat Q;

  // The eigenvalues
  Vec lambda, sdev;

  // The dimensions of the data
  size_t m, n, nModes;
};

#endif
