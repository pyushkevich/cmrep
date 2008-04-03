#include "PrincipalComponents.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

bool ReadMatrixFile(vnl_matrix<double> &mat, const char *file)
{
  FILE *fin = fopen(file,"rt");
  if(fin == NULL) return false;

  size_t r = 0, c = 0;
  fscanf(fin, "# MATRIX %d %d\n", &r, &c);
  if(r == 0 || c == 0)
    { fclose(fin); return false; }

  mat.set_size(r, c);
  for(size_t i = 0; i < r; i++)
    for(size_t j = 0; j < c; j++)
      fscanf(fin, "%lg", &mat[i][j]);

  fclose(fin);
  return true;
}

PrincipalComponents::PrincipalComponents(const Mat &A)
{
  size_t i, j;
  m = A.columns();
  n = A.rows();

  // Compute the mean
  mu.set_size(m); mu.fill(0.0);
  for(i = 0; i < n; i++) for(j = 0; j < m; j++)
    mu[j] += A(i,j);
  mu /= n;

  // Compute the zero mean matrix
  Mat Z(n, m);
  for(i = 0; i < n; i++) for(j = 0; j < m; j++)
    Z(i,j) = A(i,j) - mu[j];

  // Split, depending on if m > n or not
  if(m > n)
    {
    // Compute the covariance matrix Z Zt
    Mat Zt = Z.transpose();
    Mat K = ( Z * Zt ) / n;

    // Get the eigensystem of K
    vnl_symmetric_eigensystem<double> eig(K);

    // Take the eigenvalues
    lambda.set_size(n); Q.set_size(n, m); 
    lambda.fill(0.0); Q.fill(0.0);
    
    for(i = 0; i < n; i++) 
      {
      lambda[i] = eig.get_eigenvalue(n - 1 - i);
      if(lambda[i] > 0.0)
        Q.set_row(i, (Zt * eig.get_eigenvector(n - 1 - i)).normalize());
      else break;
      }
    }
  else
    {
    // Compute the real covariance matrix At Z
    Mat K = (Z.transpose() * Z) / n;
    vnl_symmetric_eigensystem<double> eig(K);

    // Copy the eigenvalues and vectors
    lambda.set_size(m); Q.set_size(m,m); 
    lambda.fill(0.0); Q.fill(0.0);
    
    for(i = 0; i < m; i++) 
      {
      lambda[i] = eig.get_eigenvalue(m - 1 - i);
      if(lambda[i] > 0.0)
        Q.set_row(i, eig.get_eigenvector(n - 1 - i));
      else break;
      }
    }

  // Finish up
  nModes = i;
  sdev = lambda.apply(sqrt);
}

/*
PrincipalComponents::LeaveOneOutAnalysis(const Mat &A)
{
  // Repeat for each element of A
  unsigned int i, n = A.rows();
  for(i = 0; i < n; i++)
    {
    // Create a matrix with row i taken out
    Mat B(A.rows()-1, A.columns());
    if(i > 0) B.update(0,0,A.get_n_rows(0, i));
    if(i < n-1) B.update(i,0,A.get_n_rows(i+1,n-i-1));

    // Compute the PCA of B
    PrincipalComponents p(B);

    // Project the left out guy on the principal components
    p.

    }
}
*/

