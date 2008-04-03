#include "Procrustes.h"
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_trace.h>
#include <vnl/vnl_rotation_matrix.h>
#include <vnl/vnl_math.h>
#include <iostream>

using namespace std;

typedef vnl_matrix<double> Mat;
typedef vnl_vector<double> Vec;

/**
 * Subtract the mean of the rows from each row in the matrix
 */
void SubtractMean(Mat &A, Vec &meanA)
{
  // Get the data dimensions
  size_t k = A.cols(), p = A.rows(), i, j;

  // Initialize the mean
  meanA.set_size(k);
  meanA.fill(0.0);

  // Compute the mean
  for(i = 0; i < p; i++) for(j = 0; j < k; j++)
    meanA[j] += A[i][j];
  meanA /= p;

  // Subtract the mean
  for(i = 0; i < p; i++) for(j = 0; j < k; j++)
    A[i][j] = A[i][j] - meanA[j];
}

/**
 * Compute procrustes match from A to B 
 */
void ExtendedProcrustesAnalysis(const Mat &A, const Mat &B, 
  Mat &Ahat, Mat &R, Vec &t, double &s)
{
  // Get the data dimensions
  size_t k = A.cols(), p = A.rows();

  // Subtract out the mean from matrix A
  Mat A0 = A;
  Vec meanA;
  SubtractMean(A0, meanA);

  // Compute the matrix that should be SVD'd
  Mat S = A0.transpose() * B;
  
  // Compute the SVD of S
  vnl_svd<double> svd(S);

  // Compute the rotation matrix
  R = svd.U() * svd.V().transpose();

  // Compute the scale factor
  s = vnl_trace(R.transpose() * S) / vnl_trace(A0.transpose() * A);

  // Rotate A0 into B0
  Ahat = s * A * R;

  // Compute the translation factor
  Mat D = B - Ahat;
  SubtractMean(D, t);

  // Translate into B
  for(size_t i = 0; i < p; i++) for(size_t j = 0; j < k; j++)
    Ahat[i][j] += t[j];

  // That's it!
  cout << "." << flush;
}

/** Compute the Procrustes error of a set of shapes */
double ComputeProcrustesError(size_t m, Mat *A)
{
  size_t k = A[0].cols(), p = A[0].rows(), i, j;

  // Compute the overall error
  Mat E(k, k, 0.0);
  for(i = 0; i < m; i++) for(j = i+1; j < m; j++)
    {
    Mat D = A[i] - A[j];
    E += D.transpose() * D;
    }

  return sqrt(vnl_trace(E) / (m * m * p));
}

/**
 * Compute the generalized procrustes alignment of a set of shapes
 */
void GeneralizedProcrustesAnalysis(size_t m, Mat *A, Mat *R, Vec *t, double *s)
{
  // Get the data dimensions
  size_t k = A[0].cols(), p = A[0].rows(), i;

  cerr << "Procrustes matrix: " << endl << A[2] << endl;

  // Make a copy of all the original objects
  Mat *B = new Mat[m];
  for(i = 0; i < m; i++)
    B[i] = A[i];

  // Report the initial procrustes error
  cout << "Initial Procrustes error: " << ComputeProcrustesError(m,B) << endl;
  
  // Start with the initial Procrustes mean
  Mat C(A[0]);

  // Continue until the Procrustes mean converges
  bool done = false;
  int niter = 0;
  while(!done && niter++ < 40)
    {
    // Map every matrix onto the centroid
    for(i = 0; i < m; i++)
      {
      ExtendedProcrustesAnalysis(B[i], C, A[i], R[i], t[i], s[i]);
      }

    // Compute the new Procrustes mean
    Mat C1(p, k, 0.0);
    for(i = 0; i < m; i++) C1 += A[i];
    C1 /= m;

    // Align the mean to shape A0
    Mat R0(k, k, 0.0); Vec t0(k, 0.0); double s0 = 0;
    ExtendedProcrustesAnalysis(C1, B[0], C1, R0, t0, s0);

    // Report
    cout << "*" << endl;
    cout << "Current Procrustes Error: " << ComputeProcrustesError(m, A) << endl;
    cout << "Convergence Factor: " << (C1 - C).array_inf_norm() << endl;

    // Compare the procrustes means
    if((C1 - C).array_inf_norm() < 1.0e-8)
      done = true;
    C = C1;
    }

  // Clean up
  delete[] B;
}

/** Test the procrustes method */
void TestProcrustes(const Mat &A)
{
  // Get the data dimensions
  const size_t m = 100;
  size_t k = A.cols(), p = A.rows(), i, j;

  // Generate a list of rotations, translations and scales
  Mat R[m], R0[m];
  Vec t[m], t0[m];
  double s[m], s0[m];
  Mat B[m];

  // Compute the random rotations translations and scales
  for(i = 0; i < m; i++)
    {
    // Initialize the rotation matrix and an identity
    Mat I(k, k, 0.0); I.fill_diagonal(0.0);
    
    R[i].set_size(k, k); 
    R[i].fill(0.0);
    R[i].fill_diagonal(0.0);

    t[i].set_size(k);

    // Compute the rotation matrix and the translation for each coordinate
    for(j = 0; j < k; j++)
      {
      double aRot = rand() * vnl_math::pi * 0.5 / RAND_MAX;
      Vec xRot = I.get_row(j) * aRot;
      Mat Rj(k, k);
      vnl_rotation_matrix(xRot, Rj);
      R[i] *= Rj;
      t[i][j] = rand() * 4.0 / RAND_MAX - 2.0;      
      }
    s[i] = rand() * 3.0 / RAND_MAX;

    // Create the matrix
    B[i] = s[i] * A * R[i];
    for(j = 0; j < p; j++) for(size_t q = 0; q < k; q++)
      B[i][j][q] += t[i][q];

    }

  // Compute the Procrustes alignment
  GeneralizedProcrustesAnalysis(m, B, R0, t0, s0);
}
