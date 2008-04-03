#ifndef __Procrustes_h_
#define __Procrustes_h_

#include <vnl/vnl_matrix.h>

/**
 * Compute procrustes match from A to B 
 */
void ExtendedProcrustesAnalysis(
  const vnl_matrix<double> &A, const vnl_matrix<double> &B, 
  vnl_matrix<double> &Ahat, vnl_matrix<double> &R,
  vnl_vector<double> &t, double &s);

void GeneralizedProcrustesAnalysis(size_t m, 
  vnl_matrix<double> *A, vnl_matrix<double> *R,
  vnl_vector<double> *t, double *s);

void TestProcrustes(const vnl_matrix<double> &A);

#endif
