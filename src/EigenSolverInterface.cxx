#include "EigenSolverInterface.h"
#include <iostream>
#include <cstdlib>
#include <cstring>

#include <Eigen/Sparse>

using namespace std;

void 
EigenSolverInterface
::SymbolicFactorization(size_t n, int *idxRows, int *idxCols, double *xMatrix)
{
  // Create a sparse matrix map
  int nnz = idxRows[n];
  m_SparseMatrix = SparseMap(n, n, nnz, idxRows, idxCols, xMatrix);
}

void 
EigenSolverInterface
::SymbolicFactorization(const ImmutableSparseMatrix<double> &mat)
{

}

void 
EigenSolverInterface
::NumericFactorization(const double *xMatrix)
{

}

void 
EigenSolverInterface
::Solve(double *xRhs, double *xSoln)
{
}

void                              
EigenSolverInterface
::Solve(size_t nRHS, double *xRhs, double *xSoln)
{
}
