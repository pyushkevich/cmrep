#ifndef __PardisoInterface_h_
#define __PardisoInterface_h_

#include <iostream>
#include "SparseMatrix.h"

class GenericRealPARDISO
{
public:
  // Factor the system for arbitrary right hand sides and matrices of the same
  // non-zer element structure
  void SymbolicFactorization(size_t n, int *idxRows, int *idxCols, double *xMatrix);

  // Perform symbolic factorization given a matrix
  void SymbolicFactorization(const ImmutableSparseMatrix<double> &mat);

  // Factor the system for a specific matrix, but arbitrary right hand side
  void NumericFactorization(const double *xMatrix);

  // Numeric factorization using sparse matrix datatype
  void NumericFactorization(const ImmutableSparseMatrix<double> &mat)
    { NumericFactorization(mat.GetSparseData()); }

  // Solve the system for the given right hand side, solution in xSoln
  void Solve(double *xRhs, double *xSoln);

  // Solve the system for a number of right hand sides, if the second vector
  // is NULL, will solve in-place
  void Solve(size_t nRHS, double *xRhs, double *xSoln);

protected:
  
  // Constructor, takes the problem type
  GenericRealPARDISO(int mtype);

  // Destructor
  virtual ~GenericRealPARDISO();

  // Reset the index arrays()
  void ResetIndices();

private:
  /** Internal data for PARDISO */
  size_t PT[64];
  int MTYPE;
  int IPARM[64];

  // Storage for data in intermediate steps
  int n, *idxRows, *idxCols;
  const double *xMatrix;
  bool flagPardisoCalled, flagOwnIndexArrays;
};

class UnsymmetricRealPARDISO : public GenericRealPARDISO
{
public:
  // Initialize the solver 
  UnsymmetricRealPARDISO() : GenericRealPARDISO(11) {};
};

class SymmetricPositiveDefiniteRealPARDISO : public GenericRealPARDISO
{
public:
  // Initialize the solver 
  SymmetricPositiveDefiniteRealPARDISO() : GenericRealPARDISO(2) {};
};

#endif //__PardisoInterface_h_
