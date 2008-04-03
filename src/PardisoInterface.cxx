#include "PardisoInterface.h"
#include <iostream>

using namespace std;

// BLAS/PARDISO references
#ifdef PARDISO_FOUND

extern "C" {
  void pardisoinit_(size_t *, int *, int *);
  void pardiso_(size_t *, int *, int *, int *, int *, int *, double *, int *, int *, 
    int *, int *, int *, int *, double *, double *, int*);
}

#else

void pardisoinit_(size_t *, int *, int *)
  { cerr << "Pardiso unavailable; exiting." << endl; exit(-1); }

void pardiso_(size_t *, int *, int *, int *, int *, int *, double *, int *, int *, 
    int *, int *, int *, int *, double *, double *, int*)
  { cerr << "Pardiso unavailable; exiting." << endl; exit(-1); }

#endif


GenericRealPARDISO::GenericRealPARDISO(int type)
{
  // Set the type of matrix to unsymmetric real
  MTYPE = type; 

  // Clear the parameter array
  memset(IPARM, 0, sizeof(int) * 64);

  // Initialize PARDISO to default values
  pardisoinit_(PT,&MTYPE,IPARM);

  // Specify the number of processors on the system (1)
  IPARM[2] = 1;

  flagPardisoCalled = false;

  // Whether the index arrays are owned
  flagOwnIndexArrays = false;
}

void 
GenericRealPARDISO
::ResetIndices()
{
  if(flagOwnIndexArrays)
    { delete idxRows; delete idxCols; }
  flagOwnIndexArrays = false;
}

void 
GenericRealPARDISO
::SymbolicFactorization(size_t n, int *idxRows, int *idxCols, double *xMatrix)
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = 11, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    xMatrix, idxRows, idxCols,
    NULL, &NRHS, IPARM, &MSGLVL, NULL, NULL, &ERROR);

  // Record the parameter for next phase
  ResetIndices();
  this->idxCols = idxCols;
  this->idxRows = idxRows;
  this->n = n;

  // Set the flag so we know that pardiso was launched before
  flagPardisoCalled = true;
}

void 
GenericRealPARDISO
::SymbolicFactorization(const ImmutableSparseMatrix<double> &mat)
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = 11, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 

  // Init the index arrays
  ResetIndices();

  // We are going to own the indices
  flagOwnIndexArrays = true;

  // The arrays have to be incremented by one before calling PARDISO
  idxCols = new int[mat.GetNumberOfRows() + 1];
  for(size_t i = 0; i <= mat.GetNumberOfRows(); i++)
    idxRows[i] = (int)(1 + mat.GetRowIndex()[i]);

  idxCols = new int[mat.GetNumberOfColumns()];
  for(size_t j = 0; j <= mat.GetNumberOfColumns(); j++)
    idxCols[j] = (int)(1 + mat.GetColIndex()[j]);
  
  // Perform the symbolic factorization phase
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    const_cast<double *>(mat.GetSparseData()), idxRows, idxCols,
    NULL, &NRHS, IPARM, &MSGLVL, NULL, NULL, &ERROR);

  // Record the parameter for next phase
  this->n = mat.GetNumberOfRows();

  // Set the flag so we know that pardiso was launched before
  flagPardisoCalled = true;
}

void 
GenericRealPARDISO
::NumericFactorization(const double *xMatrix)
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = 22, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    const_cast<double *>(xMatrix), idxRows, idxCols,
    NULL, &NRHS, IPARM, &MSGLVL, NULL, NULL, &ERROR);

  // Record the parameter for next phase
  this->xMatrix = xMatrix;
}

void 
GenericRealPARDISO
::Solve(double *xRhs, double *xSoln)
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = 33, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    const_cast<double *>(xMatrix), idxRows, idxCols,
    NULL, &NRHS, IPARM, &MSGLVL, xRhs, xSoln, &ERROR);
}

void 
GenericRealPARDISO
::Solve(size_t nRHS, double *xRhs, double *xSoln)
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = 33, N = n, NRHS = nRHS, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
    const_cast<double *>(xMatrix), idxRows, idxCols,
    NULL, &NRHS, IPARM, &MSGLVL, xRhs, xSoln, &ERROR);
}

GenericRealPARDISO::
~GenericRealPARDISO()
{
  // Set the various parameters
  int MAXFCT = 1, MNUM = 1, PHASE = -1, N = n, NRHS = 1, MSGLVL = 0, ERROR = 0; 
  
  // Perform the symbolic factorization phase
  if(flagPardisoCalled)
    pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &PHASE, &N, 
      const_cast<double *>(xMatrix), idxRows, idxCols,
      NULL, &NRHS, IPARM, &MSGLVL, NULL, NULL, &ERROR);

  // Reset the arrays
  ResetIndices();
}

