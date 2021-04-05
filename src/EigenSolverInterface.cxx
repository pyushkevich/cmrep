#include "EigenSolverInterface.h"
#include <iostream>
#include <cstdlib>
#include <cstring>

#include <Eigen/Sparse>

using namespace std;

template <class TIndex>
class EigenSolverInterfaceInternal
{
public:

  void SetMatrix(size_t n, TIndex *idxRows, TIndex *idxCols, double *xMatrix)
  {
    int nnz = idxRows[n];
    m_SparseMatrix = new SparseMap(n, n, nnz, idxRows, idxCols, xMatrix);
  }

  void UpdateMatrix(double *xMatrix)
  {
    assert(m_SparseMatrix);
    if(m_SparseMatrix->valuePtr() != xMatrix)
      {
      SparseMap *temp = new SparseMap(m_SparseMatrix->rows(), m_SparseMatrix->cols(), m_SparseMatrix->nonZeros(),
                                      m_SparseMatrix->outerIndexPtr(), m_SparseMatrix->innerIndexPtr(),
                                      const_cast<double *>(xMatrix));
      delete m_SparseMatrix;
      m_SparseMatrix = temp;
      }
  }

  void Compute()
  {
    m_Solver.compute(*m_SparseMatrix);
  }

  void Solve(unsigned int nRHS, double *xRHS, double *xSoln)
  {
    unsigned int n = m_SparseMatrix->cols();
    for(unsigned int i = 0; i < nRHS; i++)
      {
      Eigen::Map<Eigen::VectorXd> mapRHS(xRHS + i * n, n), mapSoln(xSoln + i * n, n);
      mapSoln = m_Solver.solve(mapRHS);
      }
  }

  ~EigenSolverInterfaceInternal()
  {
    if(m_SparseMatrix)
      delete m_SparseMatrix;
  }

private:
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor, TIndex> SparseType;
  typedef Eigen::Map<SparseType> SparseMap;
  SparseMap *m_SparseMatrix = nullptr;

  Eigen::BiCGSTAB<SparseType, Eigen::IdentityPreconditioner> m_Solver;
};


void 
EigenSolverInterface
::SymbolicFactorization(size_t n, int *idxRows, int *idxCols, double *xMatrix)
{
  if(m_InternalSolver_Int)
    delete m_InternalSolver_Int;

  m_InternalSolver_Int = new EigenSolverInterfaceInternal<int>;
  m_InternalSolver_Int->SetMatrix(n, idxRows, idxCols, xMatrix);
}

void 
EigenSolverInterface
::SymbolicFactorization(const ImmutableSparseMatrix<double> &mat)
{
  if(m_InternalSolver_SizeT)
    delete m_InternalSolver_SizeT;

  m_InternalSolver_SizeT = new EigenSolverInterfaceInternal<long>;
  m_InternalSolver_SizeT->SetMatrix(mat.GetNumberOfRows(),
                                    reinterpret_cast<long *>(const_cast<size_t *>(mat.GetRowIndex())),
                                    reinterpret_cast<long *>(const_cast<size_t *>(mat.GetColIndex())),
                                    const_cast<double *>(mat.GetSparseData()));
}

void 
EigenSolverInterface
::NumericFactorization(const double *xMatrix)
{
  if(m_InternalSolver_Int)
    {
    m_InternalSolver_Int->UpdateMatrix(const_cast<double *>(xMatrix));
    m_InternalSolver_Int->Compute();
    }
  else if(m_InternalSolver_SizeT)
    {
    m_InternalSolver_SizeT->UpdateMatrix(const_cast<double *>(xMatrix));
    m_InternalSolver_SizeT->Compute();
    }
}

void 
EigenSolverInterface
::Solve(double *xRhs, double *xSoln)
{
  Solve(1, xRhs, xSoln);
}

void                              
EigenSolverInterface
::Solve(size_t nRHS, double *xRhs, double *xSoln)
{
  if(m_InternalSolver_Int)
    {
    m_InternalSolver_Int->Solve(nRHS, xRhs, xSoln);
    }
  else if(m_InternalSolver_SizeT)
    {
    m_InternalSolver_SizeT->Solve(nRHS, xRhs, xSoln);
    }
}

EigenSolverInterface::~EigenSolverInterface()
{
  if(m_InternalSolver_Int)
    delete m_InternalSolver_Int;
  if(m_InternalSolver_SizeT)
    delete m_InternalSolver_SizeT;
}
