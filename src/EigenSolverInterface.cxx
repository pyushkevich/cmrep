#include "EigenSolverInterface.h"
#include <iostream>
#include <cstdlib>
#include <cstring>

#include <Eigen/Sparse>

#ifdef HAVE_MKL
#include <Eigen/PardisoSupport>
#endif

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

// TODO: handle symmetric matrices!
#ifdef HAVE_MKL
  Eigen::PardisoLU<SparseType> m_Solver;
#else
  Eigen::BiCGSTAB<SparseType, Eigen::IdentityPreconditioner> m_Solver;
#endif
};

EigenSolverInterface
::EigenSolverInterface(ProblemType ptype)
: m_Type(ptype)
{
  m_InternalSolver = new EigenSolverInterfaceInternal<int>;
}

void
EigenSolverInterface
::ResetIndices()
{
  if(m_RowIndex) 
    {
    delete [] m_RowIndex;
    m_RowIndex = nullptr;
    }
  if(m_ColIndex)
    {
    delete [] m_RowIndex;
    m_ColIndex = nullptr;
    }
}

void 
EigenSolverInterface
::SymbolicFactorization(size_t n, int *idxRows, int *idxCols, double *xMatrix)
{
  ResetIndices();
  m_InternalSolver->SetMatrix(n, idxRows, idxCols, xMatrix);
}

void 
EigenSolverInterface
::SymbolicFactorization(const ImmutableSparseMatrix<double> &mat)
{
  ResetIndices();
  int nr = mat.GetNumberOfRows() + 1;
  int nc = mat.GetNumberOfSparseValues();
  m_RowIndex = new int[nr];
  m_ColIndex = new int[nc];

  for(unsigned int r = 0; r < nr; r++)
    m_RowIndex[r] = (int) mat.GetRowIndex()[r];
  for(unsigned int c = 0; c < nc; c++)
    m_ColIndex[c] = (int) mat.GetColIndex()[c];

  m_InternalSolver->SetMatrix(mat.GetNumberOfRows(),
      m_RowIndex, m_ColIndex,
      const_cast<double *>(mat.GetSparseData()));
}

void 
EigenSolverInterface
::NumericFactorization(const double *xMatrix)
{
  m_InternalSolver->UpdateMatrix(const_cast<double *>(xMatrix));
  m_InternalSolver->Compute();
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
  m_InternalSolver->Solve(nRHS, xRhs, xSoln);
}

EigenSolverInterface::~EigenSolverInterface()
{
  ResetIndices();
  delete m_InternalSolver;
}
