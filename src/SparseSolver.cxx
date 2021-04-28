#include "SparseSolver.h"
#include "MedialException.h"

#ifdef HAVE_PARDISO

#include "PardisoInterface.h"

SparseSolver* 
SparseSolver
::MakeSolver(bool symmetric)
{
  if(symmetric)
    return new SymmetricPositiveDefiniteRealPARDISO();
  else 
    return new UnsymmetricRealPARDISO();
}

#elif HAVE_TAUCS

#include "TaucsInterface.h"

SparseSolver* 
SparseSolver
::MakeSolver(bool symmetric)
{
  return new TaucsSolverInterface(symmetric);
}

#elif HAVE_EIGEN

#include "EigenSolverInterface.h"

SparseSolver* 
SparseSolver
::MakeSolver(bool symmetric)
{
  return new EigenSolverInterface(symmetric ? 
    EigenSolverInterface::SPD : 
    EigenSolverInterface::UNSYMMETRIC);
}

#else

SparseSolver* 
SparseSolver
::MakeSolver(bool symmetric)
{
  throw MedialModelException("The sparse solver has not been configured. Use PARDISO or TAUCS");
}

#endif

