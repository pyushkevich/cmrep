#ifndef __TestSolver_h_
#define __TestSolver_h_

#include "GenericMedialModel.h"
#include "CoefficientMapping.h"

class MedialOptimizationProblem;

// This method tests gradient computation. The third argument is the starting
// point in parameter space. By default it's at the origin.
int TestGradientComputation(
  GenericMedialModel *xSolver, 
  CoefficientMapping *xMask, 
  vnl_vector<double> P0,
  int nRandVariations = 0);

inline int TestGradientComputation(
  GenericMedialModel *xSolver, 
  CoefficientMapping *xMask, 
  int nRandVariations = 0)
{
  return TestGradientComputation(xSolver, xMask, 
    vnl_vector<double>(xMask->GetNumberOfParameters(),0.0), nRandVariations);
}

int TestOptimizerGradientComputation(
  MedialOptimizationProblem &mop, 
  CoefficientMapping &xMapping,
  GenericMedialModel *xSolver,
  double xStepSize, const char *nm_term, const char *nm_map);

#endif
