#ifndef __DiffeomorphicEnergyTerm_h_
#define __DiffeomorphicEnergyTerm_h_

#include "OptimizationTerms.h"

/**
 * This penalty term is used to ensure that the transformation from
 * cm-rep parameter space to physical space is diffeomorphic. It uses
 * somewhat of a brute-force approach: test no boundary point in the 
 * model lies inside one of the spheres (M, R) on the medial axis. I 
 * don't currently have a formal proof, but I believe that if this 
 * condition is satisified, we can guarantee that the model is Blum.
 *
 * This is an O(n^2) algorithm, which is unfortunate, but the computation
 * involved in evaluation and gradient computation is reasonably fast, since
 * it does not involve any expensive operations, just some multiplications
 * and additions. 
 *
 * It also is possible to speed up the computation using a kd-tree or a
 * similar algorithm. However, it's not clear just how much of a boost
 * we would be able to get that way. Anyway, something we can figure out
 * later in the game
 */
class DiffeomorphicEnergyTerm : public EnergyTerm
{
public:
  
  DiffeomorphicEnergyTerm(GenericMedialModel *model);
  ~DiffeomorphicEnergyTerm();

  // Compute the energy term 
  double ComputeEnergy(SolutionData *data)
    { return UnifiedComputeEnergy(data, false); }

  // Initialize gradient computation and return the value of the solution
  // at the current state
  double BeginGradientComputation(SolutionData *SCenter)
    { return UnifiedComputeEnergy(SCenter, true); }

  // Compute the partial derivative of the energy function
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Finish gradient computation, remove all temporary data
  void EndGradientComputation() {};

  // Print a verbose report
  void PrintReport(ostream &sout);

  // Print a short name
  string GetShortName()
    { return "DIFFEO"; }

  // Pass in parameters using a registry object
  void SetParameters(Registry &r) {}

private:

  // Struct used to store information for every intersection
  typedef vnl_vector_fixed<double, 3> Vec;
  struct Intersection {
    // a, b, c such that the derivative of penalty is given by
    // Fij' = a . (Xi' - Xj') + b * R';
    Vec a;
    double b;
  };

  // Sparse matrix storage for proximity informatino
  typedef ImmutableSparseArray<Intersection> DistanceMatrix;
  typedef DistanceMatrix::STLSourceType MutableSparseMatrix;
  DistanceMatrix xDistMat;
  MutableSparseMatrix xMuteMat;

  // Pointer to the model
  GenericMedialModel *xModel;

  // Accumulators
  StatisticsAccumulator saPenalty;

  double UnifiedComputeEnergy(SolutionData *S, bool flagGradient);
};

#endif
