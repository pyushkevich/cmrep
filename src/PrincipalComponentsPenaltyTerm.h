#ifndef __PrincipalComponentsPenaltyTerm_h_
#define __PrincipalComponentsPenaltyTerm_h_

#include "OptimizationTerms.h"

class PrincipalComponentsPenaltyTerm : public EnergyTerm
{
public:
  
  PrincipalComponentsPenaltyTerm(GenericMedialModel *model);
  ~PrincipalComponentsPenaltyTerm();

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
    { return "PCAPEN"; }

  // Pass in parameters using a registry object
  void SetParameters(Registry &r);

private:
  
  double UnifiedComputeEnergy(SolutionData *S, bool gradient_mode);
  
};


#endif
