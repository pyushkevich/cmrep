#ifndef __JacobianDistortionPenaltyTerm_h_
#define __JacobianDistortionPenaltyTerm_h_

#include "OptimizationTerms.h"
#include "ScriptInterface.h"
#include "MeshTraversal.h"


class MedialJacobianDistortionPenaltyTerm : public EnergyTerm
{
public:
  
  MedialJacobianDistortionPenaltyTerm(GenericMedialModel *model);
  ~MedialJacobianDistortionPenaltyTerm();

  // Compute the energy
  double ComputeEnergy(SolutionData *S)
    { return UnifiedComputeEnergy(S, false); }
  
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
    { return "MEDREG"; }

  // Pass in parameters using a registry object
  void SetParameters(Registry &r);

private:

  // Reference medial PDE
  medialpde::MedialPDE *mpdeRef;

  // Pointer to the model and the reference model
  GenericMedialModel *xModel, *xRefModel;

  // Solution data object for the reference model
  SolutionData *sdRef;

  // This object will be used to compute the gradient on the 
  // surface of the reference model
  MeshGradientComputer xGradComp;

  // Precomputable data arrays
  typedef vnl_vector_fixed<double, 3> Vec;
  vnl_vector<double> xLogRefAE, xLogAERatio, dLogAERatio;
  Vec *xRefMed;

  StatisticsAccumulator saGradient;


  double UnifiedComputeEnergy(SolutionData *S, bool flagGradient);


  

};




class BoundaryJacobianDistortionPenaltyTerm : public EnergyTerm
{
public:
  
  BoundaryJacobianDistortionPenaltyTerm(GenericMedialModel *model);
  ~BoundaryJacobianDistortionPenaltyTerm();

  // Compute the energy
  double ComputeEnergy(SolutionData *S)
    { return UnifiedComputeEnergy(S, false); }
  
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
    { return "BNDREG"; }

  // Pass in parameters using a registry object
  void SetParameters(Registry &r);

private:

  // Reference medial PDE
  medialpde::MedialPDE *mpdeRef;

  // Pointer to the model and the reference model
  GenericMedialModel *xModel, *xRefModel;

  // Solution data object for the reference model
  SolutionData *sdRef;

  // This object will be used to compute the gradient on the 
  // surface of the reference model
  MeshGradientComputer xGradComp;

  // Precomputable data arrays
  typedef vnl_vector_fixed<double, 3> Vec;
  vnl_vector<double> xLogRefAE, xLogAERatio, dLogAERatio;
  Vec *xRefBnd;

  StatisticsAccumulator saGradient;


  double UnifiedComputeEnergy(SolutionData *S, bool flagGradient);


  

};

#endif
