#include "JacobianDistortionPenaltyTerm.h"
#include "ScriptInterface.h"

#include "JacobianDistortionPenaltyTerm.h"
#include "ScriptInterface.h"

MedialJacobianDistortionPenaltyTerm
::MedialJacobianDistortionPenaltyTerm(GenericMedialModel *model)
{
  // Get the medial model
  this->xModel = model;
  this->mpdeRef = NULL;
  this->xRefModel = NULL;
  this->sdRef = NULL;

  // Initialize precomputable arrays
  xLogRefAE.set_size(xModel->GetNumberOfAtoms());
  xLogAERatio.set_size(xModel->GetNumberOfAtoms());
  dLogAERatio.set_size(xModel->GetNumberOfAtoms());
  this->xRefMed = new Vec[xModel->GetNumberOfAtoms()];
}

MedialJacobianDistortionPenaltyTerm
::~MedialJacobianDistortionPenaltyTerm()
{
  if(mpdeRef)
    delete mpdeRef;

  if(sdRef)
    delete sdRef;

  delete xRefMed;
}

void
MedialJacobianDistortionPenaltyTerm
::SetParameters(Registry &r)
{
  // Get the filename of the reference model
  string fnref = r["ReferenceModel"][""];

  // If there is no reference model, crap out
  if(fnref == "") 
    throw string("Missing ReferenceModel in MedialJacobianDistortionPenaltyTerm");

  // If the model is loaded, clear it
  if(mpdeRef)
    delete mpdeRef;

  // Load the reference model
  mpdeRef = new medialpde::MedialPDE(fnref.c_str());

  // Check that the reference model matches
  xRefModel = mpdeRef->GetMedialModel();
  if(
    (xModel->GetNumberOfAtoms() != xRefModel->GetNumberOfAtoms()) || 
    (xModel->GetNumberOfTriangles() != xRefModel->GetNumberOfTriangles()))
    {
    throw string("Mismatch between reference model and current model");
    }

  // Compute the area element along the boundary of the model
  if(sdRef) delete sdRef;
  sdRef = new SolutionData(
    xRefModel->GetIterationContext(), xRefModel->GetAtomArray());
  sdRef->ComputeIntegrationWeights();

  // Initialize the gradient computer
  xGradComp.SetMesh(xRefModel->GetIterationContext()->GetMedialMesh());

  // Fill out the precomputable values 
  for(MedialAtomIterator it(xRefModel->GetIterationContext()); 
    !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    xLogRefAE[i] = log(sdRef->xMedialWeights[i]);
    xRefMed[i] = xRefModel->GetAtomArray()[i].X;
    }
}

double 
MedialJacobianDistortionPenaltyTerm
::UnifiedComputeEnergy(SolutionData *S, bool flagGradient)
{
  // Reset the accumulator
  saGradient.Reset();

  // We need the integration weights
  S->ComputeIntegrationWeights();

  // Evaluate log(A / Aref) at each point in the model
  for(size_t i = 0; i < xRefModel->GetNumberOfAtoms(); i++)
    xLogAERatio[i] = log(S->xMedialWeights[i]) - xLogRefAE[i];

  // Compute the gradient using the gradient computer
  xGradComp.ComputeGradient(
    xRefMed, xLogAERatio.data_block(), flagGradient);

  // Integrate the gradient magnitude over the reference model
  double xIntegral = 0.0;
  for(size_t i = 0; i < xRefModel->GetNumberOfAtoms(); i++)
    {
    MeshGradientComputer::Vec grd = xGradComp.gradF[i];
    double grdmagsq = grd.squared_magnitude();
    saGradient.Update(grdmagsq);
    xIntegral += sdRef->xMedialWeights[i] * grdmagsq;
    }

  return xIntegral / sdRef->xMedialArea;  
}

double
MedialJacobianDistortionPenaltyTerm
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // We will need boundary weights on dS
  dS->ComputeIntegrationWeights();

  // Evaluate log(A / Aref) at each point in the model
  for(size_t i = 0; i < xRefModel->GetNumberOfAtoms(); i++)
    dLogAERatio[i] = dS->xMedialWeights[i] / S->xMedialWeights[i];

  // Set the integral accumulator
  double dIntegral = 0.0;

  // For each point, now multiply by the Jacobian
  for(MedialAtomIterator it(S->xAtomGrid); !it.IsAtEnd(); ++it)
    {   
    // Ignore atoms that are not affected by the variation
    size_t i = it.GetIndex();
    if(S->xAtoms[i].order > 1)
      continue;

    // Compute the derivative of the gradient with respect to var-n
    Vec dGradF(0.0);
    for(MeshGradientComputer::SparseMatF::RowIterator rit = 
      xGradComp.jacF.Row(i); !rit.IsAtEnd(); ++rit)
      {
      dGradF += rit.Value() * dLogAERatio[rit.Column()];
      }

    // Accumulate
    dIntegral += 
      2 * sdRef->xMedialWeights[i] * dot_product(dGradF, xGradComp.gradF[i]);
    }

  return dIntegral / sdRef->xMedialArea;
}

void
MedialJacobianDistortionPenaltyTerm
::PrintReport(ostream &sout)
{
  sout << "  MedialJacobianDistortionPenaltyTerm" << endl;
  sout << "    Distortion Gradient Mean    : " << saGradient.GetMean() << endl;
  sout << "    Distortion Gradient Max     : " << saGradient.GetMax() << endl;
}







BoundaryJacobianDistortionPenaltyTerm
::BoundaryJacobianDistortionPenaltyTerm(GenericMedialModel *model)
{
  // Get the medial model
  this->xModel = model;
  this->mpdeRef = NULL;
  this->xRefModel = NULL;
  this->sdRef = NULL;

  // Initialize precomputable arrays
  xLogRefAE.set_size(xModel->GetNumberOfBoundaryPoints());
  xLogAERatio.set_size(xModel->GetNumberOfBoundaryPoints());
  dLogAERatio.set_size(xModel->GetNumberOfBoundaryPoints());
  this->xRefBnd = new Vec[xModel->GetNumberOfBoundaryPoints()];
}

BoundaryJacobianDistortionPenaltyTerm
::~BoundaryJacobianDistortionPenaltyTerm()
{
  if(mpdeRef)
    delete mpdeRef;

  if(sdRef)
    delete sdRef;

  delete xRefBnd;
}

void
BoundaryJacobianDistortionPenaltyTerm
::SetParameters(Registry &r)
{
  // Get the filename of the reference model
  string fnref = r["ReferenceModel"][""];

  // If there is no reference model, crap out
  if(fnref == "") 
    throw string("Missing ReferenceModel in BoundaryJacobianDistortionPenaltyTerm");

  // If the model is loaded, clear it
  if(mpdeRef)
    delete mpdeRef;

  // Load the reference model
  mpdeRef = new medialpde::MedialPDE(fnref.c_str());

  // Check that the reference model matches
  xRefModel = mpdeRef->GetMedialModel();
  if(
    (xModel->GetNumberOfAtoms() != xRefModel->GetNumberOfAtoms()) || 
    (xModel->GetNumberOfTriangles() != xRefModel->GetNumberOfTriangles()))
    {
    throw string("Mismatch between reference model and current model");
    }

  // Compute the area element along the boundary of the model
  if(sdRef) delete sdRef;
  sdRef = new SolutionData(
    xRefModel->GetIterationContext(), xRefModel->GetAtomArray());
  sdRef->ComputeIntegrationWeights();

  // Initialize the gradient computer
  xGradComp.SetMesh(xRefModel->GetIterationContext()->GetBoundaryMesh());

  // Fill out the precomputable values 
  for(MedialBoundaryPointIterator bit(xRefModel->GetIterationContext()); 
    !bit.IsAtEnd(); ++bit)
    {
    size_t i = bit.GetIndex();
    xLogRefAE[i] = log(sdRef->xBoundaryWeights[i]);
    xRefBnd[i] = GetBoundaryPoint(bit, xRefModel->GetAtomArray()).X;
    }
}

double 
BoundaryJacobianDistortionPenaltyTerm
::UnifiedComputeEnergy(SolutionData *S, bool flagGradient)
{
  // Reset the accumulator
  saGradient.Reset();

  // We need the integration weights
  S->ComputeIntegrationWeights();

  // Evaluate log(A / Aref) at each point in the model
  for(size_t i = 0; i < xRefModel->GetNumberOfBoundaryPoints(); i++)
    xLogAERatio[i] = log(S->xBoundaryWeights[i]) - xLogRefAE[i];

  // Compute the gradient using the gradient computer
  xGradComp.ComputeGradient(
    xRefBnd, xLogAERatio.data_block(), flagGradient);

  // Integrate the gradient magnitude over the reference model
  double xIntegral = 0.0;
  for(size_t i = 0; i < xRefModel->GetNumberOfBoundaryPoints(); i++)
    {
    MeshGradientComputer::Vec grd = xGradComp.gradF[i];
    double grdmagsq = grd.squared_magnitude();
    saGradient.Update(grdmagsq);
    xIntegral += sdRef->xBoundaryWeights[i] * grdmagsq;
    }

  return xIntegral / sdRef->xBoundaryArea;  
}

double
BoundaryJacobianDistortionPenaltyTerm
::ComputePartialDerivative(
  SolutionData *S, PartialDerivativeSolutionData *dS)
{
  // We will need boundary weights on dS
  dS->ComputeIntegrationWeights();

  // Evaluate log(A / Aref) at each point in the model
  for(size_t i = 0; i < xRefModel->GetNumberOfBoundaryPoints(); i++)
    dLogAERatio[i] = dS->xBoundaryWeights[i] / S->xBoundaryWeights[i];

  // Set the integral accumulator
  double dIntegral = 0.0;

  // For each point, now multiply by the Jacobian
  for(MedialBoundaryPointIterator bit(S->xAtomGrid); !bit.IsAtEnd(); ++bit)
    {   
    // Ignore atoms that are not affected by the variation
    if(S->xAtoms[bit.GetAtomIndex()].order > 2)
      continue;

    // Compute the derivative of the gradient with respect to var-n
    size_t i = bit.GetIndex();
    Vec dGradF(0.0);
    for(MeshGradientComputer::SparseMatF::RowIterator rit = 
      xGradComp.jacF.Row(i); !rit.IsAtEnd(); ++rit)
      {
      dGradF += rit.Value() * dLogAERatio[rit.Column()];
      }

    // Accumulate
    dIntegral += 
      2 * sdRef->xBoundaryWeights[i] * dot_product(dGradF, xGradComp.gradF[i]);
    }

  return dIntegral / sdRef->xBoundaryArea;
}

void
BoundaryJacobianDistortionPenaltyTerm
::PrintReport(ostream &sout)
{
  sout << "  BounaryJacobianDistortionPenaltyTerm" << endl;
  sout << "    Distortion Gradient Mean    : " << saGradient.GetMean() << endl;
  sout << "    Distortion Gradient Max     : " << saGradient.GetMax() << endl;
}
