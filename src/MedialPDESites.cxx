#include "MedialPDESites.h"
#include <iostream>

using namespace std;


double FDInternalSite::ComputeEquation(const Mat &Y)
{
  double Fu, Fv, Fuu, Fuv, Fvv;
  
  // Compute the partial derivatives of F with respect to u and v
  xMask->ComputeTwoJet(Y, Fu, Fv, Fuu, Fuv, Fvv);

  // Return the sum weighted by geometry
  return Cuu * Fuu + Cuv * Fuv + Cvv * Fvv + Cu * Fu + Cv * Fv - rho;
};

/** Compute the derivatives for Newton's method */
void FDInternalSite
::ComputeDerivative(const Mat &, double *A, unsigned int iRow)
{
  // Because the equation is linear, there is no need to use values of 
  // Phi. We just copy the weightings
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivative[i];
}

void FDInternalSite::SetGeometry(GeometryDescriptor *g, double rho)
{
  // Store the rho
  this->rho = rho;

  // Compute the factors applied to partials of Phi at the site
  Cuu = g->xContravariantTensor[0][0];
  Cvv = g->xContravariantTensor[1][1];
  Cuv = g->xContravariantTensor[1][0] * 2.0;
  Cu = - (
    g->xContravariantTensor[0][0] * g->xChristoffelSecond[0][0][0] + 
    g->xContravariantTensor[0][1] * g->xChristoffelSecond[0][1][0] * 2.0 +
    g->xContravariantTensor[1][1] * g->xChristoffelSecond[1][1][0]);
  Cv = - (
    g->xContravariantTensor[0][0] * g->xChristoffelSecond[0][0][1] + 
    g->xContravariantTensor[0][1] * g->xChristoffelSecond[0][1][1] * 2.0 +
    g->xContravariantTensor[1][1] * g->xChristoffelSecond[1][1][1]);

  // Compute the derivative with respect to each cite
  for(size_t i = 0; i < xMask->Size(); i++)
    {
    // Get the contribution of node i to each partial derivative
    double Wu, Wv, Wuu, Wuv, Wvv;
    xMask->GetNodeComponentInTwoJet(i, Wu, Wv, Wuu, Wuv, Wvv);
    
    // Add them with weighting to produce the derivatives
    xDerivative[i] = Cu * Wu + Cv * Wv + Cuu * Wuu + Cuv * Wuv + Cvv * Wvv;
    }
}

/**
 * Compute the matrix component of the variational derivative
 */
void 
FDInternalSite::
ComputeVariationalDerivativeMatrix(
  const Mat &Y, double *A, MedialAtom *xAtom)
{
  // The A[i] elements are computed in the same way as for Newton's method
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivative[i];

  // Compute the partials of phi
  xMask->ComputeTwoJet(Y, Fu, Fv, Fuu, Fuv, Fvv);
}

/** 
 * Compute the variational derivative 
 */
double
FDInternalSite::
ComputeVariationalDerivativeRHS(
  const Mat &Y, MedialAtom *xAtom, MedialAtom *dAtom)
{
  // First get shorthands for the relevant vectors
  SMLVec3d &Nu = dAtom->Xu; SMLVec3d &Xu = xAtom->Xu;
  SMLVec3d &Nv = dAtom->Xv; SMLVec3d &Xv = xAtom->Xv;
  SMLVec3d &Nuu = dAtom->Xuu; SMLVec3d &Xuu = xAtom->Xuu;
  SMLVec3d &Nuv = dAtom->Xuv; SMLVec3d &Xuv = xAtom->Xuv;
  SMLVec3d &Nvv = dAtom->Xvv; SMLVec3d &Xvv = xAtom->Xvv;

  // Get shorthand for the differential geometric operators
  double (*G2)[2] = xAtom->G.xContravariantTensor;
  double (*K)[2][2] = xAtom->G.xChristoffelFirst;
  double (*K2)[2][2] = xAtom->G.xChristoffelSecond;
  
  // The derivative of the contravariant tensor
  double (*z)[2] = dAtom->G.xContravariantTensor; 

  // Compute the derivative of the Christoffel symbols of the second kind
  double NuXuu = dot(Nu, Xuu), NvXuu = dot(Nv, Xuu);
  double NuXuv = dot(Nu, Xuv), NvXuv = dot(Nv, Xuv);
  double NuXvv = dot(Nu, Xvv), NvXvv = dot(Nv, Xvv);
  double NuuXu = dot(Nuu, Xu), NuuXv = dot(Nuu, Xv);
  double NuvXu = dot(Nuv, Xu), NuvXv = dot(Nuv, Xv);
  double NvvXu = dot(Nvv, Xu), NvvXv = dot(Nvv, Xv);
  
  // The derivative of Christofel symbols of second kind  
  double (*Q)[2][2] = dAtom->G.xChristoffelSecond;
  Q[0][0][0] = 
    z[0][0] * K[0][0][0] + G2[0][0] * ( NuXuu + NuuXu ) +
    z[0][1] * K[0][0][1] + G2[0][1] * ( NvXuu + NuuXv );
  Q[0][0][1] = 
    z[1][0] * K[0][0][0] + G2[1][0] * ( NuXuu + NuuXu ) +
    z[1][1] * K[0][0][1] + G2[1][1] * ( NvXuu + NuuXv );
  
  Q[0][1][0] = Q[1][0][0] = 
    z[0][0] * K[1][0][0] + G2[0][0] * ( NuXuv + NuvXu ) +
    z[0][1] * K[1][0][1] + G2[0][1] * ( NvXuv + NuvXv );
  Q[0][1][1] = Q[1][0][1] = 
    z[1][0] * K[1][0][0] + G2[1][0] * ( NuXuv + NuvXu ) +
    z[1][1] * K[1][0][1] + G2[1][1] * ( NvXuv + NuvXv );

  Q[1][1][0] = 
    z[0][0] * K[1][1][0] + G2[0][0] * ( NuXvv + NvvXu ) +
    z[0][1] * K[1][1][1] + G2[0][1] * ( NvXvv + NvvXv );
  Q[1][1][1] = 
    z[1][0] * K[1][1][0] + G2[1][0] * ( NuXvv + NvvXu ) +
    z[1][1] * K[1][1][1] + G2[1][1] * ( NvXvv + NvvXv );

  // Compute the derivative (this should be the derivative)
  return dAtom->xLapR - ( 
    z[0][0] * (Fuu - K2[0][0][0] * Fu - K2[0][0][1] * Fv) - 
    G2[0][0] * (Q[0][0][0] * Fu + Q[0][0][1] * Fv) +
    z[1][1] * (Fvv - K2[1][1][0] * Fu - K2[1][1][1] * Fv) - 
    G2[1][1] * (Q[1][1][0] * Fu + Q[1][1][1] * Fv) +
    2.0 * ( 
      z[0][1] * (Fuv - K2[0][1][0] * Fu - K2[0][1][1] * Fv) - 
      G2[0][1] * (Q[0][1][0] * Fu + Q[0][1][1] * Fv) ));
}


/** 
 * Compute the variational derivative with respect to N
 */


void 
FDInternalSite::PrintReport()
{
	cout << "Internal Site at " << xMask->GetLocationU() 
		<< "," << xMask->GetLocationV() << endl;
	cout << "LB(F) = " 
		<< Cuu << " * Fuu + " 
		<< Cuv << " * Fuv + "
		<< Cvv << " * Fvv + "
		<< Cu  << " * Fu + "
		<< Cv  << " * Fv " << endl;
}

double FDBorderSite::ComputeEquation(const Mat &Y)
{
  double F, Fu, Fv;
  
  // Compute the partial derivatives of F with respect to u and v
  F = xMask->ComputeOneJet(Y.data_block(), Fu, Fv);

  // Compute the generalized gradient magnitude using precomputed values
  return CuCu * Fu * Fu + CuCv * Fu * Fv + CvCv * Fv * Fv - 4.0 * F;
}

/** Compute the derivatives for Newton's method */
void FDBorderSite
::ComputeDerivative(const Mat &Y, double *A, unsigned int iRow)
{
  double F, Fu, Fv;
  
  // Compute the partial derivatives of F with respect to u and v
  F = xMask->ComputeOneJet(Y.data_block(), Fu, Fv);

  // Because the equation is linear, there is no need to use values of 
  // Phi. We just copy the weightings
  for(size_t i = 0; i < xMask->Size(); i++)
    A[i] = xDerivativeFu[i] * Fu + xDerivativeFv[i] * Fv + xDerivativeF[i];
}

void
FDBorderSite
::ComputeVariationalDerivativeMatrix(
  const Mat &Y, double *A, MedialAtom *xAtom)
{
  // First, let's compute grad phi
  F = xMask->ComputeOneJet(Y.data_block(), Fu, Fv);

  // Next, compute the weights for Hu and Hv and H
  double (*G1)[2] = xAtom->G.xCovariantTensor;
  double (*G2)[2] = xAtom->G.xContravariantTensor;

  // Set the entries in the derivative matrix. Should be the same as
  // in the usual case...
  for(size_t i = 0; i < xMask->Size(); i++)
    {
    // Get the contribution of node i to the jet of H
    double Wu, Wv, W;
    W = xMask->GetNodeComponentInOneJet(i, Wu, Wv);
    
    // Compute the equation
    A[i] = -4.0 * W + 
      2.0 * (G2[0][0] * Wu + G2[0][1] * Wv) * Fu + 
      2.0 * (G2[1][0] * Wu + G2[1][1] * Wv) * Fv;
    }
}


double
FDBorderSite
::ComputeVariationalDerivativeRHS(
  const Mat &Y, MedialAtom *xAtom, MedialAtom *dAtom)
{
  // The derivative of the contravariant tensor
  double (*z)[2] = dAtom->G.xContravariantTensor; 

  // Compute the right hand side
  return - (z[0][0] * Fu * Fu + 2.0 * z[0][1] * Fu * Fv + z[1][1] * Fv * Fv);
}


/** Initialize the border site */
void FDBorderSite::SetGeometry(GeometryDescriptor *g, double)
{
  // Compute the weights of the terms in the equation
  CuCu = g->xContravariantTensor[0][0];
  CuCv = g->xContravariantTensor[0][1] * 2.0;
  CvCv = g->xContravariantTensor[1][1];

  // Compute the derivatives of Fu and Fv with respect to each node
  for(size_t i = 0; i < xMask->Size(); i++)
    {
    double W, Wu, Wv;
    
    // Get the contribution of node i to the partials 
    W = xMask->GetNodeComponentInOneJet(i, Wu, Wv);
    
    // The derivative of Fu with respect to neighbor i
    xDerivativeFu[i] = 2.0 * CuCu * Wu + CuCv * Wv;
      
    // The derivative of Fv with respect to neighbor i
    xDerivativeFv[i] = 2.0 * CvCv * Wv + CuCv * Wu;
      
    // The derivative of F with respect to neighbor i
    xDerivativeF[i] = -4.0 * W;
    }
}

void 
FDBorderSite::PrintReport()
{
	cout << "Border Site at " << xMask->GetLocationU() 
		<< "," << xMask->GetLocationV() << endl;
	cout << "LB(F) = " 
		<< CuCu << " * Fu * Fu + " 
		<< CuCv << " * Fu * Fv + "
		<< CvCv << " * Fv * Fv + " << endl;
}

