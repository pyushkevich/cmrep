#include "GeometryDescriptor.h"

using namespace std;

void
GeometryDescriptor
::SetOneJet(double *X, double *Xu, double *Xv)
{
  // Compute the covariant tensor
  xCovariantTensor[0][0] = Xu[0] * Xu[0] + Xu[1] * Xu[1] + Xu[2] * Xu[2];
  xCovariantTensor[1][1] = Xv[0] * Xv[0] + Xv[1] * Xv[1] + Xv[2] * Xv[2];
  xCovariantTensor[1][0] = Xu[0] * Xv[0] + Xu[1] * Xv[1] + Xu[2] * Xv[2];
  xCovariantTensor[0][1] = xCovariantTensor[1][0];

  // Compute the determinant of the covariant tensor
  g = xCovariantTensor[0][0] * xCovariantTensor[1][1] 
    - xCovariantTensor[0][1] * xCovariantTensor[0][1];
  gInv = 1.0 / g;

  // Compute the contravariant tensor
  xContravariantTensor[0][0] = gInv * xCovariantTensor[1][1];
  xContravariantTensor[1][1] = gInv * xCovariantTensor[0][0];
  xContravariantTensor[0][1] = - gInv * xCovariantTensor[1][0];
  xContravariantTensor[1][0] = xContravariantTensor[0][1];
}

void
GeometryDescriptor
::SetTwoJet(double *X, double *Xu, double *Xv, double *Xuu, double *Xuv, double *Xvv)
{
  // Set the 1-jet
  this->SetOneJet(X, Xu, Xv);

  // Compute the Christoffel symbols of the first kind
  xChristoffelFirst[0][0][0] = Xuu[0] * Xu[0] + Xuu[1] * Xu[1] + Xuu[2] * Xu[2];
  xChristoffelFirst[0][0][1] = Xuu[0] * Xv[0] + Xuu[1] * Xv[1] + Xuu[2] * Xv[2];
  xChristoffelFirst[0][1][0] = Xuv[0] * Xu[0] + Xuv[1] * Xu[1] + Xuv[2] * Xu[2];
  xChristoffelFirst[0][1][1] = Xuv[0] * Xv[0] + Xuv[1] * Xv[1] + Xuv[2] * Xv[2];
  xChristoffelFirst[1][1][0] = Xvv[0] * Xu[0] + Xvv[1] * Xu[1] + Xvv[2] * Xu[2];
  xChristoffelFirst[1][1][1] = Xvv[0] * Xv[0] + Xvv[1] * Xv[1] + Xvv[2] * Xv[2];
  xChristoffelFirst[1][0][0] = xChristoffelFirst[0][1][0];
  xChristoffelFirst[1][0][1] = xChristoffelFirst[0][1][1];

  // Compute the Christoffel symbols of the second kind
  size_t i, j, k;
  for(i = 0; i < 2; i++) for(j = 0; j < 2; j++) for(k = 0; k < 2; k++)
  	xChristoffelSecond[i][j][k] = 
	  	xContravariantTensor[k][0] * xChristoffelFirst[i][j][0] + 
	  	xContravariantTensor[k][1] * xChristoffelFirst[i][j][1]; 	  	 
}

void 
GeometryDescriptor
::PrintSelf(ostream &str)
{
  str << "CovariantTensor : {{" 
    << xCovariantTensor[0][0] << "," 
    << xCovariantTensor[0][1] << "}, {"
    << xCovariantTensor[1][0] << "," 
    << xCovariantTensor[1][1] << "}}" << endl;
  
  str << "ContravariantTensor : {{" 
    << xContravariantTensor[0][0] << "," 
    << xContravariantTensor[0][1] << "}, {"
    << xContravariantTensor[1][0] << "," 
    << xContravariantTensor[1][1] << "}}" << endl;

  str << "ChristoffelFirst : {{{" 
    << xChristoffelFirst[0][0][0] << ","
    << xChristoffelFirst[0][0][1] << "}, {"
    << xChristoffelFirst[0][1][0] << ","
    << xChristoffelFirst[0][1][1] << "}}, {{"
    << xChristoffelFirst[1][0][0] << ","
    << xChristoffelFirst[1][0][1] << "}, {"
    << xChristoffelFirst[1][1][0] << ","
    << xChristoffelFirst[1][1][1] << "}}}" << endl;

  str << "ChristoffelSecond : {{{" 
    << xChristoffelSecond[0][0][0] << ","
    << xChristoffelSecond[0][0][1] << "}, {"
    << xChristoffelSecond[0][1][0] << ","
    << xChristoffelSecond[0][1][1] << "}}, {{"
    << xChristoffelSecond[1][0][0] << ","
    << xChristoffelSecond[1][0][1] << "}, {"
    << xChristoffelSecond[1][1][0] << ","
    << xChristoffelSecond[1][1][1] << "}}}" << endl;
}

