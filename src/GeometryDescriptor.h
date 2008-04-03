#ifndef __GeometryDescriptor_h_
#define __GeometryDescriptor_h_

#include <iostream>

struct GeometryDescriptor
{
  double xCovariantTensor[2][2];
  double xContravariantTensor[2][2];
  double xChristoffelFirst[2][2][2];

  // xChristoffelSecond[i][j][k] => \Gamma_{ij}^k
  double xChristoffelSecond[2][2][2];

  // The determinant of the covariant tensor and its inverse
  double g, gInv;

  // Initialize the descriptor using first derivatives (sets metric tensor)
  void SetOneJet(double *X, double *Xu, double *Xv);

  // Initialize the descriptor using second derivatives (sets Christoffels)
  void SetTwoJet(double *X, double *Xu, double *Xv, double *Xuu, double *Xuv, double *Xvv);

  // Compute the squared gradient magnitude of some field 
  double SquaredGradientMagnitude(double Fu, double Fv)
    {
    // Use the least number of operations (5 mult, 3 add)
    double z = xContravariantTensor[0][1] + xContravariantTensor[0][1];
    return
      Fu * (xContravariantTensor[0][0] * Fu + z * Fv) + 
      Fv * Fv * xContravariantTensor[1][1]; 
    }

  // Dump information out
  void PrintSelf(std::ostream &str);
};

#endif
