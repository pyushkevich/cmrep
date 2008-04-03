#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

#include "BasisFunctions2D.h"
#include "BasisFunctions2D.txx"

double CosineBasisFunction::Evaluate(double u, size_t k, size_t d)
{
  // Shift u so we don't have zero derivatives
  double b = 0.23525435234532, a = 1.0 - 2.0 * b;

  // Compute the basis at the point
  double T1 = vnl_math::pi * k * a, T2 = vnl_math::pi * k * (a * u + b);

  switch(d) 
    {
    case 0: return cos(T2);
    case 1: return - T1 * sin(T2);
    case 2: return - T1 * T1 * cos(T2);
    default: return 0.0;
    }
}

FourierSurface::FourierSurface(size_t ncu, size_t ncv) : 
  FourierSurfaceBase(ncu, ncv),
  xC2FDescriptor(ncu, ncv)
{
}


FourierSurface::FourierSurface(const FourierSurface &source) : 
  FourierSurfaceBase(source.ncu, source.ncv), 
  xC2FDescriptor(source.ncu, source.ncv)
{
  // Copy the raw coefficients
  this->SetCoefficientArray(source.GetCoefficientArray());
}

/**
 * Get a coefficient mask with given numbers of coefficients in X and Rho
 */
vector<size_t>  
FourierSurface
::GetCoefficientSubset( size_t ncuX, size_t ncvX, size_t ncuRho, size_t ncvRho)
{
  size_t i, j;

  // Construct a vector of coefficients that will be included
  vector<size_t> iSelect;

  // Add the X components
  for(i = 0; i < ncuX; i++) for(j = 0; j < ncvX; j++)
    {
    iSelect.push_back(GetCoefficientIndex(i, j, 0));
    iSelect.push_back(GetCoefficientIndex(i, j, 1));
    iSelect.push_back(GetCoefficientIndex(i, j, 2));
    }

  // Add the Rho components
  for(i = 0; i < ncuRho; i++) for(j = 0; j < ncvRho; j++)
    iSelect.push_back(GetCoefficientIndex(i, j, 3));

  // Return the coefficient vector
  return iSelect;
}


IHyperSurface2D *
FourierSurface
::GetVariationSurface(const double *xCoeff)
{
  // Create a new surface with the same parameters as this one
  FourierSurface *xVariation = new FourierSurface(ncu, ncv);

  // Set the components of the variation to the components passed in; 
  // we can do this because the function computed here is linear in the 
  // coefficients
  xVariation->SetCoefficientArray(xCoeff);

  // Copy the evaluation grid
  VectorType uu, vv;
  GetEvaluationGrid(uu, vv);
  xVariation->SetEvaluationGrid(uu, vv);

  // Return this variation  
  return xVariation;
}

void
FourierSurface
::ReleaseVariationSurface(IHyperSurface2D *xSurface)
{
  delete xSurface;
}







template class GenericBasisRepresentation2D< 4, 3,        
         CosineBasisFunction, CosineBasisFunction>;
