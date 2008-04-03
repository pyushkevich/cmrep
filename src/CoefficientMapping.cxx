#include "CoefficientMapping.h"
#include "PrincipalComponents.h"
#include "SubdivisionMedialModel.h"

PCACoefficientMapping
::PCACoefficientMapping(PrincipalComponents *pca, size_t nModes)
{ 
  this->pca = pca;
  this->n = nModes;
  this->m = pca->GetMean().size(); 
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::Apply(const Vec &C, const Vec &P)
{
  return C + pca->MapToFeatureSpaceZeroMean(P);
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP) 
{ 
  return pca->MapToFeatureSpaceZeroMean(varP); 
}

PCACoefficientMapping::Vec 
PCACoefficientMapping
::ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
{ 
  return varC; 
}

PCAPlusAffineCoefficientMapping
::PCAPlusAffineCoefficientMapping(
  GenericMedialModel *model, PrincipalComponents *pca, size_t nModes) :
  CompositionCoefficientMapping(
    new AffineTransformCoefficientMapping(model),
    new PCACoefficientMapping(pca, nModes))
{
}

PCAPlusAffineCoefficientMapping
::~PCAPlusAffineCoefficientMapping()
{
  delete this->f;
  delete this->g;
}
  

ReflectionCoefficientMapping
::ReflectionCoefficientMapping(
  SubdivisionMedialModel *model,
  SMLVec3d &p, double b)
{
  this->p = p;
  this->b = b;
  this->n = model->GetNumberOfCoefficients();
  this->m = model->GetNumberOfCoefficients();

  // Assuming 4 coefficients per node (standard for this code)
  this->k = m >> 2;

  // Resize opp array
  opp.resize(k);

  // Compute the opposites map. This means for each point (u,v) in the model
  // to find the point closest to (u,-v). We do this by brute force here
  for(size_t i = 0; i < k; i++)
    {
    double ui = model->GetCoefficientU()[i];
    double vi = model->GetCoefficientV()[i];
    double dist;
    size_t iopp;
    for(size_t j = 0; j < k; j++)
      {
      double uj = model->GetCoefficientU()[j];
      double vj = model->GetCoefficientV()[j];
      double d = (ui-uj)*(ui-uj) + (vi+vj)*(vi+vj);
      if(j == 0 || dist > d)
        { iopp = j; dist = d; }
      }
    opp[i] = iopp;
    // printf("Map (%f,%f) to (%f,%f)\n",ui,vi,
    //   model->GetCoefficientU()[opp[i]], model->GetCoefficientV()[opp[i]]);
    }
}

ReflectionCoefficientMapping::Vec
ReflectionCoefficientMapping
::Apply(const Vec &C, const Vec &P)
{
  Vec X = C+P, Y(n,0.0);
  for(size_t i = 0; i < k; i++)
    {
    // Get the vector component
    SMLVec3d Xi = X.extract(3,i<<2);
    SMLVec3d Xj = X.extract(3,opp[i] << 2);
    SMLVec3d Yi = 0.5 * (Xi + Xj) - (dot_product(p,Xj) - b) * p;
    double Ri = X[(i<<2)+3], Rj = X[(opp[i] << 2) + 3];
    double Qi = 0.5 * (Ri + Rj);
    
    // Set the output
    Y.update(Yi, i<<2);
    Y[(i<<2)+3] = Qi;
    }
  return Y;
}

/** Compute the variation J_T(P) * v_P given C, P and v_P */
ReflectionCoefficientMapping::Vec 
ReflectionCoefficientMapping
::ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
{
  Vec dY(n, 0.0);
  for(size_t i = 0; i < k; i++)
    {
    SMLVec3d Vi = varC.extract(3,i<<2);
    SMLVec3d Vj = varC.extract(3,opp[i]<<2);
    SMLVec3d dYi = 0.5 * (Vi + Vj) - dot_product(Vj,p) * p;
    double dQi = 0.5 * (varC[(i<<2)+3] + varC[(opp[i]<<2)+3]);
    dY.update(dYi,i<<2);
    dY[(i<<2)+3] = dQi;
    }
  return dY;
}

ReflectionCoefficientMapping::Vec 
ReflectionCoefficientMapping
::ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP)
{
  Vec dY(n, 0.0);
  for(size_t i = 0; i < k; i++)
    {
    SMLVec3d Vi = varP.extract(3,i<<2);
    SMLVec3d Vj = varP.extract(3,opp[i]<<2);
    SMLVec3d dYi = 0.5 * (Vi + Vj) - dot_product(Vj,p) * p;
    double dQi = 0.5 * (varP[(i<<2)+3] + varP[(opp[i]<<2)+3]);
    dY.update(dYi,i<<2);
    dY[(i<<2)+3] = dQi;
    }
  return dY;
}

