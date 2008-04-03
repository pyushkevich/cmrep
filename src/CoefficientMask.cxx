#include "CoefficientMask.h"
#include "BasisFunctions2D.h"
#include "PrincipalComponents.h"

/******************************************************************************
 * SelectionMedialCoefficientMask 
 *****************************************************************************/

SelectionMedialCoefficientMask::SelectionMedialCoefficientMask(
  ICoefficientSettable *source, const vector<size_t> &mask)
: xMask(mask), 
  xRawCoefficients(mask.size()), xMaskCoefficients(mask.size())
{
  xSource = source;
  nCoeff = mask.size();
}

SelectionMedialCoefficientMask::SelectionMedialCoefficientMask(
  ICoefficientSettable *source, size_t iMaskSize, const size_t *iMask)
: xRawCoefficients(iMaskSize), xMaskCoefficients(iMaskSize)
{
  xSource = source;
  nCoeff = iMaskSize;
  xMask.reserve(iMaskSize);
  for(size_t i=0;i<iMaskSize;i++)
    xMask.push_back(iMask[i]);
}

/** Get the coefficient array */
double *SelectionMedialCoefficientMask::GetCoefficientArray()
{
  double *xRaw = xSource->GetCoefficientArray();
  for(size_t i=0; i < nCoeff; i++)
    { xMaskCoefficients[i] = xRaw[ xMask[i] ]; }
  return &(xMaskCoefficients[0]);
}

/** Set the coefficient array */
void SelectionMedialCoefficientMask::SetCoefficientArray(const double *xData)
{
  double *xRaw = xSource->GetCoefficientArray();
  for(size_t i=0; i < nCoeff; i++)
    xRaw[ xMask[i] ] = xData[i];
  xSource->SetCoefficientArray(xRaw);
}

double SelectionMedialCoefficientMask::GetCoefficient(size_t i) const
{
  return xSource->GetCoefficientArray()[xMask[i]];
}

void SelectionMedialCoefficientMask::SetCoefficient(size_t i, double x)
{ 
  xSource->SetCoefficient( xMask[i], x ); 
}

SelectionMedialCoefficientMask::Vec
SelectionMedialCoefficientMask
::GetVariationPerComponent(size_t iComponent)
{
  return xSource->GetVariationPerComponent(xMask[iComponent]);
}


IHyperSurface2D *
SelectionMedialCoefficientMask
::GetVariationSurface(const double *xData)
{
  // Create an array of coefficients, set unknowns to zeros
  vnl_vector<double> xFullData(xSource->GetNumberOfCoefficients(), 0.0);
  for(size_t i=0; i < nCoeff; i++)
    xFullData[ xMask[i] ] = xData[i];
  
  // Call the source's method
  return xSource->GetVariationPerVariation(xFullData.data_block());
}

/******************************************************************************
 * PassThroughCoefficientMask 
 *****************************************************************************/

PassThroughCoefficientMask::PassThroughCoefficientMask(
  ICoefficientSettable *source) 
: xSource(source) 
{
}

double *PassThroughCoefficientMask::GetCoefficientArray()
{ 
  return xSource->GetCoefficientArray(); 
}

void PassThroughCoefficientMask::SetCoefficientArray(const double *xData)
{ 
  xSource->SetCoefficientArray(xData); 
}

size_t PassThroughCoefficientMask::GetNumberOfCoefficients() const
{ 
  return xSource->GetNumberOfCoefficients(); 
}

double PassThroughCoefficientMask::GetCoefficient(size_t i) const
{ 
  return xSource->GetCoefficient(i); 
}

void PassThroughCoefficientMask::SetCoefficient(size_t i, double x)
{ 
  xSource->SetCoefficient(i, x); 
}

PassThroughCoefficientMask::Vec
PassThroughCoefficientMask
::GetVariationPerComponent(size_t iComponent)
{
  return xSource->GetVariationPerComponent(iComponent);
}

PassThroughCoefficientMask::Vec
PassThroughCoefficientMask
::GetVariationPerVariation(const double *xVariation)
{
  return xSource->GetVariationPerVariation(xVariation);
}

/******************************************************************************
 * AffineComponentSurface 
 *****************************************************************************/

class AffineComponentSurface : public IHyperSurface2D
{
public:
  // Typed definitions
  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  // Constructor: pass in matrix derivatives
  AffineComponentSurface(
    IHyperSurface2D *xSource, size_t nDim,
    const MatrixType &A, const VectorType &b, const VectorType &c)
    {
    this->A = A;
    this->b = b;
    this->c = c;
    this->xSource = xSource;
    this->nDim = nDim;
    }

  // Apply the affine transform (this is in the same way as in FourierSurface)
  void ApplyAffineTransform(VectorType &X, VectorType &Y, bool flagDeriv)
    {
    if(flagDeriv)
      { X[0] -= c[0]; X[1] -= c[1]; X[2] -= c[2]; }
    
    Y[0] = A[0][0] * X[0] + A[0][1] * X[1] + A[0][2] * X[2];
    Y[1] = A[1][0] * X[0] + A[1][1] * X[1] + A[1][2] * X[2];
    Y[2] = A[2][0] * X[0] + A[2][1] * X[1] + A[2][2] * X[2];
    Y[3] = X[3] * A[3][3];

    if(flagDeriv)
      { Y[0] += b[0]; Y[1] += b[1]; Y[2] += b[2]; }
    }

  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  void Evaluate(double u, double v, double *x)
    {
    VectorType X(nDim), Y(nDim);
    xSource->Evaluate(u, v, X.data_block());
    ApplyAffineTransform(X, Y, true);
    Y.copy_out(x);
    }
 
  /** Evaluate a partial derivative of the function at the u, v coordinate.
   * The first component is c0, number of components nc */
  void EvaluateDerivative(
    double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x)
    {
    VectorType X(nDim), Y(nDim);

    // Compute the surface point
    xSource->EvaluateDerivative(u, v, ou, ov, 0, nDim, X.data_block());
    
    // Apply the affine transform
    ApplyAffineTransform(X, Y, ou + ov == 0);

    // Copy the values into x
    for(size_t i = 0; i < nc; i++)
      x[i] =  Y[i + c0];
    }

  /** Evaluate the function at a grid point */
  void EvaluateAtGridIndex(
    size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x)
    {
    VectorType X(nDim), Y(nDim);

    // Compute the surface point
    xSource->EvaluateAtGridIndex(iu, iv, ou, ov, 0, nDim, X.data_block());
    
    // Apply the affine transform
    ApplyAffineTransform(X, Y, ou + ov == 0);

    // Copy the values into x
    for(size_t i = 0; i < nc; i++)
      x[i] =  Y[i + c0];
    }

  /** Get the evaluation grid parameters */
  void GetEvaluationGrid(VectorType &uu, VectorType &vv) const
    { return xSource->GetEvaluationGrid(uu, vv); }

  /** Get the number of dimensions of the surface */
  size_t GetNumberOfDimensions()
    { return xSource->GetNumberOfDimensions(); }

private:
  friend class AffineTransformCoefficientMask;
  friend class AffineAndPCACoefficientMask;
  size_t nDim;
  IHyperSurface2D *xSource;
  MatrixType A;
  VectorType b, c;
};

/******************************************************************************
 * AdditiveComponentSurface
 *****************************************************************************/

/**
 * This class represents a sum of two surfaces.
 */
class AdditiveComponentSurface : public IHyperSurface2D
{
public:
  
  // Constructor: pass in the surfaces
  AdditiveComponentSurface(IHyperSurface2D *xSurf1, IHyperSurface2D *xSurf2)
    {
    this->surf1 = xSurf1;
    this->surf2 = xSurf2;
    this->nDim = xSurf1->GetNumberOfDimensions();
    }

  /** Evaluate the function at a particular u, v coordinate. The return value
   * is a vector of doubles */
  void Evaluate(double u, double v, double *x)
    {
    vnl_vector<double> X1(nDim), X2(nDim), X;
    surf1->Evaluate(u, v, X1.data_block());
    surf2->Evaluate(u, v, X2.data_block());
    X = X1 + X2;
    X.copy_out(x);
    }
 
  /** Evaluate a partial derivative of the function at the u, v coordinate.
   * The first component is c0, number of components nc */
  void EvaluateDerivative(
    double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x)
    {
    vnl_vector<double> X1(nc), X2(nc), X;
    surf1->EvaluateDerivative(u, v, ou, ov, c0, nc, X1.data_block());
    surf2->EvaluateDerivative(u, v, ou, ov, c0, nc, X2.data_block());
    X = X1 + X2;
    X.copy_out(x);
    }

  /** Evaluate the function at a grid point */
  void EvaluateAtGridIndex(
    size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x)
    {
    vnl_vector<double> X1(nc), X2(nc), X;
    surf1->EvaluateAtGridIndex(iu, iv, ou, ov, c0, nc, X1.data_block());
    surf2->EvaluateAtGridIndex(iu, iv, ou, ov, c0, nc, X2.data_block());
    X = X1 + X2;
    X.copy_out(x);
    }

  /** Get the evaluation grid parameters */
  void GetEvaluationGrid(VectorType &uu, VectorType &vv) const
    { return surf1->GetEvaluationGrid(uu, vv); }

  /** Get the number of dimensions of the surface */
  size_t GetNumberOfDimensions()
    { return surf1->GetNumberOfDimensions(); }

private:
  IHyperSurface2D *surf1, *surf2;
  size_t nDim;

  friend class AffineAndPCACoefficientMask;
};


/******************************************************************************
 * AffineTransformCoefficientMask 
 *****************************************************************************/

AffineTransformCoefficientMask
::AffineTransformCoefficientMask(IBasisRepresentation2D *surface)
: xOriginalCoefficients(
    surface->GetCoefficientArray(), surface->GetNumberOfCoefficients()), 
  xSurface(surface), 
  nDim(surface->GetNumberOfDimensions()),
  A(nDim, nDim), b(nDim), xData(nDim * nDim + nDim)
{ 
  // Initialize the affine transform
  A.set_identity(); b.fill(0);
  c = surface->GetCenterOfRotation();

  // Copy it to the data vector
  A.copy_out(xData.data_block());

  iIndexB = nDim * nDim; 
  nCoeff = iIndexB + nDim;
  b.copy_out(xData.data_block() + iIndexB);
}

void AffineTransformCoefficientMask::SetCoefficientArray(const double *inData)
{ 
  xData.copy_in(inData);
  A.copy_in(inData);
  b.copy_in(inData + iIndexB);

  xSurface->SetCoefficientArray(GetCoefficientsBeforeTransform());
  xSurface->ApplyAffineTransform(A, b, c);
}

IHyperSurface2D *
AffineTransformCoefficientMask
::GetComponentSurface(size_t iCoefficient)
{
  // Create the derivative coefficient vector
  VectorType xIndex(nCoeff, 0.0);
  xIndex[iCoefficient] = 1.0;
  
  // Create the special surface
  IHyperSurface2D *xSurfaceRaw 
    = xSurface->GetVariationSurface(GetCoefficientsBeforeTransform());
  return CreateAffineVariation(xSurfaceRaw, xIndex.data_block());
}

AffineTransformCoefficientMask::Vec
AffineTransformCoefficientMask
::GetVariationPerComponent(size_t iComponent)
{
  // Create the derivative coefficient vector
  VectorType xIndex(nCoeff, 0.0);
  xIndex[iCoefficient] = 1.0;
  
  // Create the special surface
  IHyperSurface2D *xSurfaceRaw 
    = xSurface->GetVariationSurface(GetCoefficientsBeforeTransform());
  return CreateAffineVariation(xSurfaceRaw, xIndex.data_block());
}

IHyperSurface2D *
AffineTransformCoefficientMask
::GetVariationSurface(const double *xData)
{
  // Create the special surface
  IHyperSurface2D *xSurfaceRaw 
    = xSurface->GetVariationSurface(GetCoefficientsBeforeTransform());
  return CreateAffineVariation(xSurfaceRaw, xData);
}

IHyperSurface2D *
AffineTransformCoefficientMask
::CreateAffineVariation(
  IHyperSurface2D *xSurfaceVariation, const double *xAffineVariation)
{
  // Compute the matrices that are to be applied
  MatrixType dA(nDim, nDim);
  VectorType db(nDim);

  dA.copy_in(xAffineVariation);
  db.copy_in(xAffineVariation + iIndexB);

  return new AffineComponentSurface(xSurfaceVariation, nDim, dA, db, c);
}
  
void
AffineTransformCoefficientMask
::ReleaseComponentSurface(IHyperSurface2D *xCompSurface)
{
  ReleaseVariationSurface(xCompSurface);
}

void
AffineTransformCoefficientMask
::ReleaseVariationSurface(IHyperSurface2D *xCompSurface)
{
  AffineComponentSurface *xAffCompSurf 
    = dynamic_cast<AffineComponentSurface *>(xCompSurface);

  delete xAffCompSurf->xSource;
  delete xAffCompSurf;
}

/******************************************************************************
 * AffineTransform3DCoefficientMask 
 *****************************************************************************/
const size_t AffineTransform3DCoefficientMask::xIndexArray[] = 
  { 0, 1, 2, 4, 5, 6, 8, 9, 10, 16, 17, 18 };

/******************************************************************************
 * PCACoefficientMask 
 *****************************************************************************/
PCACoefficientMask
::PCACoefficientMask(
  const vnl_matrix<double> &mPCA, 
  ICoefficientSettable *xSurface, size_t nModes)
: xCoeffArray(nModes, 0.0), pca(mPCA)
{
  this->xSurface = xSurface;
  this->nModes = nModes;

  // Substract the PCA mean from the original coefficients vector
  xOriginalOffsetFromMean = vnl_vector<double>(
    xSurface->GetCoefficientArray(), xSurface->GetNumberOfCoefficients()) 
    - pca.GetMean();
}

double *PCACoefficientMask::GetCoefficientArray()
{
  return xCoeffArray.data_block();
}

void PCACoefficientMask::SetCoefficientArray(const double *xData)
{
  // Save the parameter vector
  xCoeffArray.copy_in(xData);
  
  // Apply the parameter transformation
  vnl_vector<double> X = 
    xOriginalOffsetFromMean + pca.MapToFeatureSpace(xCoeffArray);

  // Set the parameters of the parent surface
  xSurface->SetCoefficientArray(X.data_block());
}

double PCACoefficientMask::GetCoefficient(size_t i) const
{
  return xCoeffArray[i];
}

void PCACoefficientMask::SetCoefficient(size_t i, double x)
{
  // Update the coefficient array
  xCoeffArray[i] = x;

  // Apply the parameter transformation
  vnl_vector<double> X = 
    xOriginalOffsetFromMean + pca.MapToFeatureSpace(xCoeffArray);

  // Set the parameters of the parent surface
  xSurface->SetCoefficientArray(X.data_block());
}

IHyperSurface2D *PCACoefficientMask::GetComponentSurface(size_t iCoefficient)
{
  // Create a coefficient vector
  vnl_vector<double> xVariationCoeff(nModes, 0.0);
  xVariationCoeff[iCoefficient] = 1.0;
  return GetVariationSurface(xVariationCoeff.data_block());
}

void PCACoefficientMask::ReleaseComponentSurface(IHyperSurface2D *xSurface)
{
  ReleaseVariationSurface(xSurface);
}

IHyperSurface2D *PCACoefficientMask::GetVariationSurface(const double *xVariation)
{
  // Map the variation into the feature space
  vnl_vector<double> xCSVariation(xVariation, nModes);
  vnl_vector<double> xFSVariation 
    = pca.MapToFeatureSpaceZeroMean(xCSVariation);
  
  // Generate a variation surface using parent's method
  return xSurface->GetVariationSurface(xFSVariation.data_block());
}

void PCACoefficientMask::ReleaseVariationSurface(IHyperSurface2D *xVarSurface)
{
  xSurface->ReleaseVariationSurface(xVarSurface);
}

/******************************************************************************
 * AffineAndPCACoefficientMask
 *****************************************************************************/
AffineAndPCACoefficientMask::AffineAndPCACoefficientMask(
  const vnl_matrix<double> &mPCA, IBasisRepresentation2D *xSurface, size_t nModes)
: AffineTransformCoefficientMask(xSurface), 
  maskPCA(mPCA, xSurface, nModes),
  xPreAffineCoeffArray(xSurface->GetCoefficientArray(),
    xSurface->GetNumberOfCoefficients())
{
  nAffine = AffineTransformCoefficientMask::GetNumberOfCoefficients();
  nPCA = nModes;
  n = nAffine + nPCA;
  
  xCoeffArray.set_size(n);
  for(size_t i=0; i < nAffine; i++)
    xCoeffArray[i] = AffineTransformCoefficientMask::GetCoefficient(i);
  for(size_t j=0; j < nPCA; j++)
    xCoeffArray[j + nAffine] = maskPCA.GetCoefficient(j);
}

void
AffineAndPCACoefficientMask
::SetCoefficientArray(const double *xData)
{
  // Fill the coefficient array
  xCoeffArray.copy_in(xData);

  // Apply the PCA mask first
  maskPCA.SetCoefficientArray(xData  + nAffine);
  xPreAffineCoeffArray.copy_in(xSurface->GetCoefficientArray());

  // Apply the Affine mask next
  AffineTransformCoefficientMask::SetCoefficientArray(xData);
}

IHyperSurface2D *
AffineAndPCACoefficientMask
::GetComponentSurface(size_t iComp)
{
  // Same as the affine transform's variation surface
  if(iComp < nAffine)
    {
    // Use the current parameters as the PCA component of variation
    IHyperSurface2D *xVariationPCA 
      = xSurface->GetVariationSurface(xPreAffineCoeffArray.data_block());

    // Create a simple variation in the affine direction
    VectorType xAffineVarArray(nAffine, 0.0);
    xAffineVarArray[iComp] = 1.0;

    // Generate a combined variation
    return CreateAffineVariation(xVariationPCA, xAffineVarArray.data_block());
    }
  else
    {
    // Get a PCA-based component surface
    IHyperSurface2D *xVariationPCA = maskPCA.GetComponentSurface(iComp - nAffine);

    // This surface will have to be affine-transformed, but with no offsets
    MatrixType A(xCoeffArray.data_block(), nDim, nDim);

    // Now, create an affine wrapper on top of this surface
    return new AffineComponentSurface(xVariationPCA, nDim, A, 
      VectorType(nDim, 0.0), VectorType(nDim, 0.0));
    }
}

void
AffineAndPCACoefficientMask
::ReleaseComponentSurface(IHyperSurface2D *xVariation)
{
  // Delete the internal component
  AffineComponentSurface * xAffineSurf 
    = dynamic_cast<AffineComponentSurface *>(xVariation);

  // TODO: technically, this is wrong, I should be calling ReleaseVariationSurface
  // depending on which parameter was used to generate this
  delete xAffineSurf->xSource;

  // Delete the surface itself
  delete xAffineSurf;
}

IHyperSurface2D *AffineAndPCACoefficientMask
::GetVariationSurface(const double *xData)
{
  // Create the component dA * S
  IHyperSurface2D *xComponentPCA 
    = xSurface->GetVariationSurface(xPreAffineCoeffArray.data_block());

  IHyperSurface2D *surf1 =  CreateAffineVariation(xComponentPCA, xData);

  // Create the component A * dS
  IHyperSurface2D *xVariationPCA = maskPCA.GetVariationSurface(xData + nAffine);

  MatrixType A(xCoeffArray.data_block(), nDim, nDim);

  IHyperSurface2D *surf2 = new AffineComponentSurface(xVariationPCA, nDim, 
    A, VectorType(nDim, 0.0), VectorType(nDim, 0.0));

  // Combine the two components
  return new AdditiveComponentSurface(surf1, surf2);
}

void AffineAndPCACoefficientMask
::ReleaseVariationSurface(IHyperSurface2D *xVarSurface)
{
  // Cast pointers to the right classes
  AdditiveComponentSurface *xSumSurface = 
    dynamic_cast<AdditiveComponentSurface *> (xVarSurface);

  AffineComponentSurface *surf1 = 
    dynamic_cast<AffineComponentSurface *> (xSumSurface->surf1);

  AffineComponentSurface *surf2 = 
    dynamic_cast<AffineComponentSurface *> (xSumSurface->surf2);

  // Deallocate the source surfaces
  xSurface->ReleaseVariationSurface(surf1->xSource);
  maskPCA.ReleaseVariationSurface(surf2->xSource);

  // Delete the wrapper surfaces
  delete surf1;
  delete surf2;
  delete xSumSurface;
}

/******************************************************************************
 * Affine3DAndPCACoefficientMask
 *****************************************************************************/

Affine3DAndPCACoefficientMask
::Affine3DAndPCACoefficientMask(const vnl_matrix<double> &mPCA, 
  IBasisRepresentation2D *xSurface, size_t nModes)
: SelectionMedialCoefficientMask(
   xWrappedMask = new AffineAndPCACoefficientMask(mPCA, xSurface, nModes), 
   nModes+12, FillIndexArray(nModes) )
{
}

size_t *Affine3DAndPCACoefficientMask
::FillIndexArray(size_t nModes)
{
  xIndexArray = new size_t[12 + nModes];

  for(size_t i = 0; i < 12; i++)
    xIndexArray[i] = AffineTransform3DCoefficientMask::xIndexArray[i];
  for(size_t j = 0; j < nModes; j++)
    xIndexArray[12 + j] = 20 + j;

  return xIndexArray;
}
