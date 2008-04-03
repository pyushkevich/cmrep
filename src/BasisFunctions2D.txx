#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_qr.h>

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::GenericBasisRepresentation2D(size_t ncu, size_t ncv)
: C(NComponents, ncu, ncv)
{
  Initialize(ncu, ncv);
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::Initialize(size_t ncu, size_t ncv)
{
  // initialize the coefficient array
  this->ncu = ncu; this->ncv = ncv;
  nc = ncu * ncv;
  ncRaw= nc * NComponents;

  // initialize the evaluation arrays
  uEval.resize(ncu); vEval.resize(ncv);
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::SetNumberOfCoefficients(size_t ncu, size_t ncv)
{
  // Make a copy of the coefficients
  Index3D C1(NComponents, ncu, ncv);
  for(size_t i=0; i < ncu; i++) for(size_t j=0; j < ncv; j++)
    if(i < this->ncu && j < this->ncv)
      for(size_t k=0; k < NComponents; k++)
        C1(k,i,j) = C(k,i,j);

  C.resize(NComponents, ncu, ncv);
  C = C1;

  // initialize the coefficient array
  this->ncu = ncu; this->ncv = ncv;
  nc = ncu * ncv;
  ncRaw= nc * NComponents;

  // initialize the evaluation arrays
  uEval.resize(ncu); vEval.resize(ncv);
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::EvaluateDerivative(double u, double v, size_t ou, size_t ov, size_t c0, size_t nc, double *x)
{
  // No grid is available - must evaluate each basis function
  size_t iu, iv, k;
  for(iu = 0; iu < ncu; iu++) uEval[iu] = fu.Evaluate(u, iu, ou);
  for(iv = 0; iv < ncv; iv++) vEval[iv] = fv.Evaluate(v, iv, ov);

  // Clear the output array
  for(k = 0; k < nc; k++) x[k] = 0.0;

  // Compute the cross product of the bases
  unsigned int iOffset = 0;
  for(iv = 0; iv < ncv; iv++) for(iu = 0; iu < ncu; iu++) 
    {
    // Compute the product of the basis functions
    double Fuv = uEval[iu] * vEval[iv];

    // Get the coefficient corresponding to iu and iv    
    for(k = 0; k < nc; k++)
      x[k] += C[iOffset + k + c0] * Fuv;
    
    iOffset += NComponents;
    }
}


template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::PrintReport()
{
  size_t iu, iv, c0;

  for(c0 = 0; c0 < NComponents; c0++)
    {
    unsigned int iOffset = 0;
    cout << "F_" << c0 << "[u,v] = "; 
    for(iv = 0; iv < ncv; iv++) for(iu = 0; iu < ncu; iu++) 
      {
      if(iOffset > 0) cout << " + ";
      cout << C[iOffset + c0] << " * f[u," << iu << "] * f[v," << iv << "]";
      iOffset += NComponents;
      }
    cout << endl;
    }
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::SetEvaluationGrid(const VectorType &uu, const VectorType &vv)
{
  // Get the size of the grid
  size_t nu = uu.size(), nv = vv.size();

  // Initialize the evaluation grid to the specified dimensions
  uGrid.resize(nu * ncu * NOrder);
  vGrid.resize(nv * ncv * NOrder);
  uGridValues.resize(nu); 
  vGridValues.resize(nv);

  // Precompute the U-grid
  size_t iGrid = 0;
  for(size_t iu = 0; iu < nu; iu++)
    {
    uGridValues[iu] = uu[iu];
    for(size_t icu = 0; icu < ncu; icu++)
      for(size_t p = 0; p < NOrder; p++)
        uGrid[iGrid++] = fu.Evaluate(uGridValues[iu], icu, p);
    }
  
  // Precompute the V-grid
  iGrid = 0;
  for(size_t iv = 0; iv < nv; iv++)
    {
    vGridValues[iv] = vv[iv];
    for(size_t icv = 0; icv < ncv; icv++)
      for(size_t p = 0; p < NOrder; p++)
        vGrid[iGrid++] = fv.Evaluate(vGridValues[iv], icv, p);
    }
} 

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::GetEvaluationGrid(VectorType &uu, VectorType &vv) const
{
  uu.set_size(uGridValues.size());
  vv.set_size(vGridValues.size());
  for(size_t iu = 0; iu < uu.size(); iu++) 
    uu[iu] = uGridValues[iu];
  for(size_t iv = 0; iv < vv.size(); iv++)
    vv[iv] = vGridValues[iv];
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::EvaluateAtGridIndex(size_t iu, size_t iv, size_t ou, size_t ov, size_t c0, size_t nc, double *x)
{
  // Get a slice into the U grid. Indexing over this slice gets consecutive
  // basis functions of order ou at site iu
  // slice sGridU = 
  valarray<double> sGridU = uGrid[slice(NOrder * ncu * iu + ou, ncu, NOrder)]; 
  valarray<double> sGridV = vGrid[slice(NOrder * ncv * iv + ov, ncv, NOrder)];

  // Clear the output array
  size_t icu, icv, k;
  for(k = 0; k < nc; k++) x[k] = 0.0;

  // Compute the cross product of the bases
  size_t iOffset = 0;
  for(icv = 0; icv < ncv; icv++) for(icu = 0; icu < ncu; icu++) 
    {
    // Compute the product of the basis functions
    double Fuv = sGridU[icu] * sGridV[icv];

    // Get the coefficient corresponding to iu and iv
    for(k = 0; k < nc; k++)
      x[k] += C[iOffset + k + c0] * Fuv;

    iOffset += NComponents;
    }
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::FitToData(size_t n, size_t iComponent, double *uu, double *vv, double *xx)
{
  // The number of coefficients may not exceed the number of data points being
  // fitted (otherwise the fit will be nonsensical)
  size_t nu = ncu; size_t nv = ncv;
  while(nu * nv > n)
    { nu--; nv--; }

  // Get the number of unknown coefficients
  size_t iu, iv, i = 0, j, k;
  size_t nUnkowns = nu * nv;

  cout << "Fitting data with " << nUnkowns << " unknowns to " << n << " points." << endl;

  // Compute the basis for each point
	vnl_matrix<double> Z(nUnkowns, n);
  for(iu = 0, i = 0; iu < nu; iu++) for(iv = 0; iv < nv; iv++, i++)
    for(k = 0; k < n; k++) 
      Z[i][k] = fu.Evaluate(uu[k], iu, 0) * fv.Evaluate(vv[k], iv, 0); 


  // Allocate the matrix A and vector b
	vnl_matrix<double> A(nUnkowns, nUnkowns);
	vnl_vector<double> b(nUnkowns);
  
  // Set the elements of A and b
  unsigned int offset = 0;
  for(j = 0; j < nUnkowns; j++)
    {
    // Compute the b vector
    b[j] = 0;
    for(k = 0; k < n; k++) 
      b[j] += xx[k] * Z[j][k]; 

    // Compute the elements of the A matrix
    for(i = 0; i < nUnkowns; i++)
      {
      A[j][i] = 0.0;
      for(k = 0; k < n; k++) 
        A[j][i] += Z[i][k] * Z[j][k];
      ++offset;
      }
    }

  // Solve the system Ax = b (LU decomposition)
	vnl_qr<double> qr(A);
	vnl_vector<double> y = qr.solve(b);

  // Solution has been placed into B. Map it to the coefficients
  i = 0;
  for(iu = 0; iu < nu; iu++) for(iv = 0; iv < nv; iv++)
    C(iComponent, iu, iv) = y[i++];
}


template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::SaveToRegistry(Registry &R)
{
  // Store the main information
  R["Size.U"] << ncu;
  R["Size.V"] << ncv;
  R["Components"] << NComponents;

  // Store each coefficient
  for(size_t i = 0; i < ncu; i++) for(size_t j = 0; j < ncv; j++)
      for(size_t k = 0; k < NComponents; k++)
        R [ R.Key("Coefficient[%d][%d][%d]", i, j, k) ] << C(k, i, j);
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
bool
GenericBasisRepresentation2D<NComponents, NOrder, BasisFunctionU, BasisFunctionV>
::ReadFromRegistry(Registry &R)
{
  size_t mu = R["Size.U"][0];
  size_t mv = R["Size.V"][0];
  size_t mk = R["Components"][0];

  // Check the parameters
  if(mu == 0 || mv == 0 || mk != NComponents) return false;

  // Expand the number of coefficients
  if(mu < ncu) mu = ncu; if(mv < ncv) mv = ncv;

  // Reinitialize the coefficient array
  C.resize(NComponents, mu, mv);
  Initialize(mu, mv);
  
  // Load the coefficients
  for(size_t i = 0; i < mu; i++) for(size_t j = 0; j < mv; j++)
    for(size_t k = 0; k < NComponents; k++)
      C(k, i, j) = R[ R.Key("Coefficient[%d][%d][%d]", i, j, k) ][0.0];
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
size_t
GenericBasisRepresentation2D<NComponents,NOrder,BasisFunctionU,BasisFunctionV>
::SingleFunctionAdapter
::GetNumberOfDimensions()
{
  return xParent->GetNumberOfDimensions();
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents,NOrder,BasisFunctionU,BasisFunctionV>
::SingleFunctionAdapter
::Evaluate(double u, double v, double *x)
{
  EvaluateDerivative( u, v, 0, 0, 0, NComponents, x); 
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents,NOrder,BasisFunctionU,BasisFunctionV>
::SingleFunctionAdapter
::EvaluateDerivative( double u, double v, size_t ou, size_t ov, 
  size_t c0, size_t nc, double *x)
{
  // Check that the component is in range
  for(size_t k = 0; k < nc; k++) x[k] = 0.0;
  if(c0 <= iComp && iComp < c0 + nc)
    {
    x[iComp - c0] = xParent->fu.Evaluate(u, icu, ou) * 
      xParent->fv.Evaluate(v, icv, ov);
    }
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents,NOrder,BasisFunctionU,BasisFunctionV>
::SingleFunctionAdapter
::EvaluateAtGridIndex(size_t iu, size_t iv, size_t ou, size_t ov, 
  size_t c0, size_t nc, double *x)
{
  // Check that the component is in range
  for(size_t k = 0; k < nc; k++) x[k] = 0.0;
  if(c0 <= iComp && iComp < c0 + nc)
    {
    x[iComp - c0] = 
      xParent->uGrid[ NOrder * ( xParent->ncu * iu + icu ) + ou ] *
      xParent->vGrid[ NOrder * ( xParent->ncv * iv + icv ) + ov ];
    }
}

template< size_t NComponents, size_t NOrder, typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents,NOrder,BasisFunctionU,BasisFunctionV>
::SingleFunctionAdapter
::GetEvaluationGrid(VectorType &uu, VectorType &vv) const 
{
  xParent->GetEvaluationGrid(uu, vv);
}

template<
  size_t NComponents, size_t NOrder, 
  typename BasisFunctionU, typename BasisFunctionV >
IHyperSurface2D *
GenericBasisRepresentation2D<NComponents,NOrder,BasisFunctionU,BasisFunctionV>
::GetComponentSurface(size_t icu, size_t icv, size_t iComp)
{
  SingleFunctionAdapter *xAdapt = new SingleFunctionAdapter;
  xAdapt->xParent = this;
  xAdapt->icu = icu;
  xAdapt->icv = icv;
  xAdapt->iComp = iComp;

  return xAdapt;
}

template<
  size_t NComponents, size_t NOrder, 
  typename BasisFunctionU, typename BasisFunctionV >
IHyperSurface2D*
GenericBasisRepresentation2D<NComponents,NOrder,BasisFunctionU,BasisFunctionV>
::GetComponentSurface(size_t iComponent)
{
  size_t icu, icv, iComp;
  GetRawCoefficientIndices(iComponent, icu, icv, iComp);
  return GetComponentSurface(icu, icv, iComp);
}

template< 
  size_t NComponents, size_t NOrder, 
  typename BasisFunctionU, typename BasisFunctionV >
void
GenericBasisRepresentation2D<NComponents,NOrder,BasisFunctionU,BasisFunctionV>
::ReleaseComponentSurface(IHyperSurface2D *xSurface)
{
  delete xSurface;
}
