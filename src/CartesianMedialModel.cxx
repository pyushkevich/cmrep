#include "CartesianMedialModel.h"
#include "MedialException.h"
#include "MedialPDEMasks.h"
#include <cmath>
#include <algorithm>
#include "smlmath.h"


using namespace std;

void SparseMultiply(
  unsigned int n, unsigned int *rowIndex, unsigned int *colIndex, 
  double *values, double *x, double *y)
{
	for(unsigned int r = 1; r <=n ; r++)
    	{
        y[r-1] = 0;
    	for(unsigned int k = rowIndex[r-1]; k < rowIndex[r]; k++) 
      		{
      		// r and c are 1-based
      		int c = colIndex[k - 1] ;	
      		double v = values[k - 1];
      		y[r-1] += v * x[c - 1];
      		}
    	}	
}

double SparseLinearTest(unsigned int n, unsigned int *rowIndex, 
  unsigned int *colIndex, double *values, double *x, double *y,double *b)
{
  SparseMultiply(n, rowIndex, colIndex, values, x, y);
  double maxdiff = 0.0;
  for(unsigned int i = 0; i < n; i++)
    if(fabs(b[i] - y[i]) > maxdiff)
      maxdiff = fabs(b[i] - y[i]);
  return maxdiff;
}

void DumpSparseMatrix(unsigned int n, unsigned int *rowIndex, 
  unsigned int *colIndex, double *values)
{
  char ch = '{';
  cout << "SparseArray[";
  for(unsigned int r = 1; r <=n ; r++)
    {
    for(unsigned int c = rowIndex[r-1]; c < rowIndex[r]; c++) 
      {
      double v = values[c - 1];
      if(v < 0.000001 && v > -0.000001) v = 0;
      cout << ch << " {" << r << "," << colIndex[c-1] << "}->" << v;
      ch = ',';
      }
      cout << endl;
    }
  cout << "}]" << endl;
}

CartesianMedialModel
::CartesianMedialModel(size_t nu, size_t nv, double xScale, size_t pu, size_t pv)
{
  // Compute the total number of grid points
  size_t tu = nu + 2 * pu;
  size_t tv = nv + 2 * pv;

  // Allocate the initialization grids
  Vec uuGrid(tu), vvGrid(tv);

  // Set the endpoints
  uuGrid[0] = vvGrid[0] = 0.0;
  uuGrid[tu-1] = vvGrid[tv-1] = 1.0;

  // Set the regularly spaced points
  double du = 1.0 / (nu - 1);
  double dv = 1.0 / (nv - 1);
  size_t iu, iv;

  double zu = du;
  for(iu = 1; iu <= nu-2; iu++, zu+=du)
    uuGrid[iu + pu] = zu;

  double zv = dv;
  for(iv = 1; iv <= nv-2; iv++, zv+=dv)
    vvGrid[iv + pv] = zv;

  // Set the irregularly spaced points
  for(iu = 0; iu < pu; iu++)
    {
    double duScaled = du * pow(xScale, iu+1.0); 
    uuGrid[pu - iu] = duScaled;
    uuGrid[(tu - 1) - (pu - iu)] = 1.0 - duScaled; 
    }
  
  for(iv = 0; iv < pv; iv++)
    {
    double dvScaled = dv * pow(xScale, iv+1.0); 
    vvGrid[pv - iv] = dvScaled;
    vvGrid[(tv - 1) - (pv - iv)] = 1.0 - dvScaled; 
    }

  // Apply random jitter
  /*
  for(iu = 1; iu <= nu-2; iu++)
    uuGrid[iu + pu] += 0.45 * (-1.0 + rand() * 2.0 / RAND_MAX) * du;

  for(iv = 1; iv <= nv-2; iv++)
    vvGrid[iv + pv] += 0.45 * (-1.0 + rand() * 2.0 / RAND_MAX) * dv;
  */

  Initialize(uuGrid, vvGrid);
}

CartesianMedialModel
::CartesianMedialModel(const Vec &uGrid, const Vec &vGrid)
{
  Initialize(uGrid, vGrid);
}

void 
CartesianMedialModel
::Initialize(const Vec &uGrid, const Vec &vGrid)
{
  // Copy the grid vectors
  this->uGrid = uGrid; this->vGrid = vGrid;

  // Record the size of the grid
  m = uGrid.size(); n = vGrid.size();

  // Set the iteration context
  if(this->xIterationContext != NULL) delete(this->xIterationContext);
  this->xIterationContext = new CartesianGridMedialIterationContext(m, n);

  // Number of sites
  nSites = m * n;

  // Initialize the mask and site arrays
  xSites.resize(nSites, NULL);
  xMasks.resize(nSites, NULL);

  // Set the size of the size index
  xSiteIndex.set_size(m, n);

  // Set up the initial guess
  xDefaultInitSoln.set_size(m, n);

  // Initialize all the masks and sites
  size_t iSite = 0, i, j;
  for(i = 0; i < m; i++) for(j = 0; j < n; j++) 
    {
    // Flags indicating adjacence to a border
    bool uBorder = (i == 0 || i == m-1), vBorder = (j == 0 || j == n-1);

    // cout << i << ", " << j << " : ";

    // Check if we are at the border
    if(uBorder || vBorder)
      {
      if(uBorder && vBorder)
        {
        // Create a corner mask
        // cout << "Corner Mask" << endl;
        xMasks[iSite] = new CornerFDMask();
        }
      else
        {
        // Create a border mask
        // cout << "Border Mask" << endl;
        xMasks[iSite] = new BorderFDMask();
        if(vBorder)
          xMasks[iSite]->TransposeMask();
        }

      // In either case, flip the mask over
      xMasks[iSite]->FlipMask(i != 0, j != 0);
      }

    // We are at an internal site
    else
      {
      // Check the difference in grid spacing
      double ddu = fabs((uGrid[i+1] - uGrid[i]) - (uGrid[i] - uGrid[i-1]));
      double ddv = fabs((vGrid[j+1] - vGrid[j]) - (vGrid[j] - vGrid[j-1]));
      // Check if the grid is uniform
      if(ddu < 1.0e-13 && ddv < 1.0e-13)
        {
        // Create a uniform mask
        // cout << "Uniform Mask" << endl;
        xMasks[iSite] = new UniformFDMask();
        }
      else
        {
        // Create a non-uniform mask, and orient it correctly
        // cout << "Nonuniform Mask" << endl;
        xMasks[iSite] = new DoublyNonuniformFDMask();
        xMasks[iSite]->FlipMask( i > m/2, j > n/2 );
        }
      }

    // Set the location of the mask and specify the grid
    xMasks[iSite]->SetLocation(i, j);
    xMasks[iSite]->SortNodes();
    xMasks[iSite]->ComputeWeights(uGrid.data_block(), vGrid.data_block());
    xMasks[iSite]->OptimizeWeights(uGrid.size(), vGrid.size());

    // Compute the initial solution value as the distance from the nearest
    // edge
    // double uMin = i > m/2 ? uGrid[m-(i+1)] - uGrid[m-(i+2)] : uGrid[i+1] - uGrid[i];
    // double vMin = j > n/2 ? vGrid[n-(j+1)] - vGrid[n-(j+2)] : vGrid[j+1] - vGrid[j];
    
    // xDefaultInitSoln[i][j] = sqrt( uMin * uMin + vMin * vMin );
    xDefaultInitSoln[i][j] = 1.0;

    // xMasks[iSite]->PrintReport();

    // Make sure that every mask is in bounds
    for(size_t k = 0; k < xMasks[iSite]->Size(); k++)
      {
      int ku = xMasks[iSite]->GetOffsetU(k) + xMasks[iSite]->GetLocationU();
      int kv = xMasks[iSite]->GetOffsetV(k) + xMasks[iSite]->GetLocationV();
      if(ku < 0 || ku >= (int) m || kv < 0 || kv >= (int) n) 
        cout << "Site " << i << ", " << j << " is out of bounds!" << endl;
      }

    // Create a site wrapping around the mask
    if(uBorder || vBorder)
      xSites[iSite] = new FDBorderSite(xMasks[iSite]);
    else
      xSites[iSite] = new FDInternalSite(xMasks[iSite]);

    // Associate the raw offset with the u-v index
    xSiteIndex[i][j] = iSite++;
    }

  // Initialize the sparse matrix row index (use Fortran-based indexing)
  xRowIndex = new int[nSites + 1];
  xRowIndex[0] = 1;
  for(iSite = 0; iSite < nSites; iSite++)
    xRowIndex[iSite+1] = xRowIndex[iSite] + xMasks[iSite]->Size();
  
  // Initialize the column index
  nSparseEntries = xRowIndex[nSites] - 1;
  xColIndex = new int[nSparseEntries];
  
  // Initialize the sparse data
  xSparseValues = new double[nSparseEntries];
  memset(xSparseValues, 0, sizeof(double) * nSparseEntries);

  // Compute the column index and the sparse data
  for(iSite = 0; iSite < nSites; iSite++)
    {
    // Number of neighbors in this node
    size_t nEntries = xRowIndex[iSite+1] - xRowIndex[iSite];

    for(unsigned int iEntry = 0; iEntry < nEntries; iEntry++)
      {
      // Get the coordinates of the given neighbor node
      int iu = xMasks[iSite]->GetLocationU() + xMasks[iSite]->GetOffsetU(iEntry);
      int iv = xMasks[iSite]->GetLocationV() + xMasks[iSite]->GetOffsetV(iEntry);
      int iNeighbor = xSiteIndex[iu][iv];

      // Set the column entry in the sparse matrix
      xColIndex[ xRowIndex[iSite] + iEntry - 1 ] = iNeighbor + 1;
      }
    }

  // Initialize our three vectors
  b.set_size(m, n);
  y.set_size(m, n);
  eps.set_size(m, n);
  dy.set_size(m, n);
  zTest.set_size(m, n);

  // Initialize the medial atom array and atom grid
  xAtoms = new MedialAtom[nSites];

  // Compute the initial solution
  SetDefaultInitialGuess(1);

  // Set the surface to NULL
  xSurface = NULL;

  // Create a solver 
  solver = SparseSolver::MakeSolver(false);
}

CartesianMedialModel
::~CartesianMedialModel()
{
  // Delete the iteration context
  delete this->xIterationContext;

  // Delete the atoms
  delete this->xAtoms;

  // Delete all the sites and masks
  for(size_t i = 0; i < nSites; i++)
    {
    delete xMasks[i];
    delete xSites[i];
    }

  // Delete other indices
  delete xSparseValues;
  delete xRowIndex;
  delete xColIndex;

  // Delete the surface
  if(xSurface != NULL) delete xSurface;

  // Delete solver
  delete solver;
}

void
CartesianMedialModel
::SetDefaultInitialGuess(double xMagnitude)
{
  // Set the initial guess to default values provided by the sites
  xInitSoln = xMagnitude * xDefaultInitSoln;
}

void 
CartesianMedialModel
::SetSolutionAsInitialGuess()
{
  // Make the current solution the initial guess. This makes Newton converge
  // quickly when exploring near a point
  xInitSoln = y;
  flagReuseLastSolution = true;
}

void 
CartesianMedialModel
::AdoptMedialSurface(IBasisRepresentation2D *xSurface)
{
  // Delete the old surface
  if(this->xSurface != NULL && this->xSurface != xSurface) 
    delete this->xSurface;

  // Adopt the new surface
  this->xSurface = xSurface;

  // Tell the surface where it will be evaluated
  xSurface->SetEvaluationGrid(uGrid, vGrid);
}

void CartesianMedialModel::InitializeSiteGeometry()
{
  size_t i, j;
  
  // Initialize each site with the current surface properties
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    unsigned int iGrid = GetGridAtomIndex(i, j);
    unsigned int iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &xAtom = xAtoms[iGrid];

    // Set the atoms' domain coordinates
    xAtom.u = uGrid[i]; xAtom.v = vGrid[j];
    xAtom.uIndex = i; xAtom.vIndex = j;
    
    // Compute the surface jet and the laplacian    
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, xAtom.X.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, xAtom.Xu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, xAtom.Xv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, xAtom.Xuu.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, xAtom.Xuv.data_block());
    xSurface->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, xAtom.Xvv.data_block());

    // Compute the differential geometric tensors
    xAtom.ComputeDifferentialGeometry();

    // Compute the normal vector
    xAtom.ComputeNormalVector();

    // Compute the laplacian of R 
    xSurface->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &xAtom.xLapR);
    
    // Compute the solution at this point
    xSites[iSite]->SetGeometry( &xAtom.G, xAtom.xLapR);
    }
}

void
CartesianMedialModel
::ReconstructAtoms(const Mat &ySolution)
{
  // Keep track of bad atoms
  size_t nBadsF = 0, nBadsGradR = 0;

  // Keep track of invalid atoms
  bool flagAllAtomsValid = true;
  
  // Once the iterations are complete, reconstruct the atoms
  for(unsigned int i = 0; i < m; i++) for(unsigned int j = 0; j < n; j++)
    {
    // Map to a grid index and to our own index
    unsigned int iGrid = GetGridAtomIndex(i, j);
    unsigned int iLocal = xSiteIndex[i][j];
    
    // The medial atom to update
    MedialAtom &xAtom = xAtoms[iGrid];

    // The case where phi is negative is undefined
    if( ySolution[i][j] >= 0 )
      {
      // Compute the derivatives of R using finite differences
      xSites[iLocal]->ComputeEquation(ySolution);
      xAtom.F = 
        xMasks[iLocal]->ComputeOneJet(ySolution.data_block(), xAtom.Fu, xAtom.Fv);

      // Compute the boundary properties of the medial point
      xAtom.ComputeBoundaryAtoms(xSites[iLocal]->IsBorderSite());

      // Report bad atoms
      if(!xAtom.flagValid) nBadsGradR++;

      /*
      if(!xAtom.flagValid)
        cout << "Invalid atom at " << xAtom.u << ", " << xAtom.v <<
          "(1 - |gradR|^2 = " << 1.0 - xAtom.xGradRMagSqr << ")" << endl;
          */
      }
    else
      {
      // What are we supposed to do?
      nBadsF++;

      /* 
      cout << "Negative F at " << xAtom.u << ", " << xAtom.v << 
        " ( F = " << ySolution[i][j] << ")" << endl;
        */
      xAtom.R = xAtom.F = xAtom.Fu = xAtom.Fv = 0.0;
      xAtom.flagValid = false;
      
      }


    // Update the valid atoms flag
    flagAllAtomsValid &= xAtom.flagValid;
    }

  // Return value: whether all atoms are valid
  if(!flagAllAtomsValid)
    {
    std::ostringstream oss;
    oss << "Invalid atoms (" << nBadsF << " with negative F, " 
      << nBadsGradR << " with |GradR| > 1)";
    throw MedialModelException(oss.str().c_str());
    }
}

void ArrayMinMax(double *array, size_t n, double &xMin, double &xMax)
{
	xMax = 0, xMin = 1e100;
  for(size_t q = 0; q < n; q++)
  	{
    if(xMax < fabs(array[q])) xMax = fabs(array[q]);
    if(xMin > fabs(array[q])) xMin = fabs(array[q]);
    }
}

/** This computes the right hand side of the site equations, given phi=x */
double CartesianMedialModel::ComputeNewtonRHS(const Mat& x, Mat &b)
{
  size_t i, j;
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    b[i][j] = - xSites[xSiteIndex[i][j]]->ComputeEquation(x);
    if(std::isnan(b[i][j]))
      throw MedialModelException("NAN in CartesianMedialModel::ComputeNewtonRHS");
    }
  return dot_product(b, b);
}

CartesianMedialModel::Vec 
CartesianMedialModel::GetHintArray() const
{
  Vec xHint(nSites, 0.0);
  for(size_t i = 0; i < m; i++) for(size_t j = 0; j < n; j++)
    xHint[GetGridAtomIndex(i,j)] = y[i][j];
  return xHint;
}

void CartesianMedialModel::SolveOnce(const double *xHint, double delta)
{
  size_t iIter;
  double epsMax, bMax, bMagSqr;
  
  // Initialize epsilon to zero
  eps.fill(0.0);
  
  // Copy the initial solution to the current solution
  if(xHint != NULL)
    for(unsigned int i = 0; i < m; i++) for(unsigned int j = 0; j < n; j++)
      y[i][j] = xHint[GetGridAtomIndex(i,j)];

  // Compute the right hand side, put it into b
  bMagSqr = ComputeNewtonRHS(y, b);

  // We are now ready to perform the Newton loop
  for(iIter = 0; iIter < 50; iIter++)
    {
    // Compute the Jacobian matrix
    for(size_t iSite = 0; iSite < nSites; iSite++)
      xSites[iSite]->
        ComputeDerivative(y, xSparseValues + xRowIndex[iSite] - 1, iSite+1);
      
    // Perform the symbolic factorization only for the first iteration
    tSolver.Start();
    if(iIter == 0)
      solver->SymbolicFactorization(nSites, xRowIndex, xColIndex, xSparseValues);

    // Compute the Jacobian inverse
    solver->NumericFactorization(xSparseValues);
    solver->Solve(b.data_block(), eps.data_block());
    tSolver.Stop();

    // A plus means solver step
    // cout << "+";

    // Advance y to the new Newton position
    y += eps; 
       
    // Compute the solution at the full Newton step
    double bMagSqrTest = ComputeNewtonRHS(y, b);

    // Perform backtracking if necessary (this is boneheaded backtracking but
    // it appears to work). TODO: Implement lnsrch from NRC
    double lambda = 0.5;
    while(bMagSqrTest > bMagSqr && lambda > 1e-4)
      {
      // Go back along eps
      y -= lambda * eps; lambda *= 0.5;

      // Compute the right hand side again
      bMagSqrTest = ComputeNewtonRHS(y, b);

      // A "-" means backtrack step
      // cout << "-";
      }

    // Store the new b magnitude value
    bMagSqr = bMagSqrTest;
    // cout << " " << bMagSqr << " ";

    // Get the largest error (eps)
    epsMax = eps.array_inf_norm();
    bMax = b.array_inf_norm();

    // Print the statistics
    // cout << "-----------" << endl;
    // cout << "Step " << iIter << ": " << endl;
    // cout << "  Largest Epsilon: " << epsMax << endl;
    // cout << "  Largest Eqn Error: " << bMax << endl; 

    // Convergence is defined when epsilon is smaller than some threshold
    if(bMax < delta || epsMax < delta) 
      break;
    }

  // Let the user know if we can't converge on the root
  if(bMax > delta && epsMax > delta)
    {
    cout << "  *** CONVERGENCE FAILURE *** ";
    cout << " epsMax = " << epsMax << "; bMax = " << bMax << endl;
    throw MedialModelException("Cartesian Medial PDE failed to converge");
    }
}

void
CartesianMedialModel
::Solve(const double *xHint, double delta)
{
  // Intialize the sites
  InitializeSiteGeometry();

  // Just solve once - no initialization tricks!
  SolveOnce(xHint, delta);
  
  // Reconstruct the medial atoms
  ReconstructAtoms(y);
}

void
CartesianMedialModel
::SetVariationalBasis(const Mat &xBasis)
{
  // Allocate the array of terms linearly dependent on the variation
  xVariationalBasis = VariationalBasisRep(xBasis.rows(), VariationRep(n*m));

  // Loop over the variations
  for(size_t var = 0; var < xBasis.rows(); var++)
    {
    // The current variation
    Vec xVariation = xBasis.get_row(var);

    // Get a variation surface corresponding to this variation
    IHyperSurface2D *fnVariation = xSurface->GetVariationSurface(xVariation.data_block());

    // Evaluate the variation over the atom grid
    for(size_t i = 0; i < m; i++) for(size_t j = 0; j < n; j++)
      {
      // Get the index of the site
      size_t iGrid = GetGridAtomIndex(i, j);

      // Access the medial atom underneath
      VariationalBasisAtomData &vbad = xVariationalBasis[var][iGrid];

      // Evaluate the variation and its derivatives
      fnVariation->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, vbad.X.data_block());
      fnVariation->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, vbad.Xu.data_block());
      fnVariation->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, vbad.Xv.data_block());
      fnVariation->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, vbad.Xuu.data_block());
      fnVariation->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, vbad.Xuv.data_block());
      fnVariation->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, vbad.Xvv.data_block());
      fnVariation->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &vbad.xLapR);
      }

    // Release the variation surface
    xSurface->ReleaseVariationSurface(fnVariation);
    }
}

void
CartesianMedialModel
::BeginGradientComputation()
{
  size_t i, j;

  // Reset the array of derivative terms
  xTempDerivativeTerms = vector<MedialAtom::DerivativeTerms>(nSites);

  // Prepare for the gradient computation
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = GetGridAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Get the medial atom
    MedialAtom &xAtom = xAtoms[iGrid];

    // Initialize the atom gradient terms
    xAtom.ComputeCommonDerivativeTerms(xTempDerivativeTerms[iSite]);

    // Compute the matrix for linear solver
    xSites[iSite]->ComputeVariationalDerivativeMatrix(
      y, xSparseValues + xRowIndex[iSite] - 1, &xAtom);
    }

  // Factorize the matrix so that the linear system associated with each
  // variational derivative can be solved instantly
  clock_t t0 = clock();
  tSolver.Start();
  solver->NumericFactorization(xSparseValues);
  tSolver.Stop();
  cout << " [FCT: " << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;
}

void
CartesianMedialModel
::ComputeAtomVariationalDerivative(size_t ivar, MedialAtom *dAtoms)
{
  size_t i, j; 

  // Compute the right hand side for the derivative computation
  Vec rhs(nSites, 0.0), soln(nSites, 0.0);
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = GetGridAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Get the reference to the atom's derivative
    MedialAtom &a  = xAtoms[iGrid];
    MedialAtom &da = dAtoms[iGrid];

    // Get the reference to the way this variation affects this atom
    VariationalBasisAtomData &vbad = xVariationalBasis[ivar][iGrid];

    // Copy the derivatives of X-jet into da
    da.X = vbad.X; da.Xu = vbad.Xu; da.Xv = vbad.Xv;
    da.Xuu = vbad.Xuu; da.Xuv = vbad.Xuv; da.Xvv = vbad.Xvv;
    da.xLapR = vbad.xLapR;

    // Set the atoms' domain coordinates
    da.u = uGrid[i]; da.v = vGrid[j];
    da.uIndex = i; da.vIndex = j;

    // Compute the derivative of the atom's metric tensor
    a.ComputeMetricTensorDerivatives(da);
    a.ComputeChristoffelDerivatives(da);

    // Compute the right hand side
    rhs[iSite] = xSites[iSite]->ComputeVariationalDerivativeRHS(y, &a, &da);
    }

  // Solve for the derivative of phi 
  solver->Solve(rhs.data_block(), soln.data_block());

  // Access the phi matrix
  vnl_matrix_ref<double> phi(m, n, soln.data_block());
    
  // Compute phi-dependent atom quantities
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = GetGridAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &da = dAtoms[iGrid];

    // Compute the gradient of phi for the new atom
    da.F = xMasks[iSite]->ComputeOneJet(phi.data_block(), da.Fu, da.Fv);

    // Compute the rest of the atom derivative
    xAtoms[iGrid].ComputeBoundaryAtomDerivatives(da, xTempDerivativeTerms[iSite]);
    }
}

/*
void
CartesianMedialModel
::PrepareAtomsForVariationalDerivative(const Vec &xVariation, MedialAtom *dAtoms) const
{
  // Get a variation surface corresponding to this variation
  IHyperSurface2D *fnVariation = xSurface->GetVariationSurface(xVariation.data_block());

  // Evaluate the variation over the atom grid
  for(size_t i = 0; i < m; i++) for(size_t j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = GetGridAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &dAtom = dAtoms[iGrid];

    // Set the atoms' domain coordinates
    dAtom.u = uGrid[i]; dAtom.v = vGrid[j];
    dAtom.uIndex = i; dAtom.vIndex = j;

    // Evaluate the variation and its derivatives
    fnVariation->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, dAtom.X.data_block());
    fnVariation->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, dAtom.Xu.data_block());
    fnVariation->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, dAtom.Xv.data_block());
    fnVariation->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, dAtom.Xuu.data_block());
    fnVariation->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, dAtom.Xuv.data_block());
    fnVariation->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, dAtom.Xvv.data_block());
    fnVariation->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &dAtom.xLapR);
    }

  // Release the variation surface
  xSurface->ReleaseVariationSurface(fnVariation);
}

void
CartesianMedialModel
::ComputeVariationalGradient(vector<MedialAtom *> &dAtomArray)
{
  size_t i, j, k, N = dAtomArray.size();

  // Create an array of right-hand sides
  Mat rhs(N, nSites), soln(N, nSites);

  // Create an array of helper structs for atom computation
  vector<MedialAtom::DerivativeTerms> agt(nSites);

  // Prepare for the gradient computation
  double t0 = clock();
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = GetGridAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Get the medial atom
    MedialAtom &xAtom = xAtoms[iGrid];

    // Initialize the atom gradient terms
    xAtom.ComputeCommonDerivativeTerms(agt[iSite]);

    // Compute the matrix for linear solver
    xSites[iSite]->ComputeVariationalDerivativeMatrix(
      y, xSparseValues + xRowIndex[iSite] - 1, &xAtom);

    // For each variation compute atom's geometry gradient and RHS
    for(k = 0; k < N; k++)
      {
      // Get the reference to the atom's derivative
      MedialAtom &dAtom = dAtomArray[k][iGrid];

      // Compute the derivative of the atom's metric tensor
      xAtom.ComputeMetricTensorDerivatives(dAtom);
      xAtom.ComputeChristoffelDerivatives(dAtom);
      
      // Compute the right hand side
      rhs[k][iSite] = 
        xSites[iSite]->ComputeVariationalDerivativeRHS(y, &xAtom, &dAtom);
      }
    }
  cout << " [RHS: " << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;

  // Solve the linear system for all the RHS
  t0 = clock();
  tSolver.Start();
  xPardiso.NumericFactorization(xSparseValues);
  xPardiso.Solve(N, rhs.data_block(), soln.data_block());
  tSolver.Stop();
  cout << " [SLV: " << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;

  // For each atom, compute the boundary derivatives
  t0 = clock();
  for(k = 0; k < N; k++)
    {
    // Access the phi matrix
    vnl_matrix_ref<double> phi(m, n, soln.data_array()[k]);
    
    // Compute phi-dependent atom quantities
    for(i = 0; i < m; i++) for(j = 0; j < n; j++)
      {
      // Get the index of the site
      size_t iGrid = GetGridAtomIndex(i, j);
      size_t iSite = xSiteIndex[i][j];

      // Access the medial atom underneath
      MedialAtom &dAtom = dAtomArray[k][iGrid];

      // Compute the gradient of phi for the new atom
      dAtom.F = xMasks[iSite]->ComputeOneJet(phi.data_block(), dAtom.Fu, dAtom.Fv);

      // Compute the rest of the atom derivative
      xAtoms[iGrid].ComputeBoundaryAtomDerivatives(dAtom, agt[iSite]);
      }
    }
  cout << " [BND: " << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;
}
*/

double fnTestF01(double u, double v)
{
  return exp(u + v);
}

void fnTestF01Jet(double u, double v, double jet[])
{
  jet[0] = exp(u + v);
  jet[1] = jet[2] = jet[3] = jet[4] = jet[5] = jet[0];
}

void CartesianMedialModel
::TestFiniteDifferenceConvergence()
{
  size_t i, j;
  
  // Create a field of values Phi
  Mat F(m, n);

  // Compute the values
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    F[i][j] = fnTestF01(uGrid[i], vGrid[j]);

  Vec xMaxDiff(6, 0.0);

  // Compute the partial derivatives using all the available masks
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // The jet vectors
    Vec j1(6), j2(6);
    
    // Compute the finite difference jet
    j1[0] = xMasks[xSiteIndex[i][j]]->ComputeTwoJet(F, j1[1], j1[2], j1[3], j1[4], j1[5]);
    
    // Compute the actual derivative
    fnTestF01Jet(uGrid[i], vGrid[j], j2.data_block());

    // Compute the differences in the jets
    Vec xDiff = (j1-j2).apply(fabs);

    int kn = xSites[xSiteIndex[i][j]]->IsBorderSite() ? 3 : 6;
    for(int k = 0; k < kn; k++)
      {
      if(xMaxDiff[k] < xDiff[k])
        {
        xMaxDiff[k] = xDiff[k];
        cout << "Max changed in " << k << " to " << xDiff[k] << " at " << i << "," << j << endl;
        }
      }
    }

  cout << "CartesianMedialModel FD Test: " << xMaxDiff << endl;
}

void
CartesianMedialModel
::ReadFromRegistry(Registry &R)
{
  // Ensure that the grid type matches this model
  if(R["Grid.Type"][""] != Registry::StringType("Cartesian"))
    throw ModelIOException(
      "CartesianMedialModel can not read Grid.Type that is not Cartesian");

  // Read the grid dimensions
  unsigned int m = R["Grid.Size.U"][64];
  unsigned int n = R["Grid.Size.V"][64];

  // Try reading the 

  // Ensure that the grid size matches the model (in the original code there was
  // no support for grid resizing, and we keep it that way for compatibility)
  if(this->m != m || this->n != n)
    throw ModelIOException(
      "Grid size does not match model when reading CartesianMedialModel");

  // Store the fourier coefficient information
  if(R["SurfaceModel"][""] == Registry::StringType("Fourier"))
    {
    // Read the surface from the parameters
    xSurface->ReadFromRegistry(R.Folder("Fourier"));

    // Set its evaluation grid
    xSurface->SetEvaluationGrid(uGrid, vGrid);
    }
  else throw ModelIOException(
    "SurfaceModel is not Fourier when reading CartesianMedialModel");

  // Read the phi matrix is it's available
  size_t nSites = this->GetNumberOfAtoms();
  if(R["Grid.PhiAvailable"][false] && nSites == R.Folder("Grid.Phi").GetArraySize())
    {
    // Read the phi array
    vnl_matrix<double> xPhi(m, n);
    R.Folder("Grid.Phi").GetArray(xPhi.data_block(), 0.0);

    // Pass the phi values to atoms (and to the phi-array)
    for(size_t i = 0; i < m; i++) for(size_t j = 0; j < n; j++)
      xAtoms[GetGridAtomIndex(i,j)].F = y[i][j] = xPhi[i][j];
    }
  else
    {
    // Initialize the atoms' phi to zero
    for(size_t i = 0; i < nSites; i++)
      xAtoms[i].F = 0.0;
    }

  // Check for consistency in the atoms
  double xSumPhi = 0.0;
  for(size_t i = 0; i < nSites; i++)
    xSumPhi += xAtoms[i].F;
  cout << "READING; SubPhi = " << xSumPhi << endl;
  cout << xAtoms[123].F << endl;
  cout << xAtoms[352].F << endl;

  cout << "UGRID " << uGrid << endl;
  cout << "VGRID " << vGrid << endl;

  // Solve the PDE
  vnl_vector<double> hint = this->GetHintArray();
  try
    {
    this->ComputeAtoms(false);
    }
  catch(...)
    {
    cout << "EXCPETION" << endl;
    }
}

void CartesianMedialModel
::WriteToRegistry(Registry &R)
{
  // Store the type of the m-rep specification
  R["Grid.Type"] << "Cartesian";
  R["Grid.Size.U"] << m;
  R["Grid.Size.V"] << n;

  // Store the grid points (this is not required to read, needed for
  // non-uniform grids)
  R.Folder("Grid.Spacing.U").PutArray(m, uGrid.data_block());
  R.Folder("Grid.Spacing.V").PutArray(n, vGrid.data_block());

  // Store the phi values computed for this mpde
  Mat xPhi(m, n);
  for(size_t i = 0; i < m; i++) for(size_t j = 0; j < n; j++)
    xPhi[i][j] = xAtoms[GetGridAtomIndex(i,j)].F;

  R["Grid.PhiAvailable"] << true;
  R.Folder("Grid.Phi").PutArray(xPhi.size(), xPhi.data_block());

  // Store the fourier coefficient information
  R["SurfaceModel"] << "Fourier";
  xSurface->SaveToRegistry(R.Folder("Fourier"));
  
  // Check for consistency in the atoms
  double xSumPhi = 0.0;
  for(size_t i = 0; i < nSites; i++)
    xSumPhi += xAtoms[i].F;
  cout << "WRITING; SubPhi = " << xSumPhi << endl;
}

