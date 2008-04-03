#include "MedialPDESolver.h"
#include <cmath>
#include <algorithm>
#include "smlmath.h"

#include "MedialPDEMasks.h"

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

MedialPDESolver
::MedialPDESolver(size_t nu, size_t nv, double xScale, size_t pu, size_t pv)
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

MedialPDESolver
::MedialPDESolver(const Vec &uGrid, const Vec &vGrid)
{
  Initialize(uGrid, vGrid);
}

void 
MedialPDESolver
::Initialize(const Vec &uGrid, const Vec &vGrid)
{
  // Copy the grid vectors
  this->uGrid = uGrid; this->vGrid = vGrid;

  // Record the size of the grid
  m = uGrid.size(); n = vGrid.size();

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
      if(ku < 0 || ku >= m || kv < 0 || kv >= n) 
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
  xGrid = new CartesianMedialAtomGrid(m, n);

  // Compute the initial solution
  SetDefaultInitialGuess(1);
}

void
MedialPDESolver
::SetDefaultInitialGuess(double xMagnitude)
{
  // Set the initial guess to default values provided by the sites
  xInitSoln = xMagnitude * xDefaultInitSoln;
}

void 
MedialPDESolver
::SetSolutionAsInitialGuess()
{
  // Make the current solution the initial guess. This makes Newton converge
  // quickly when exploring near a point
  xInitSoln = y;
  flagReuseLastSolution = true;
}

void 
MedialPDESolver
::SetMedialSurface(IMutableHyperSurface2D *xSurface)
{
  // Remember the surface
  this->xSurface = xSurface;

  // Tell the surface where it will be evaluated
  xSurface->SetEvaluationGrid(uGrid, vGrid);
}

void MedialPDESolver::InitializeSiteGeometry()
{
  size_t i, j;
  
  // Initialize each site with the current surface properties
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    unsigned int iGrid = xGrid->GetAtomIndex(i, j);
    unsigned int iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &xAtom = xAtoms[iGrid];

    // Set the atoms' domain coordinates
    xAtom.u = uGrid[i]; xAtom.v = vGrid[j];
    xAtom.uIndex = i; xAtom.vIndex = j;
    double u = uGrid[i], v = vGrid[j];
    
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

bool
MedialPDESolver
::ReconstructAtoms(const Mat &ySolution)
{
  // Keep track of invalid atoms
  bool flagAllAtomsValid = true;
  
  // Once the iterations are complete, reconstruct the atoms
  for(unsigned int i = 0; i < m; i++) for(unsigned int j = 0; j < n; j++)
    {
    // Map to a grid index and to our own index
    unsigned int iGrid = xGrid->GetAtomIndex(i, j);
    unsigned int iLocal = xSiteIndex[i][j];
    
    // The medial atom to update
    MedialAtom &xAtom = xAtoms[iGrid];

    // The case where phi is negative is undefined
    if( ySolution[i][j] > 0 )
      {
      // Compute the derivatives of R using finite differences
      double xTest = xSites[iLocal]->ComputeEquation(ySolution);
      xAtom.F = 
        xMasks[iLocal]->ComputeOneJet(ySolution.data_block(), xAtom.Fu, xAtom.Fv);

      // Compute the boundary properties of the medial point
      xAtom.ComputeBoundaryAtoms(xSites[iLocal]->IsBorderSite());
      }
    else
      {
      // What are we supposed to do?
      // cout << "Negative F at " << xAtom.u << ", " << xAtom.v << endl;
      xAtom.R = xAtom.F = xAtom.Fu = xAtom.Fv = 0.0;
      xAtom.flagValid = false;
      }

    // Update the valid atoms flag
    flagAllAtomsValid &= xAtom.flagValid;
    }

  // Return value: whether all atoms are valid
  return flagAllAtomsValid;
}

double ArrayMinMax(double *array, size_t n, double &xMin, double &xMax)
{
	xMax = 0, xMin = 1e100;
  for(size_t q = 0; q < n; q++)
  	{
    if(xMax < fabs(array[q])) xMax = fabs(array[q]);
    if(xMin > fabs(array[q])) xMin = fabs(array[q]);
    }
}

/** This computes the right hand side of the site equations, given phi=x */
double MedialPDESolver::ComputeNewtonRHS(const Mat& x, Mat &b)
{
  size_t i, j;
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    b[i][j] = - xSites[xSiteIndex[i][j]]->ComputeEquation(x);
  return dot_product(b, b);
}

bool MedialPDESolver::SolveOnce(const Mat &xGuess, double delta)
{
  size_t i, j, k, iIter;
  double epsMax, bMax, bMagSqr;
  
  // Initialize epsilon to zero
  eps.fill(0.0);
  
  // Copy the initial solution to the current solution
  y = xGuess;

  // Compute the right hand side, put it into b
  bMagSqr = ComputeNewtonRHS(y, b);

  // We are now ready to perform the Newton loop
  for(iIter = 0; iIter < 20; iIter++)
    {
    // Compute the Jacobian matrix
    for(size_t iSite = 0; iSite < nSites; iSite++)
      xSites[iSite]->
        ComputeDerivative(y, xSparseValues + xRowIndex[iSite] - 1, iSite+1);
      
    // Perform the symbolic factorization only for the first iteration
    tSolver.Start();
    if(iIter == 0)
      xPardiso.SymbolicFactorization(nSites, xRowIndex, xColIndex, xSparseValues);

    // Compute the Jacobian inverse
    xPardiso.NumericFactorization(xSparseValues);
    xPardiso.Solve(b.data_block(), eps.data_block());
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
    return false;
    }
  else
    {
    // cout << endl;
    return true;
    }
}

bool
MedialPDESolver
::Solve(const Mat &xGuess, double delta)
{
  // Intialize the sites
  InitializeSiteGeometry();

  // Just solve once - no initialization tricks!
  bool flagSolved = SolveOnce(xGuess, delta);
  
  // Reconstruct the medial atoms
  bool flagValid = ReconstructAtoms(y);

  return flagValid && flagSolved;
}

struct AtomGradientTerms
{
  double x1_2R, Ru, Rv, x1_2F, Ru_R, Rv_R, Ru_2F, Rv_2F;
  double g1iRi, g2iRi;
  SMLVec3d N_2g, N_2nt, Xu_aelt, Xv_aelt;
};

void ComputeMedialAtomBoundaryDerivativeCommonTerms(
  MedialAtom *xAtom, AtomGradientTerms &agt)
{
  // Get the relevant elements of the atoms
  SMLVec3d &X = xAtom->X, &Xu = xAtom->Xu, &Xv = xAtom->Xv;
  
  // Get the elements of the first fundamental form and its derivative
  double &g11 = xAtom->G.xContravariantTensor[0][0];
  double &g12 = xAtom->G.xContravariantTensor[0][1];
  double &g22 = xAtom->G.xContravariantTensor[1][1];

  // Get the g's
  double &g = xAtom->G.g; 

  // Get the partials of Phi and its variational derivative
  double &F = xAtom->F, &Fu = xAtom->Fu, &Fv = xAtom->Fv;

  // Get the derivatives of R
  double &R = xAtom->R; 
  agt.x1_2R = 0.5 / R;
  agt.Ru = Fu * agt.x1_2R, agt.Rv = Fv * agt.x1_2R;
  agt.Ru_R = agt.Ru / R;
  agt.Rv_R = agt.Rv / R;
  agt.Ru_2F = 0.5 * agt.Ru / F;
  agt.Rv_2F = 0.5 * agt.Rv / F;
  
  // Terms used to compute the derivative of the normal vector
  agt.Xu_aelt = Xu / xAtom->aelt; agt.Xv_aelt = Xv / xAtom->aelt;
  agt.N_2g = 0.5 * xAtom->N / g;

  // We will compute several intermediate terms
  agt.g1iRi = g11 * agt.Ru + g12 * agt.Rv;
  agt.g2iRi = g12 * agt.Ru + g22 * agt.Rv;

  // Compute the plus-minus term
  agt.N_2nt = 0.5 * xAtom->N / xAtom->xNormalFactor;
}

void ComputeMedialAtomBoundaryDerivative(
  MedialAtom *xAtom, MedialAtom *dAtom, AtomGradientTerms &agt, bool isEdge)
{
  // Get the relevant elements of the atoms
  SMLVec3d &X = xAtom->X, &Xu = xAtom->Xu, &Xv = xAtom->Xv;
  SMLVec3d &Y = dAtom->X, &Yu = dAtom->Xu, &Yv = dAtom->Xv;
  
  // Get the elements of the first fundamental form and its derivative
  double &g11 = xAtom->G.xContravariantTensor[0][0];
  double &g12 = xAtom->G.xContravariantTensor[0][1];
  double &g22 = xAtom->G.xContravariantTensor[1][1];
  double &z11 = dAtom->G.xContravariantTensor[0][0];
  double &z12 = dAtom->G.xContravariantTensor[0][1];
  double &z22 = dAtom->G.xContravariantTensor[1][1];

  // Get the g's
  double &g = xAtom->G.g; double &z = dAtom->G.g;

  // Get the partials of Phi and its variational derivative
  double &F = xAtom->F, &Fu = xAtom->Fu, &Fv = xAtom->Fv;
  double &H = dAtom->F, &Hu = dAtom->Fu, &Hv = dAtom->Fv;

  // Get the derivatives of R
  double &R = xAtom->R;
  double P = H * agt.x1_2R;
  double Pu = Hu * agt.x1_2R - H * agt.Ru_2F;
  double Pv = Hv * agt.x1_2R - H * agt.Rv_2F;
  
  // This is the derivative of the normal vector
  dAtom->N =  vnl_cross_3d(agt.Xu_aelt, Yv);
  dAtom->N += vnl_cross_3d(Yu, agt.Xv_aelt);
  vmuladd(dAtom->N, agt.N_2g, -z);

  // We will compute several intermediate terms
  double z1iRi = z11 * agt.Ru + z12 * agt.Rv;
  double z2iRi = z12 * agt.Ru + z22 * agt.Rv;

  // Compute the derivative of grad R
  vmulset(dAtom->xGradR, Xu, z1iRi + g11 * Pu + g12 * Pv);
  vmuladd(dAtom->xGradR, Xv, z2iRi + g12 * Pu + g22 * Pv);
  vmuladd(dAtom->xGradR, Yu, agt.g1iRi);
  vmuladd(dAtom->xGradR, Yv, agt.g2iRi);

  // Address the edge case first
  if(isEdge) 
    {
    // The normal term vanishes
    dAtom->xBnd[0].N = dAtom->xBnd[1].N = - dAtom->xGradR;
    dAtom->xBnd[1].X = dAtom->xBnd[0].X = 
      Y + R * dAtom->xBnd[0].N + P * xAtom->xBnd[0].N;
    dAtom->xGradRMagSqr = 0.0;
    dAtom->xNormalFactor = 0.0;
    return;
    }

  // Compute the derivative of grad R . grad R
  dAtom->xGradRMagSqr = z1iRi * agt.Ru + z2iRi * agt.Rv 
    + 2.0 * (agt.g1iRi * Pu + agt.g2iRi * Pv);

  // Compute the plus-minus term
  SMLVec3d dNormalTerm;
  vmulset(dNormalTerm, dAtom->N, xAtom->xNormalFactor);
  vmuladd(dNormalTerm, agt.N_2nt, -dAtom->xGradRMagSqr);
  
  // Compute the boundary atom normals
  dAtom->xBnd[0].N = dAtom->xBnd[1].N = - dAtom->xGradR;
  dAtom->xBnd[0].N -= dNormalTerm;
  dAtom->xBnd[1].N += dNormalTerm;

  // Compute the boundary atoms
  dAtom->xBnd[0].X = dAtom->xBnd[1].X = Y;
  vmuladd(dAtom->xBnd[0].X, dAtom->xBnd[0].N, R);
  vmuladd(dAtom->xBnd[0].X, xAtom->xBnd[0].N, P);
  vmuladd(dAtom->xBnd[1].X, dAtom->xBnd[1].N, R);
  vmuladd(dAtom->xBnd[1].X, xAtom->xBnd[1].N, P);
}

void
MedialPDESolver
::PrepareAtomsForVariationalDerivative(IHyperSurface2D *xVariation, MedialAtom *dAtoms)
{
  // Evaluate the variation over the atom grid
  for(size_t i = 0; i < m; i++) for(size_t j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = xGrid->GetAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Access the medial atom underneath
    MedialAtom &dAtom = dAtoms[iGrid];

    // Set the atoms' domain coordinates
    dAtom.u = uGrid[i]; dAtom.v = vGrid[j];
    dAtom.uIndex = i; dAtom.vIndex = j;

    // Evaluate the variation and its derivatives
    xVariation->EvaluateAtGridIndex(i, j, 0, 0, 0, 3, dAtom.X.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 1, 0, 0, 3, dAtom.Xu.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 0, 1, 0, 3, dAtom.Xv.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 2, 0, 0, 3, dAtom.Xuu.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 1, 1, 0, 3, dAtom.Xuv.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 0, 2, 0, 3, dAtom.Xvv.data_block());
    xVariation->EvaluateAtGridIndex(i, j, 0, 0, 3, 1, &dAtom.xLapR);
    }
}

void
MedialPDESolver
::ComputeVariationalGradient(vector<MedialAtom *> &dAtomArray)
{
  size_t i, j, k, N = dAtomArray.size();

  // Create an array of right-hand sides
  Mat rhs(N, nSites), soln(N, nSites);

  // Create an array of helper structs for atom computation
  vector<AtomGradientTerms> agt(nSites);

  // Prepare for the gradient computation
  double t0 = clock();
  for(i = 0; i < m; i++) for(j = 0; j < n; j++)
    {
    // Get the index of the site
    size_t iGrid = xGrid->GetAtomIndex(i, j);
    size_t iSite = xSiteIndex[i][j];

    // Get the medial atom
    MedialAtom &xAtom = xAtoms[iGrid];

    // Initialize the atom gradient terms
    ComputeMedialAtomBoundaryDerivativeCommonTerms(&xAtom, agt[iSite]);

    // Compute the matrix for linear solver
    xSites[iSite]->ComputeVariationalDerivativeMatrix(
      y, xSparseValues + xRowIndex[iSite] - 1, &xAtom);

    for(k = 0; k < N; k++)
      {
      // Compute the derivatives of the metric tensor
      xAtom.ComputeMetricTensorDerivatives(dAtomArray[k][iGrid]);

      // Compute each of the right hand sides
      rhs[k][iSite] = xSites[iSite]->ComputeVariationalDerivativeRHS(
        y, &xAtom, &dAtomArray[k][iGrid]);
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
      size_t iGrid = xGrid->GetAtomIndex(i, j);
      size_t iSite = xSiteIndex[i][j];

      // Access the medial atom underneath
      MedialAtom &dAtom = dAtomArray[k][iGrid];

      // Compute the gradient of phi for the new atom
      dAtom.F = xMasks[iSite]->ComputeOneJet(phi.data_block(), dAtom.Fu, dAtom.Fv);

      // Compute the rest of the atom derivative
      ComputeMedialAtomBoundaryDerivative(
        &xAtoms[iGrid], &dAtom, agt[iSite], xSites[iSite]->IsBorderSite());
      }
    }
  cout << " [BND: " << (clock() - t0) / CLOCKS_PER_SEC << " s] " << flush;
}

double fnTestF01(double u, double v)
{
  return exp(u + v);
}

void fnTestF01Jet(double u, double v, double jet[])
{
  jet[0] = exp(u + v);
  jet[1] = jet[2] = jet[3] = jet[4] = jet[5] = jet[0];
}

void MedialPDESolver
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

  cout << "MedialPDESolver FD Test: " << xMaxDiff << endl;
}
