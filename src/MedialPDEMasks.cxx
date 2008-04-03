#include <iostream>
#include "MedialPDEMasks.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_qr.h"

#include <algorithm>

using namespace std;

template<class T>
vector<size_t> SortTwoArrays(vector<T> &a1, vector<T> &a2)
{
  // Gotta be same size
  assert(a1.size() == a2.size());
  size_t i, n = a1.size();
  
  // Create a joint structure
  typedef pair<T, T> DataPair;
  typedef pair<DataPair, size_t> IndexPair;

  // Create a list of index pairs
  vector<IndexPair> vSort(n);
  vector<size_t> vOutput(n);
  for(i = 0; i < n; i++)
    {
    vSort[i].first.first = a1[i];
    vSort[i].first.second = a2[i];
    vSort[i].second = i;
    }

  // Sort the array based on first two indices
  sort(vSort.begin(), vSort.end());

  // Update the original arrays
  for(i = 0; i < n; i++)
    {
    a1[i] = vSort[i].first.first;
    a2[i] = vSort[i].first.second;
    vOutput[i] = vSort[i].second;
    }

  return vOutput;
}

const size_t FiniteDifferenceMask::NUM_WEIGHTS = 6;

FiniteDifferenceMask::FiniteDifferenceMask(size_t nNodes)
: n(nNodes), W(nNodes, NUM_WEIGHTS, 0.0), flagTransposed(false)
{ 
  qu.resize(nNodes); 
  qv.resize(nNodes); 
  for(size_t i = 0; i < NUM_WEIGHTS; i++)
    {
    fwWeight[i] = NULL;
    fwIndex[i] = NULL;
    fwCount[i] = 0;
    }
}

void FiniteDifferenceMask
::PrintReport()
{
  // Describe how the jet is computed
  const char *strJet[] = {"F","Fu","Fv","Fuu","Fvv","Fuv"};
  for(size_t i = F00; i <= F11; i++)
    {
    bool flagNeedPlus = false;
    cout << strJet[i] << " = ";
    for(size_t j = 0; j < n; j++)
      {
      if(fabs(W[j][i]) > 1.0e-10) 
        {
        if(flagNeedPlus) cout << " + ";
        cout << W[j][i] << " * F[" << qu[j] << "," << qv[j] << "]";
        flagNeedPlus = true; 
        }
      }
    cout << endl;
    }
}

void FiniteDifferenceMask
::FlipMask(bool uFlip, bool vFlip)
{
  size_t i;

  if(uFlip)
    for(i = 0; i < n; i++)
      qu[i] = - qu[i];

  if(vFlip)
    for(i = 0; i < n; i++)
      qv[i] = - qv[i];
}

void FiniteDifferenceMask
::TransposeMask()
{
  // Swap the offsets between u and v
  for(size_t i = 0; i < n; i++)
    {
    int ku = qu[i], kv = qv[i];
    qu[i] = kv; qv[i] = ku;
    }

  // Record the transposed state
  flagTransposed = !flagTransposed;
}

FiniteDifferenceMask
::~FiniteDifferenceMask()
{
  for(size_t i = 0; i < NUM_WEIGHTS; i++)
    {
    if(fwCount[i])
      {
      delete fwIndex[i];
      delete fwWeight[i];
      }
    }
}

void
FiniteDifferenceMask
::OptimizeWeights(size_t uNodes, size_t vNodes)
{
  // Optimize weights for one/two jet computation
  for(size_t i = 0; i < NUM_WEIGHTS; i++)
    {
    // Clean up existing weights
    size_t j;
    if(fwCount[i])
      {
      delete fwIndex[i];
      delete fwWeight[i];
      fwCount[i] = 0;
      }

    // Count non-zero weights
    for(j = 0; j < n; j++)
      if(W[j][i] != 0.0)
        fwCount[i]++;

    if(fwCount[i])
      {
      fwIndex[i] = new size_t[fwCount[i]];
      fwWeight[i] = new double[fwCount[i]];

      // Assign the weights
      size_t k = 0;
      for(size_t j = 0; j < n; j++)
        if(W[j][i] != 0.0)
          {
          fwIndex[i][k] = vNodes * (iu + qu[j]) + iv + qv[j];
          fwWeight[i][k] = W[j][i];
          k++;
          }
      }
    }
  
  // Set the raw position
  iraw = vNodes * (iu) + iv;
}

double FiniteDifferenceMask
::ComputeOneJet(const double *F, double &Fu, double &Fv)
{
  Fu = 0; Fv = 0;

  size_t i;

  for(i = 0; i < fwCount[F10]; i++)
    Fu += fwWeight[F10][i] * F[fwIndex[F10][i]];
  for(i = 0; i < fwCount[F01]; i++)
    Fv += fwWeight[F01][i] * F[fwIndex[F01][i]];

  // Return the value of Y
  return F[iraw];
}

double FiniteDifferenceMask
::ComputeTwoJet(const Mat &F, double &Fu, double &Fv,
  double &Fuu, double &Fuv, double &Fvv)
{
  Fu = 0; Fv = 0; Fuu = 0; Fuv = 0; Fvv = 0;
  for(size_t i = 0; i < n; i++)
    {
    double fi = F[iu + qu[i]][iv + qv[i]];
    Fu  += W[i][F10] * fi;
    Fv  += W[i][F01] * fi;
    Fuu += W[i][F20] * fi;
    Fuv += W[i][F11] * fi;
    Fvv += W[i][F02] * fi;
    }

  // Return the value of Y
  return F[iu][iv];
}

double 
FiniteDifferenceMask
::GetNodeComponentInOneJet(size_t iNode, double &Wu, double &Wv)
{
  Wu = W[iNode][F10];
  Wv = W[iNode][F01];
  return W[iNode][F00];
}

double 
FiniteDifferenceMask
::GetNodeComponentInTwoJet(size_t iNode, 
  double &Wu, double &Wv, double &Wuu, double &Wuv, double &Wvv)
{
  Wu = W[iNode][F10];
  Wv = W[iNode][F01];
  Wuu = W[iNode][F20];
  Wuv = W[iNode][F11];
  Wvv = W[iNode][F02];
  return W[iNode][F00];
}

void 
FiniteDifferenceMask
::SortNodes()
{
  SortTwoArrays(qu, qv);
}

UniformFDMask
::UniformFDMask()
  : FiniteDifferenceMask(9)
{
  // Initialize the pointers
  qu[0] = -1; qv[0] = -1; 
  qu[1] = -1; qv[1] =  0; 
  qu[2] = -1; qv[2] =  1; 
  qu[3] =  0; qv[3] = -1; 
  qu[4] =  0; qv[4] =  0; 
  qu[5] =  0; qv[5] =  1; 
  qu[6] =  1; qv[6] = -1; 
  qu[7] =  1; qv[7] =  0; 
  qu[8] =  1; qv[8] =  1; 
}

void UniformFDMask::ComputeWeights(const double *uGrid, const double *vGrid)
{
  // Get the du and dv quantities
  du = uGrid[iu+1] - uGrid[iu];
  dv = vGrid[iv+1] - vGrid[iv];
  
  // Fu and Fv
  W[7][F10] = 0.5 / du; W[1][F10] = -0.5 / du;
  W[5][F01] = 0.5 / dv; W[3][F01] = -0.5 / dv;

  // Fuu and Fvv
  W[7][F20] = W[1][F20] = 1.0 / (du * du); W[4][F20] = -2.0 / (du * du);
  W[5][F02] = W[3][F02] = 1.0 / (dv * dv); W[4][F02] = -2.0 / (dv * dv);

  // Fuv
  W[8][F11] = W[0][F11] = 0.25 / (du * dv); 
  W[2][F11] = W[6][F11] = - W[0][F11];
}

DoublyNonuniformFDMask
::DoublyNonuniformFDMask()
  : FiniteDifferenceMask(10)
{
  qu[0] =  0;  qv[0] =  0;
  qu[1] =  1;  qv[1] =  0;
  qu[2] =  0;  qv[2] =  1;
  qu[3] = -1;  qv[3] =  1;
  qu[4] = -1;  qv[4] =  0;
  qu[5] = -1;  qv[5] = -1;
  qu[6] =  0;  qv[6] = -1;
  qu[7] =  1;  qv[7] = -1;
  qu[8] =  2;  qv[8] =  0;
  qu[9] =  0;  qv[9] =  2;
}

/*
void
SinglyNonuniformFDMask
::SetLocation(size_t u, size_t v)
{
  qu[0] = u +  0;  qv[0] = v +  0;
  qu[1] = u +  1;  qv[1] = v +  0;
  qu[2] = u +  1;  qv[2] = v +  1;
  qu[3] = u +  0;  qv[3] = v +  1;
  qu[4] = u + -1;  qv[4] = v +  1;
  qu[5] = u + -1;  qv[5] = v +  0;
  qu[6] = u + -1;  qv[6] = v + -1;
  qu[7] = u +  0;  qv[7] = v + -1;
  qu[8] = u +  1;  qv[8] = v + -1;
  qu[9] = u +  2;  qv[9] = v +  0;
}


void SinglyNonuniformFDMask
::ComputeWeights(const double *uGrid, const double *vGrid)
{ 
  // Use the generic formula to compute the weights
  ComputeWeightsGeneric(uGrid, vGrid); 

  // At this point, the weights for Fu / Fvv and Fv / Fvv are going to be 
  // messed up. Use the uniform grid settings there
  if(flagTransposed)
    {
    // Clear the rows of W
    for(size_t i = 0; i < n; i++)
      W[i][F10] = W[i][F20] = 0.0;
    
    // U is uniformly spaced
    double du = uGrid[qu[3]] - uGrid[qu[0]];

    // Set the U partial derivatives
    W[3][F10] = 0.5 / du; W[7][F10] = -0.5 / du;
    W[3][F20] = W[7][F20] = 1.0 / (du * du); W[0][F20] = -2.0 / (du * du);
    }
  else
    {
    // Clear the rows of W
    for(size_t i = 0; i < n; i++)
      W[i][F01] = W[i][F02] = 0.0;
    
    // V is uniformly spaced
    double dv = vGrid[qv[3]] - vGrid[qv[0]];

    // Set the U partial derivatives
    W[3][F01] = 0.5 / dv; W[7][F01] = -0.5 / dv;
    W[3][F02] = W[7][F02] = 1.0 / (dv * dv); W[0][F02] = -2.0 / (dv * dv);
    }
}
*/

void FiniteDifferenceMask
::ComputeWeights(const double *uGrid, const double *vGrid)
{
  // Create a matrix with all the weightings
  Mat A(10, n, 0.0);
  for(size_t i = 0; i < n; i++)
    {
    // Get the grid spacing elements
    double du = uGrid[iu + qu[i]] - uGrid[iu];
    double dv = vGrid[iv + qv[i]] - vGrid[iv];
    
    A[F00][i] = 1;
    A[F10][i] = du;
    A[F01][i] = dv;
    A[F20][i] = 0.5 * du * du;
    A[F02][i] = 0.5 * dv * dv;
    A[F11][i] = du * dv;
    A[F30][i] = 0.1666666666666667 * du * du * du;
    A[F03][i] = 0.1666666666666667 * dv * dv * dv;
    A[F21][i] = 0.5 * du * du * dv;
    A[F12][i] = 0.5 * du * dv * dv;
    }

  // Take the submatrix of size corresponding to the number of nodes and
  // compute its SVD
  vnl_qr<double> qr(A.extract(n, n));

  // Solve the system for the relevant partial derivatives
  for(int p = F00; p < F00 + n && p <= F11; p++)
    {
    // Set up the right hand (solution vector)
    vnl_vector<double> b(n, 0.0); b[p] = 1.0;

    // Solve to get the n weights
    W.set_column(p, qr.solve(b));
    }
}

void CrossFDMask
::ComputeWeights(const double *uGrid, const double *vGrid)
{
  size_t i, iCenter;

  // Generate horizontal and vertical array
  vector<size_t> uSamples, vSamples;
  for(i = 0; i < n; i++)
    {
    if(qu[i] == 0) vSamples.push_back(i);
    if(qv[i] == 0) uSamples.push_back(i);
    if(qu[i] == 0 && qv[i] == 0) iCenter = i;
    } 

  // Get the number of samples in each array
  size_t nu = uSamples.size(), nv = vSamples.size();

  // For each array, construct a matrix
  Mat AU(nu, nu, 0.0);
  for(i = 0; i < nu; i++)
    {
    double du = uGrid[iu + qu[uSamples[i]]] - uGrid[iu];
    double jFact = 1.0, duPower = 1.0;
    for(size_t j = 0; j < nu; j++)
      {
      AU[j][i] = duPower / jFact;
      duPower *= du;
      jFact *= (j+1.0);
      }
    }

  Mat AV(nv, nv, 0.0);
  for(i = 0; i < nv; i++)
    {
    double dv = vGrid[iv + qv[vSamples[i]]] - vGrid[iv];
    double jFact = 1.0, dvPower = 1.0;
    for(size_t j = 0; j < nv; j++)
      {
      AV[j][i] = dvPower / jFact;
      dvPower *= dv;
      jFact *= (j+1.0);
      }
    }

  // Invert the two matrices
  vnl_qr<double> qru(AU), qrv(AV);

  // Solve the system for the weights of the U derivatives
  int uFactors[] = {F00, F10, F20, F30};
  for(i = 1; i < nu; i++)
    {
    // Solve for the weightings
    vnl_vector<double> b(nu, 0.0); b[i] = 1;
    Vec z = qru.solve(b);

    // Stick them into the W matrix
    for(size_t j = 0; j < nu; j++)
      W[uSamples[j]][uFactors[i]] = z[j];
    }
   
  // Solve the system for the weights of the U derivatives
  int vFactors[] = {F01, F01, F02, F03};
  for(i = 1; i < nv; i++)
    {
    // Solve for the weightings
    vnl_vector<double> b(nv, 0.0); b[i] = 1;
    Vec z = qrv.solve(b);

    // Stick them into the W matrix
    for(size_t j = 0; j < nv; j++)
      W[vSamples[j]][vFactors[i]] = z[j];
    }
    
  // Take care of the center node
  W[iCenter][F00] = 1.0;
}

/** 
void DoublyNonuniformFDMask::ComputeWeights(double *uGrid, double *vGrid)
{
  // Set up the matrix with which to compute the weights
  Mat A(10, 10, 0.0);
  for(size_t i = 0; i < 10; i++)
    {
    // Get the grid spacing elements
    double du = uGrid[qu[i]] - uGrid[qu[0]];
    double dv = vGrid[qv[i]] - vGrid[qv[0]];
    
    A[0][i] = 1;
    A[1][i] = du;
    A[2][i] = dv;
    A[3][i] = 0.5 * du * du;
    A[4][i] = du * dv;
    A[5][i] = 0.5 * dv * dv;
    A[6][i] = 0.1666666666666667 * du * du * du;
    A[7][i] = 0.5 * du * du * dv;
    A[8][i] = 0.5 * du * dv * dv;
    A[9][i] = 0.1666666666666667 * dv * dv * dv;
    }

  // Set up the solution vectors
  Mat B(10, 5, 0.0);

  // Set up the right hand sides
  B[1][0] = 1; // Fu
  B[2][1] = 1; // Fv
  B[3][2] = 1; // Fuu
  B[4][3] = 1; // Fuv
  B[5][4] = 1; // Fvv

  // Solve to obtain the weights of the factors
  vnl_svd<double> svd(A);
  W = svd.solve(B);
}
*/

CornerFDMask::CornerFDMask()
  : CrossFDMask(5)
{
  qu[0] = 0; qv[0] = 0;
  qu[1] = 1; qv[1] = 0;
  qu[2] = 0; qv[2] = 1;
  qu[3] = 2; qv[3] = 0;
  qu[4] = 0; qv[4] = 2;
}

BorderFDMask::BorderFDMask()
  : CrossFDMask(5)
{
  qu[0] =  0; qv[0] =  0;
  qu[1] =  1; qv[1] =  0;
  qu[2] =  0; qv[2] =  1;
  qu[3] =  0; qv[3] = -1;
  qu[4] =  2; qv[4] =  0;
}

SimpleCornerFDMask::SimpleCornerFDMask()
  : CrossFDMask(3)
{
  qu[0] = 0;  qv[0] = 0;
  qu[1] = 1;  qv[1] = 0;
  qu[2] = 0;  qv[2] = 1;
}

SimpleBorderFDMask::SimpleBorderFDMask()
  : CrossFDMask(4)
{
  qu[0] = 0;  qv[0] = 0;
  qu[1] = 1;  qv[1] = 0;
  qu[2] = 0;  qv[2] = 1;
  qu[3] = 0;  qv[3] =-1;
}

double TestFunction(double u, double v)
{
  return u*u*u*v*v + 2*v*v*u*u - 3*u*u*v + 5*u*v - 4*u + 3*v + 5;
}

void TestFunctionJet(double u, double v, double jet[6])
{
  jet[0] = u*u*u*v*v + 2*v*v*u*u - 3*u*u*v + 5*u*v - 4*u + 3*v + 5;
  jet[1] = 3*u*u*v*v + 4*v*v*u   - 6*u*v   + 5*v   - 4            ;
  jet[2] = 2*u*u*u*v + 4*v*u*u   - 3*u*u   + 5*u         + 3      ;
  jet[3] = 6*u*v*v   + 4*v*v     - 6*v                            ;
  jet[4] = 2*u*u*u   + 4*u*u                                      ;
  jet[5] = 6*u*u*v   + 8*v*u     - 6*u     + 5                    ;
}

double TestFunction2(double u, double v)
{
  return sin(u + v) + cos(u*u - v*v);
}

void TestFunctionJet2(double u, double v, double jet[6])
{
  jet[0] = sin(u + v) + cos(u*u - v*v);
  jet[1] = cos(u + v) - 2 * u * sin(u*u - v*v);
  jet[2] = cos(u + v) + 2 * v * sin(u*u - v*v);
  jet[3] = - sin(u + v) - 2 * sin(u*u - v*v) - 4 * u * u * cos(u*u - v*v);
  jet[4] = - sin(u + v) + 2 * sin(u*u - v*v) - 4 * v * v * cos(u*u - v*v);
  jet[5] = - sin(u + v) + 4 * v * u * cos(u*u - v*v);
}

typedef vnl_vector<double> Vec;

void ComputeMaskJet(FiniteDifferenceMask *T, const Vec &uGrid, const Vec &vGrid, Vec &jet)
{
  // Set up a sample function
  vnl_matrix<double> F(uGrid.size(), vGrid.size());
  for(size_t iu = 0; iu < uGrid.size(); iu++) 
    for(size_t iv = 0; iv < vGrid.size(); iv++)
      F[iu][iv] = TestFunction2(uGrid[iu], vGrid[iv]);

  // Compute the weights in the site
  T->ComputeWeights(uGrid.data_block(), vGrid.data_block());

  // Compute the jet
  jet[0] = T->ComputeTwoJet(F, jet[1], jet[2], jet[3], jet[5], jet[4]);
}

void TestFDMask(
  FiniteDifferenceMask *T, double x, double y, const Vec &uGrid, const Vec &vGrid)
{
  // Two epsili to test accuracy
  double eps1 = 0.01, eps2 = 0.001;

  // Set the location of the mask to the middle
  T->SetLocation(uGrid.size() / 2, vGrid.size() / 2);

  // Transpose the mask - just testing
  T->TransposeMask();
  T->FlipMask(false, true);

  // Sort mask: for extra testing
  T->SortNodes();
  
  // Shift and scale the grid
  Vec uGrid1 = x + uGrid * eps1, vGrid1 = y + vGrid * eps1;
  Vec uGrid2 = x + uGrid * eps2, vGrid2 = y + vGrid * eps2;

  // cout << "SITE REPORT : " << endl;
  // T1.PrintReport();
  // cout << endl;
  //


  // Compute the jet and compare to correct answer
  vnl_vector<double> jet1(6), jet2(6), jetTrue(6);

  // Compute the first jet
  ComputeMaskJet(T, uGrid1, vGrid1, jet1);
  ComputeMaskJet(T, uGrid2, vGrid2, jet2);

  // Compute the true jet 
  TestFunctionJet2(x, y, jetTrue.data_block());

  // Compute the ratio of the errors for each term
  Vec xDiffRatio = 
    (jet1 - jetTrue).apply(fabs).apply(log10) - 
    (jet2 - jetTrue).apply(fabs).apply(log10);

  // Compute the error order estimate
  Vec alhpa1 = (jet1 - jetTrue).apply(fabs).apply(log) / log(eps1);
  Vec alhpa2 = (jet2 - jetTrue).apply(fabs).apply(log) / log(eps2);

  // Report the error
  cout << "Order of Error (eps1) : " << alhpa1 << endl;
  cout << "Order of Error (eps2) : " << alhpa2 << endl;
  cout << "Convergence Order     : " << xDiffRatio << endl << endl;
}

void TestFDStuff()
{
  int n = 7;
  double xUniformGrid[] = {-3, -2, -1, 0, 1, 2, 3};
  double xNonUniformGrid[] = {-4, -3, -2, 0, 1, 3, 4};
  Vec xUniGrid(xUniformGrid, 6), xNonGrid(xNonUniformGrid, 6);

  // Test the uniform mask
  UniformFDMask mUniform;
  cout << "TESTING: Uniform Mask" << endl;
  TestFDMask(&mUniform, 0.3, 0.6, xUniGrid, xUniGrid);

  // Test the singly nonuniform mask
  // SinglyNonuniformFDMask mSingly;
  // cout << "TESTING: Singly Uniform Mask" << endl;
  // TestFDMask(&mSingly, 0.3, 0.6, xNonGrid, xUniGrid);

  // Test the doubly nonuniform mask
  DoublyNonuniformFDMask mDoubly;
  cout << "TESTING: Doubly Uniform Mask" << endl;
  TestFDMask(&mDoubly, 0.3, 0.6, xNonGrid, xNonGrid);

  // Test the doubly nonuniform mask
  DoublyNonuniformFDMask mDoubly2;
  cout << "TESTING: Doubly Uniform Mask (Mixed Uniformity) " << endl;
  TestFDMask(&mDoubly2, 0.3, 0.6, xNonGrid, xUniGrid);

  // Test the border mask
  BorderFDMask mBorder;
  cout << "TESTING: Border Mask" << endl;
  TestFDMask(&mBorder, 0.3, 0.6, xNonGrid, xNonGrid);

  // Test the border mask
  CornerFDMask mCorner;
  cout << "TESTING: Corner Mask" << endl;
  TestFDMask(&mCorner, 0.3, 0.6, xNonGrid, xNonGrid);

  // Test the border mask
  SimpleBorderFDMask mSimpleBorder;
  cout << "TESTING: SimpleBorder Mask" << endl;
  TestFDMask(&mSimpleBorder, 0.3, 0.6, xUniGrid, xUniGrid);

  // Test the border mask
  SimpleCornerFDMask mSimpleCorner;
  cout << "TESTING: SimpleCorner Mask" << endl;
  TestFDMask(&mSimpleCorner, 0.3, 0.6, xUniGrid, xUniGrid);
}
