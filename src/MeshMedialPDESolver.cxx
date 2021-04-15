#include "MeshMedialPDESolver.h"
#include "MedialAtomGrid.h"
#include <iomanip>
#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_random.h>
#include <algorithm>
#include <fstream>

using namespace std;

#define TRANSFER_R2
#undef TRANSFER_LOGR2

void tensolve_mpde_solve(MeshMedialPDESolver *solver, double *x);

double levenberg_marquardt(
  void *handle,
  int n,
  double *x_init,
  double *f_init,
  ImmutableSparseMatrix<double> &J,
  void (*fnResidual)(void *, int n, double *x, double *fx),
  void (*fnJacobian)(void *, int n, double *x, ImmutableSparseMatrix<double> &J),
  double *x_out,
  int n_iter = 100,
  double eps1 = 1.0e-8,
  double eps2 = 1.0e-12,
  double tau = 1.0e-3)
{
  // Typedefs 
  typedef vnl_vector<double> Vec;
  typedef ImmutableSparseMatrix<double> SMat;
  
  // Initial solution and updating solution
  Vec x(x_init, n);

  // Function value vector
  Vec fx(f_init, n);

  // Squared norm of fx
  double normsq_fx = fx.squared_magnitude();

  // Initial parameters
  size_t i, k = 0;
  double nu = 2;

  // Initialize and compute the matrix A
  SMat A;
  SMat::InitializeATA(A, J);
  SMat::ComputeATA(A, J);

  // We also need a pair of offset row and column arrays for calling pardiso
  vnl_vector<int> xPardisoRowIndex(n + 1);
  vnl_vector<int> xPardisoColIndex(A.GetNumberOfSparseValues());
  for(i = 0; i < xPardisoRowIndex.size(); i++)
    xPardisoRowIndex[i] = 1 + A.GetRowIndex()[i];
  for(i = 0; i < xPardisoColIndex.size(); i++)
    xPardisoColIndex[i] = 1 + A.GetColIndex()[i];

  // Initialize the PARDISO interface
  SparseSolver *xSolver = SparseSolver::MakeSolver(true);
  xSolver->SymbolicFactorization(A);

  // Compute the vector g
  Vec minus_g = -J.MultiplyTransposeByVector(fx);

  // Vector representing the LM step
  Vec h(n, 0.0);

  // Set found to false
  bool found = false;

  // Set mu
  double maxA = 0;
  for(i = 0; i < A.GetNumberOfSparseValues(); i++)
    if(maxA < A.GetSparseData()[i])
      maxA = A.GetSparseData()[i];
  printf("Max(A) = %g, mu = %g\n", maxA, tau * maxA);
  double mu = tau * maxA;

  // Iterate
  while(!found && k++ < (size_t) n_iter)
    {
    // Let A = J' * J + mu I
    for(i = 0; i < (size_t) n; i++)
      A.GetSparseData()[A.GetRowIndex()[i]] += mu;
   
    // Solve the linear system (J' * J + mu I) x = -g
    xSolver->NumericFactorization(A.GetSparseData());
    xSolver->Solve(minus_g.data_block(), h.data_block());
    delete xSolver;

    // Compute the norm of x and h
    double norm_h = h.magnitude();
    double norm_x = x.magnitude();

    // Report the solution
    bool verbose = true;
    if(verbose)
      printf(
        "   -> LM Iter: %3lu    F = %8.4e    fmax = %8.4e    fmean = %8.4e"
        "    mu = %8.4e    "
        "|h| = %8.4e    |x| = %8.4e\n",
        (long unsigned) k, normsq_fx / 2, fx.inf_norm(), 
       dot_product(fx, vnl_vector<double>(n, 1.0/n)),
        mu, norm_h, norm_x); 
    else
      cout << "." << flush;
   
    // Check for convergence
    if(norm_h <= eps2 * (norm_x + eps2)) found = true;
    else
      {
      // Update the x vector
      Vec xnew = x + h;

      // Compute the new function value
      Vec fxnew(n, 0.0); 
      fnResidual(handle, n, xnew.data_block(), fxnew.data_block());

      // Compute J . h + fx
      Vec q = fx + J.MultiplyByVector(h);
      
      // Compute rho
      double normsq_fxnew = fxnew.squared_magnitude();
      double normsq_q = q.squared_magnitude();
      double rho = (normsq_fx - normsq_fxnew) / (normsq_fx - normsq_q);

      // Step update (or not)
      if(rho > 0)
        {
        // Accept the new step
        x = xnew;
        fx = fxnew;
        normsq_fx = normsq_fxnew;

        // Recompute the Jacobian
        fnJacobian(handle, n, x.data_block(), J);

        // Recompute A
        SMat::ComputeATA(A, J);

        // Recompute g
        minus_g = -J.MultiplyTransposeByVector(fx);

        // Get the infinity norm of g
        double norminf_g = minus_g.inf_norm();
        found = (norminf_g <= eps1);

        // Recompute mu
        double tmp = 2.0 * rho - 1.0;
        mu *= std::max(1.0 / 3.0, 1.0 - tmp * tmp * tmp);
        nu = 2;
        }
      else
        {
        mu *= nu;
        nu *= 2.0;
        }
      } // If convergence
    } // iteration loop

  printf("* %8.4e\n", normsq_fx / 2);

  // Copy the solution to x out
  x.copy_out(x_out);

  return normsq_fx / 2;
}

#define SE(x) if(x != NULL) {delete x; x = NULL; }

template <class T>
void reset_ptr(T* &x)
{
  if(x != NULL)
    { delete[] x; x = NULL; }
}

/* 
void
MeshMedialPDESolver
::LMSolve()
{
  // Solves the equation using LM method
  int k = 0;
  double nu = 2;

  // The initial solution
  vnl_vector<double> x0 = xSoln;
}
*/

/*
void 
MeshMedialPDESolver
::ComputeJacobianConditionNumber()
{
  size_t i, j;

  // Compute the inverse
  double xNormInvSqr = 0.0;
  Vec ei(topology->nVertices, 0.0), ai(topology->nVertices, 0.0);
  for(i = 0; i < topology->nVertices; i++)
    {
    // Get a column of the inverse
    ei[i] = 1.0; xPardiso.Solve(ei.data_block(), ai.data_block()); ei[i] = 0.0;

    // Print "+" to indicate a solver step
    xNormInvSqr += ai.squared_magnitude();
    }

  // Compute the norm of the matrix itself
  double xNormSqr = 0.0;
  for(j = 0; j < A.GetNumberOfSparseValues(); j++)
    xNormSqr += A.GetSparseData()[j] * A.GetSparseData()[j];

  // Return the ratio of the norms
  cout << "CONDITION NUMBER: " << sqrt(xNormSqr * xNormInvSqr);
}
*/

void
MeshMedialPDESolver
::SetMeshTopology(MeshLevel *topology, MedialAtom *inputAtoms)
{
  // We must first set the dimensions of the matrix A. The finite difference
  // equations involving the LBO are specified at internal vertices in the
  // mesh and involve all neighbors of the vertex. Equations involving the
  // gradient use Loop's approximation, which will also use the entire
  // neigborhood of the vertex (except the vertex itself, but the finite
  // difference equation does use it)

  // Clean up storage
  Reset();

  // Store the mesh
  this->topology = topology;

  // Shorthand for number of vertices, etc
  size_t n = topology->nVertices;

  // Sparse matrix defining the mesh structure
  typedef TriangleMesh::NeighborMatrix NType;
  const TriangleMesh::NeighborMatrix &Nt = topology->GetNeighborMatrix();

  // This code will generate the sparse matrix for the linear system. The
  // system consists of a block matrix
  // [ A  U1 ] [ phi   ] = [ 0   ]
  // [ U2 0  ] [ omega ]   [ g1  ]
  // [ 0  A  ]             [ rho ]
  // [ N  0  ]             [ g2  ]

  // Source matrix (STL), is 2n by 2n
  SparseMat::STLSourceType S(2 * n);

  // Loop over the vertices
  for(size_t i = 0; i < n; i++)
    {
    // Create a list of all the neighbors and the vertex itself
    vector<int> ring;

    // Add the vertex itself to the ring
    ring.push_back(i);
    
    // Add each neighbor
    for(NType::ConstRowIterator row = Nt.Row(i); !row.IsAtEnd(); ++row)
      ring.push_back(row.Column());

    // Sort the list, because Pardiso requires it to be that way
    sort(ring.begin(), ring.end());

    if(topology->IsVertexInternal(i))
      {
      for(vector<int>::iterator ir = ring.begin(); ir != ring.end(); ++ir)
        {
        size_t j = *ir;

        // Add entry in A(phi)
        S[i].push_back(make_pair(j, 0.0));

        // Add entry in A(omega)
        S[n+i].push_back(make_pair(n+j, 0.0));
        }                     
      S[i].push_back(make_pair(n + i, -1.0));
      }
    else
      {
      for(vector<int>::iterator ir = ring.begin(); ir != ring.end(); ++ir)
        {
        size_t j = *ir;

        // Add entry in N
        S[i+n].push_back(make_pair(j, 0.0));
        }
      S[i].push_back(make_pair(i, 1.0));
      }
    }

  // Now create the sparse matrix
  M.SetFromSTL(S, 2 * n);

  // Perform PARDISO symbolic factorization
  // xPardiso.SymbolicFactorization(M);

  // For each of the 'interesting' quadrants in M (A-phi, A-omega, N)
  // we need to create an index mapping the topology neighborhood 
  // structure into the sparse matrix.
  xIndexAPhi.Initialize(Nt.GetNumberOfRows(), Nt.GetNumberOfSparseValues());
  xIndexAOmega.Initialize(Nt.GetNumberOfRows(), Nt.GetNumberOfSparseValues());
  xIndexN.Initialize(Nt.GetNumberOfRows(), Nt.GetNumberOfSparseValues());

  for(size_t i = 0; i < n; i++)
    {
    // Go over the row in the A-phi quadrant
    for(SparseMat::RowIterator rM = M.Row(i); !rM.IsAtEnd(); ++rM)
      {
      size_t j = rM.Column();
      if(i == j)
        xIndexAPhi.xSelfIndex[i] = rM.SparseIndex();
      else if(j < n)
        xIndexAPhi.xNbrIndex[Nt.FindEntryIndex(i, j)] = rM.SparseIndex();
      }

    // Go over the row in the bottom of the matrix
    for(SparseMat::RowIterator rM = M.Row(i+n); !rM.IsAtEnd(); ++rM)
      {
      size_t j = rM.Column();

      // We are in the N part
      if(j < n)
        {
        if(i == j)
          xIndexN.xSelfIndex[i] = rM.SparseIndex();
        else
          xIndexN.xNbrIndex[Nt.FindEntryIndex(i, j)] = rM.SparseIndex();
        }
      // We are in the A-omega part
      else
        {
        if(i == j - n)
          xIndexAOmega.xSelfIndex[i] = rM.SparseIndex();
        else
          xIndexAOmega.xNbrIndex[Nt.FindEntryIndex(i, j - n)] = rM.SparseIndex();
        }
      }
    }

  // Initialize right hand side and solution arrays
  xRHS.set_size(2 * n);
  xSolution.set_size(2 * n);

  // Initialize the triangle area array and other such arrays
  xTriangleGeom = new TriangleGeom[topology->triangles.size()];
  xVertexGeom = new VertexGeom[topology->nVertices];

  // Initialize the tangent weights
  xLoopScheme.SetMesh(topology);

  // Copy the pointer to the atoms. The atoms are passed in to this method
  // and are managed by the parent
  xAtoms = inputAtoms;

  // Initialize the weight arrays for gradient computation
  WX = new SMLVec3d[M.GetNumberOfSparseValues()];
  Wfu = new NeumannDerivWeights[n];

  // Allocate vector for derivative computation
  xTempDerivativeTerms = TempDerivativeTermsArray(topology->nVertices);
}

MeshMedialPDESolver
::MeshMedialPDESolver()
{
  xTriangleGeom = NULL;
  xVertexGeom = NULL;
  topology = NULL;
  xAtoms = NULL;
  WX = NULL;
  Wfu = NULL;
  xSolver = SparseSolver::MakeSolver(false);
}

MeshMedialPDESolver
::~MeshMedialPDESolver()
{
  Reset();
  delete xSolver;
}

void
MeshMedialPDESolver
::Reset()
{
  reset_ptr(xTriangleGeom);
  reset_ptr(xVertexGeom);
  reset_ptr(WX);
  reset_ptr(Wfu);
}

void
MeshMedialPDESolver
::ComputeMeshGeometry(bool flagGradient)
{
  size_t i;
  
  // Clear the vertex fan values that will be accumulated later on
  for(i = 0; i < topology->nVertices; i++)
    xVertexGeom[i].xFanArea = 0.0;
  
  // First, we precompute some triangle-wise measurements, specifically the
  // area of each of the mesh triangles. Unfortunately, this requires sqrt...
  for(i = 0; i < topology->triangles.size(); i++)
    {
    // References to the triangle data
    TriangleGeom &TG = xTriangleGeom[i];
    
    // Get the triangle and its three vertices
    Triangle &t = topology->triangles[i];
    const SMLVec3d &A = xAtoms[t.vertices[0]].X;
    const SMLVec3d &B = xAtoms[t.vertices[1]].X;
    const SMLVec3d &C = xAtoms[t.vertices[2]].X;

    // Compute edge vectors
    SMLVec3d BC = B - C, CA = C - A, AB = A - B;

    // Get the squared lengths of the three segments
    double a = BC.squared_magnitude();
    double b = CA.squared_magnitude();
    double c = AB.squared_magnitude();
    
    // Compute the area of the triange (use lengths)
    TG.xArea = 0.25 * sqrt(
      (a + a - b) * b + (b + b - c) * c + (c + c - a) * a); 

    // Compute the terms of the form (a + b - c)
    double qc = a + b - c, qb = a + c - b, qa = b + c - a;

    // Compute the cotangent of each angle
    double xOneOverFourArea = 0.25 / TG.xArea; 
    TG.xCotangent[0] = xOneOverFourArea * qa;
    TG.xCotangent[1] = xOneOverFourArea * qb;
    TG.xCotangent[2] = xOneOverFourArea * qc;

    // Add the weights to the fan-weight array
    xVertexGeom[t.vertices[0]].xFanArea += xTriangleGeom[i].xArea;
    xVertexGeom[t.vertices[1]].xFanArea += xTriangleGeom[i].xArea;
    xVertexGeom[t.vertices[2]].xFanArea += xTriangleGeom[i].xArea;
      
    // Now compute the derivatives of these terms with respect to each vertex
    if(flagGradient) 
      {
      // Compute the terms of the form qa * (B - C)
      SMLVec3d pa = qa * BC, pb = qb * CA, pc = qc * AB;

      // Compute the derivatives of the area
      double xOneOverEightArea = 0.5 * xOneOverFourArea;
      TG.xAreaGrad[0] = xOneOverEightArea * (pc - pb);
      TG.xAreaGrad[1] = xOneOverEightArea * (pa - pc);
      TG.xAreaGrad[2] = xOneOverEightArea * (pb - pa);

      // Compute intermediates for the derivatives of the cotangent
      double xOneOverTwoArea = xOneOverFourArea + xOneOverFourArea;
      SMLVec3d wAB = AB * xOneOverTwoArea;
      SMLVec3d wBC = BC * xOneOverTwoArea;
      SMLVec3d wCA = CA * xOneOverTwoArea;

      // Compute quantities Cot[A] / Area
      double xOneOverArea = xOneOverTwoArea + xOneOverTwoArea;
      double xCOA_A = TG.xCotangent[0] * xOneOverArea;
      double xCOA_B = TG.xCotangent[1] * xOneOverArea;
      double xCOA_C = TG.xCotangent[2] * xOneOverArea;

      // Compute the cotangent derivatives
      TG.xCotGrad[0][0] =  wAB - wCA - TG.xAreaGrad[0] * xCOA_A;
      TG.xCotGrad[0][1] =        wCA - TG.xAreaGrad[1] * xCOA_A;
      TG.xCotGrad[0][2] = -wAB       - TG.xAreaGrad[2] * xCOA_A;
      
      TG.xCotGrad[1][0] = -wBC       - TG.xAreaGrad[0] * xCOA_B;
      TG.xCotGrad[1][1] =  wBC - wAB - TG.xAreaGrad[1] * xCOA_B;
      TG.xCotGrad[1][2] =        wAB - TG.xAreaGrad[2] * xCOA_B;

      TG.xCotGrad[2][0] =        wBC - TG.xAreaGrad[0] * xCOA_C;
      TG.xCotGrad[2][1] = -wCA       - TG.xAreaGrad[1] * xCOA_C;
      TG.xCotGrad[2][2] =  wCA - wBC - TG.xAreaGrad[2] * xCOA_C;
      }

    // cout << "Triangle geometry at " << i << ": ";
    // cout << "Area = " << xTriangleGeom[i].xArea;
    // cout << " Cot = " << xTriangleGeom[i].xCotangent << endl;
    }

  // Now, compute the tangent vectors for all the points. Tangent vectors are 
  // given by Loop's formula and require iteration around vertices.
  for(i = 0; i < topology->nVertices; i++)
    {
    // Get a reference to the geometry object
    MedialAtom &a = xAtoms[i];

    // Clear the tangent for this vertex
    a.Xu = xLoopScheme.Xu(i, xAtoms);
    a.Xv = xLoopScheme.Xv(i, xAtoms);

    // Set the one get of the geometry descriptor
    a.G.SetOneJet(a.X.data_block(), a.Xu.data_block(), a.Xv.data_block());

    // Print a warning if angle between a.Xu and a.Xv is too small
    double g11 = a.G.xCovariantTensor[0][0];
    double g22 = a.G.xCovariantTensor[1][1];
    double g12 = a.G.xCovariantTensor[0][1];
    double g11_times_g12 = g11 * g22;
    double sin_alpha_sqr = (g11_times_g12 - g12 * g12) / g11_times_g12;
    if(sin_alpha_sqr < 3.0459e-04)
      {
      cerr << "WARNING! Xu and Xv almost parallel in atom " << i << endl;
      }
    
    // Compute the normal vector
    a.ComputeNormalVector();
    } 

#if defined(TRANSFER_R2)

  for(i = 0; i < topology->nVertices; i++)
    if(!topology->IsVertexInternal(i))
      xAtoms[i].F = xAtoms[i].R * xAtoms[i].R;
  for(i = 0; i < topology->nVertices; i++)
    if(!topology->IsVertexInternal(i))
      xAtoms[i].Fv = xLoopScheme.Fv(i, xAtoms);

#elif defined(TRANSFER_LOGR2)

  for(i = 0; i < topology->nVertices; i++)
    if(!topology->IsVertexInternal(i))
      xAtoms[i].F = log(xAtoms[i].R * xAtoms[i].R);

  for(i = 0; i < topology->nVertices; i++)
    if(!topology->IsVertexInternal(i))
      xAtoms[i].Fv = xLoopScheme.Fv(i, xAtoms);

#endif
}

void
MeshMedialPDESolver
::FillSparseMatrix(bool flagInputChange)
{
  // At this point, the structure of the matrix A has been specified. We have
  // to specify the values. This is done one vertex at a time.
  for(size_t i = 0; i < topology->nVertices; i++)
    {
    EdgeWalkAroundVertex it(topology, i);
    if(!it.IsOpen())
      {
      // If there is no input change since the last call (only phi has
      // changed) we don't do any computation here
      if(flagInputChange == false) continue;

      // V is an internal, LBO vertex. It, and each of its neighbors are
      // involved in the finite difference equation. Since this is a
      // differential operator, we can use the fact that all weights must add
      // up to zero 

      // Accumulator for the sum of the weights around the vertex 
      double w_accum = 0.0;

      // The scaling factor applied to cotangent of each angle
      double scale = 1.5 / xVertexGeom[i].xFanArea;

      // Add weight to each vertex
      for( ; !it.IsAtEnd(); ++it)
        {
        // What is the weight of eps_j? Under Desbrun's discretization given
        // in (2.10) of Xu's paper, coefficient of f(p_j) is 
        //   C * (Cot[\alpha_ij] + Cot[\beta_ij], 
        // where \alpha and \beta are angles facing the edge ij and C is equal
        // to 1.5 / A(p_i), the sum of triangle areas around p_i

        // Compute the weight associated with this neighbor vertex
        double cota = xTriangleGeom[it.TriangleAhead()].xCotangent[
          it.OppositeVertexIndexInTriangleAhead()];
        double cotb = xTriangleGeom[it.TriangleBehind()].xCotangent[
          it.OppositeVertexIndexInTriangleBehind()];
        double weight = scale * (cota + cotb);

        // Set the weight for the moving vertex
        size_t idxMesh = it.GetPositionInMeshSparseArray();
        M.GetSparseData()[xIndexAPhi.xNbrIndex[idxMesh]] = weight;
        M.GetSparseData()[xIndexAOmega.xNbrIndex[idxMesh]] = weight;

        // Update accumulator
        w_accum += weight;
        }

      // Set the diagonal entry in A
      M.GetSparseData()[xIndexAPhi.xSelfIndex[i]] = -w_accum;
      M.GetSparseData()[xIndexAOmega.xSelfIndex[i]] = -w_accum;
      }
    else
      {
      // The matrix will only have the expression for Phi_v, which is the
      // weights of the tangent vectors. Another option would be to put here
      // the expression for Phi_,nu. But let's leave that for now.
      M.GetSparseData()[xIndexN.xSelfIndex[i]] = xLoopScheme.GetOwnWeight(0, i);
      for(EdgeWalkAroundVertex it(topology, i); !it.IsAtEnd(); ++it)
        {
        size_t idxMesh = it.GetPositionInMeshSparseArray();
        M.GetSparseData()[xIndexN.xNbrIndex[idxMesh]] = 
          xLoopScheme.GetNeighborWeight(0, it);
        }
      }
    }
}

void
MeshMedialPDESolver
::FillRHS()
{
  size_t n = topology->nVertices;

  // Fill the right hand side. It's a vector [0 g1 rho g2], so it should be easy
  // to compute it
  for(size_t i = 0; i < n; i++)
    {
    if(topology->IsVertexInternal(i))
      {
      xRHS[i] = 0.0;
      xRHS[n + i] = xAtoms[i].xLapR;
      }
    else
      {
#if defined(TRANSFER_R2)

      // The top part (Dirichlet condition) is trivial
      xRHS[i] = xAtoms[i].R * xAtoms[i].R;

      // This is the hairy part, where we compute
      //   ( - fv * g12 - sqrt(4 * f * g11 - fv * fv * gInv)) / g11
      // Equivalently, using covariant tensor 
      //   ( fv * g12 - sqrt( (4 * f * g22 - fv * fv ) * g ) / g22
      double f = xAtoms[i].F;
      double fv = xAtoms[i].Fv;
      double g12 = xAtoms[i].G.xCovariantTensor[0][1];
      double g22 = xAtoms[i].G.xCovariantTensor[1][1];
      double g = xAtoms[i].G.g;

      // We combine this code with computing the terms for the derivative
      // of the boundary condition. This may have a little overhead, but 
      // it should not be serious. It's nicer to have this code here because
      // all the intermediate variables can be defined only once
      double S2 = (4 * f * g22 - fv * fv );
      double S1 = sqrt( S2 * g );
      double g22inv = 1.0 / g22;
      double S1inv = 1.0 / S1;

      // This is the Neumann condition
      double fu = ( fv * g12 - S1 ) * g22inv;

      // Here are the derivative terms
      Wfu[i].wgt_g22 = (-fu - 2 * f * g * S1inv) * g22inv;
      Wfu[i].wgt_fv  = (g12 + fv * S1inv * g) * g22inv;
      Wfu[i].wgt_g12 = (fv * g22inv);
      Wfu[i].wgt_g = - 0.5 * S2 * S1inv * g22inv;
      Wfu[i].wgt_f = - 2 * S1inv * g;

      // We will save intermediate terms for later derivative computation
      xRHS[n + i] = fu;

#elif defined(TRANSFER_LOGR2)

      // The top part (Dirichlet condition) is trivial
      xRHS[i] = log(xAtoms[i].R * xAtoms[i].R);

      // This is the hairy part, where we compute
      //   ( - fv * g12 - sqrt(4 * f * g11 - fv * fv * gInv)) / g11
      // Equivalently, using covariant tensor 
      //   ( fv * g12 - sqrt( (4 * f * g22 - fv * fv ) * g ) / g22
      double f = xAtoms[i].F;
      double fv = xAtoms[i].Fv;
      double g12 = xAtoms[i].G.xCovariantTensor[0][1];
      double g22 = xAtoms[i].G.xCovariantTensor[1][1];
      double g = xAtoms[i].G.g;

      // We combine this code with computing the terms for the derivative
      // of the boundary condition. This may have a little overhead, but 
      // it should not be serious. It's nicer to have this code here because
      // all the intermediate variables can be defined only once
      double ef = exp(-f);
      double S2 = (4 * ef * g22 - fv * fv );
      double S1 = sqrt( S2 * g );
      double g22inv = 1.0 / g22;
      double S1inv = 1.0 / S1;

      // This is the Neumann condition
      double fu = ( fv * g12 - S1 ) * g22inv;

      // Here are the derivative terms
      Wfu[i].wgt_g22 = (-fu - 2 * f * g * S1inv) * g22inv;
      Wfu[i].wgt_fv  = (g12 + fv * S1inv * g) * g22inv;
      Wfu[i].wgt_g12 = (fv * g22inv);
      Wfu[i].wgt_g = - 0.5 * S2 * S1inv * g22inv;
      Wfu[i].wgt_f = 2 * ef * S1inv * g;

      // We will save intermediate terms for later derivative computation
      xRHS[n + i] = fu;

#endif

      }
    }
}

/*
vnl_vector<double>
MeshMedialPDESolver
::ComputeLBO(const double *phi)
{
  // Compute the RHS
  ComputeMeshGeometry(false);
  FillNewtonMatrix(phi, true);

  // Return the LBO
  vnl_vector<double> lbo(topology->nVertices, 0.0);

  // Compute the lbo at internal vertices
  for(size_t i = 0; i < topology->nVertices; i++)
    {
    if(topology->IsVertexInternal(i))
      {
      // To compute the laplacian of phi, simply use the weights in the corresponding
      // row of sparse matrix A.
      for(size_t j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        {
        lbo[i] += A.GetSparseData()[j] * phi[A.GetColIndex()[j]];
        }
      }
    }

  return lbo;
}
*/

void
MeshMedialPDESolver
::ComputeMedialAtoms(const double *soln)
{
  // Loop over all vertices in the mesh
  for(size_t i = 0; i < topology->nVertices; i++)
    {
    // Get the medial atom corresponding to a
    MedialAtom &a = xAtoms[i];

#if defined(TRANSFER_R2)

    // Set the phi in the atom
    a.F = soln[i];

    // Compute the partials of phi in the tangent directions
    a.Fu = xLoopScheme.GetPartialDerivative(0, i, soln);
    a.Fv = xLoopScheme.GetPartialDerivative(1, i, soln);

    // Compute the boundary atoms from this geometric information
    a.ComputeBoundaryAtoms(!topology->IsVertexInternal(i));


#elif defined(TRANSFER_LOGR2)

    // Set the phi in the atom
    a.F = soln[i];

    // Compute the partials of phi in the tangent directions
    a.Fu = xLoopScheme.GetPartialDerivative(0, i, soln);
    a.Fv = xLoopScheme.GetPartialDerivative(1, i, soln);

    // Compute the partials in R
    a.R = exp(a.F / 2);
    a.Ru = 0.5 * a.R * a.Fu;
    a.Fv = 0.5 * a.R * a.Fv;

    // Compute the boundary atoms from this geometric information
    a.ComputeBoundaryAtomsUsingR(!topology->IsVertexInternal(i));

#endif 

    // Set the 'orig' values in atoms to equal their computed values
    // (the 'orig' values are from BruteForceSubdivisionMedialModel)
    a.xGradRMagSqrOrig = a.xGradRMagSqr;
    a.xNormalFactorOrig = a.xNormalFactor;

    // Check that gradR is reasonable
    if(!a.flagCrest && a.xGradRMagSqr > 0.99)
      {
      // Issue a verbal warning (because this can lead to failed optimization)
      cerr << "WARNING: |gradR| = " << a.xGradRMagSqr << " (too close to 1) in non-crest atom " << i << endl;
      if(a.xGradRMagSqr > 1.0)
        throw MedialModelException("Exception: |gradR| > 1 in non-crest atom");
      }

    else if(a.flagCrest && a.Rs2 > 0.99)
      {
      // Issue a verbal warning (because this can lead to failed optimization)
      cerr << "WARNING: |dR/dS|^2 = " << a.Rs2 << " (too close to 1) in crest atom " << i << endl;
      if(a.Rs2 >= 1.0)
        throw MedialModelException("Exception: |dR/dS| >= 1 in crest atom");
      }
    }
}


void
MeshMedialPDESolver
::SolveEquation(bool flagGradient, bool flagAllowErrors)
{
  // cout << "MeshMedialPDESolver::SolveEquation()" << endl;

  // Compute the mesh geometry
  ComputeMeshGeometry(flagGradient);

  // Compute the sparse matrix
  FillSparseMatrix(true);

  // Compute the right hand side
  FillRHS();

  // Check validity (TODO: remove later)
  for(size_t i = 0; i < xRHS.size(); i++)
    if(std::isinf(xRHS[i]) || std::isnan(xRHS[i]))
      throw MedialModelException("PDE r.h.s. infinite or nan");

  // Check validity (TODO: remove later)
  for(size_t i = 0; i < M.GetNumberOfSparseValues(); i++)
    if(!std::isfinite(M.GetSparseData()[i]))
      throw MedialModelException("PDE sparse matrix infinite or nan");

  // Use pardiso to solve the problem
  xSolver->SetVerbose(false);
  xSolver->SymbolicFactorization(M);
  xSolver->NumericFactorization(M);
  xSolver->Solve(xRHS.data_block(), xSolution.data_block());

  // Check the accuracy of the solution 
  if(!flagAllowErrors)
    {
    static const double MAX_RESIDUAL = 1e-7;
    double residual = (M.MultiplyByVector(xSolution) - xRHS).inf_norm();
    if(residual > MAX_RESIDUAL)
      {
      cerr << "Excessive residual from PDE solver: max(|A*x-b|) = " << residual << endl;
      throw MedialModelException("PDE solver did not solve linear system");
      }
    }

  // The solution must be positive, otherwise this is an ill-posed problem
  if(!flagAllowErrors)
    {
    for(size_t i = 0; i < topology->nVertices; i++)
      {
      if(!std::isfinite(xSolution[i]))
        {
        // Dump the matrix for examination
        ofstream ofs("sparsematdump.txt");
        M.PrintSelfMathematica(ofs);
        ofs << xRHS << endl;
        ofs.close();

        throw MedialModelException("PDE solution infinite or nan");
        }
      else if(xSolution[i] < 0.)
        {
#if defined(TRANSFER_R2)

        throw MedialModelException("PDE solution negative");

#elif defined(TRANSFER_LOGR2)

#endif
        }
      }
    }

  // Compute the medial atoms
  ComputeMedialAtoms(xSolution.data_block());
}



void
MeshMedialPDESolver
::ComputeRHSGradientMatrix()
{
  // Fill the elements of matrix W
  size_t n = topology->nVertices;

  // At this point, the structure of the matrix W has been specified. We have
  // to specify the values. This is done one vertex at a time.
  for(size_t i = 0; i < n; i++)
    {
    // Split over vertex type
    if(topology->IsVertexInternal(i))
      {
      // Compute the laplacian of phi (equal to omega)
      double xLapPhi = xSolution[i + n];
      double xLapOmega = xAtoms[i].xLapR;

      // Scaling factor for cotangent derivatives
      double xScale = 1.5 / xVertexGeom[i].xFanArea;

      // Reference to the weight for the fixed vertex
      SMLVec3d &wiiPhi = WX[xIndexAPhi.xSelfIndex[i]];
      SMLVec3d &wiiOmega = WX[xIndexAOmega.xSelfIndex[i]];
      wiiPhi.fill(0.0);
      wiiOmega.fill(0.0);

      // Get the value of phi at the center
      double phiFixed = xSolution[i];
      double omegaFixed = xSolution[i+n];

      // Go over the neighbors
      for(EdgeWalkAroundVertex it(topology, i) ; !it.IsAtEnd(); ++it)
        {
        // The weight vector for the moving vertex
        size_t idxMesh = it.GetPositionInMeshSparseArray();
        size_t idxAPhi = xIndexAPhi.xNbrIndex[idxMesh];
        size_t idxAOmega = xIndexAOmega.xNbrIndex[idxMesh];

        SMLVec3d &wijPhi = WX[idxAPhi];
        SMLVec3d &wijOmega = WX[idxAOmega];

        // Get the triangle index in front and behind
        size_t ta = it.TriangleAhead(), tb = it.TriangleBehind();
        size_t va = it.MovingVertexIndexInTriangleAhead();
        size_t vb = it.MovingVertexIndexInTriangleBehind();
        size_t qa = it.OppositeVertexIndexInTriangleAhead();
        size_t qb = it.OppositeVertexIndexInTriangleBehind();

        // Get the phi's at the surrounding vertices
        double phiAhead = xSolution[it.VertexIdAhead()];
        double phiBehind = xSolution[it.VertexIdBehind()];
        double phiMoving = xSolution[it.MovingVertexId()];

        double omegaAhead = xSolution[n + it.VertexIdAhead()];
        double omegaBehind = xSolution[n + it.VertexIdBehind()];
        double omegaMoving = xSolution[n + it.MovingVertexId()];

        SMLVec3d wMF = xTriangleGeom[ta].xCotGrad[qa][va] + xTriangleGeom[tb].xCotGrad[qb][vb];
        SMLVec3d wAF = xTriangleGeom[ta].xCotGrad[va][va];
        SMLVec3d wBF = xTriangleGeom[tb].xCotGrad[vb][vb];
        SMLVec3d dA = 
          (xTriangleGeom[ta].xAreaGrad[va] + xTriangleGeom[tb].xAreaGrad[vb])
          / xVertexGeom[i].xFanArea;

        // Compute the weights
        wijPhi = xScale * (
          (phiMoving - phiFixed) * wMF + 
          (phiAhead  - phiFixed) * wAF + 
          (phiBehind - phiFixed) * wBF )
          - xLapPhi * dA;
        wiiPhi -= wijPhi;

        // Compute the weights
        wijOmega = xScale * (
          (omegaMoving - omegaFixed) * wMF + 
          (omegaAhead  - omegaFixed) * wAF + 
          (omegaBehind - omegaFixed) * wBF )
          - xLapOmega * dA;
        wiiOmega -= wijOmega;
        }
      }
    }
}

void
MeshMedialPDESolver
::BeginGradientComputation()
{
  size_t i;

  // Compute the weight matrix for the right hand side computations
  this->ComputeRHSGradientMatrix();

  // Precompute common terms for atom derivatives  
  for(i = 0; i < topology->nVertices; i++)
    xAtoms[i].ComputeCommonDerivativeTerms(xTempDerivativeTerms[i]);
}

void
MeshMedialPDESolver
::ComputeAtomVariationalDerivative(MedialAtom *dAtoms)
{
  size_t i;

  const TriangleMesh::NeighborMatrix &NB = topology->GetNeighborMatrix();
  size_t n = topology->nVertices;

  // Each variational derivative is a linear system Ax = b. The matrix A
  // should be set correctly from the last call to Solve() and all we have to
  // change are the vectors B.
  vnl_vector<double> rhs(2 * n, 0.0);
  vnl_vector<double> soln(2 * n, 0.0);

  // Compute F at boundary atoms, reset at edge atoms
  for(i = 0; i < n; i++) 
    {

#if defined(TRANSFER_R2)

    dAtoms[i].F = topology->IsVertexInternal(i) ? 
      0.0 : 2.0 * xAtoms[i].R * dAtoms[i].R;

#elif defined(TRANSFER_LOGR2)
    
    dAtoms[i].F = topology->IsVertexInternal(i) ? 
      0.0 : 2.0 * dAtoms[i].R / xAtoms[i].R;

#endif
    }

  // Repeat for each atom
  for(i = 0; i < n; i++) 
    {

    // First, we need some derivatives
    // TODO: All of these can be precomputed in SetVariationalBasis
    MedialAtom &da = dAtoms[i];
    da.Xu = xLoopScheme.Xu(i, dAtoms);
    da.Xv = xLoopScheme.Xv(i, dAtoms);
    da.Xuu = da.Xuv = da.Xvv = SMLVec3d(0.0); 
    da.Fv = xLoopScheme.Fv(i, dAtoms);

    // Compute the metric tensor derivatives of the atom
    xAtoms[i].ComputeMetricTensorDerivatives(dAtoms[i]);

    // Compute the right hand side expression
    if(topology->IsVertexInternal(i))
      {
      // For internal vertices, the derivative wrt X is encoded in WX
      rhs[i]   = dot_product(dAtoms[i].X, WX[xIndexAPhi.xSelfIndex[i]]);
      SMLVec3d wx = WX[xIndexAOmega.xSelfIndex[i]];
      rhs[i+n] = dot_product(dAtoms[i].X, WX[xIndexAOmega.xSelfIndex[i]]);

      // Go over the neighbors of i
      for(size_t k = NB.GetRowIndex()[i]; k < NB.GetRowIndex()[i+1]; k++)
        {
        size_t j = NB.GetColIndex()[k];
        SMLVec3d &dX = dAtoms[j].X;
        rhs[i]   += dot_product(dX, WX[xIndexAPhi.xNbrIndex[k]]);
        rhs[i+n] += dot_product(dX, WX[xIndexAOmega.xNbrIndex[k]]);        
        }

      // There is also a derivative wrt rho
      rhs[i+n] -= dAtoms[i].xLapR;
      }

    // For edge vertices things are a bit more complicated
    else
      {
      // The top part is simple, of course (derivative of Tau)
      rhs[i] = - dAtoms[i].F;
      
      // The bottom part is the derivative of that hideous square root
      rhs[i+n] = - (
        Wfu[i].wgt_g12 * dAtoms[i].G.xCovariantTensor[0][1] + 
        Wfu[i].wgt_g22 * dAtoms[i].G.xCovariantTensor[1][1] +
        Wfu[i].wgt_g   * dAtoms[i].G.g +
        Wfu[i].wgt_f   * dAtoms[i].F +
        Wfu[i].wgt_fv  * dAtoms[i].Fv);
      }
    }

  // Solve the partial differential equation (dPhi/dVar)
  xSolver->Solve(rhs.data_block(), soln.data_block());

  

  /* // ( a little test code )

  // At this point, test the gradient computation
  double eps = 1.0e-4;
  MedialAtom *m1 = new MedialAtom[n];
  MedialAtom *m2 = new MedialAtom[n];
  for(size_t i = 0; i < n; i++)
    {
    m1[i].X     = xAtoms[i].X     - eps * dAtoms[i].X;
    m1[i].R     = xAtoms[i].R     - eps * dAtoms[i].R;
    m1[i].xLapR = xAtoms[i].xLapR - eps * dAtoms[i].xLapR;
    m1[i].F = m1[i].R * m1[i].R;

    m2[i].X     = xAtoms[i].X     + eps * dAtoms[i].X;
    m2[i].R     = xAtoms[i].R     + eps * dAtoms[i].R;
    m2[i].xLapR = xAtoms[i].xLapR + eps * dAtoms[i].xLapR;
    m2[i].F = m2[i].R * m2[i].R;
    }

  MeshMedialPDESolver s1;
  s1.SetMeshTopology(topology, m1);
  s1.SolveEquation();

  MeshMedialPDESolver s2;
  s2.SetMeshTopology(topology, m2);
  s2.SolveEquation();

  Vec q1 = s1.M.MultiplyByVector(xSolution) - s1.xRHS;
  Vec q2 = s2.M.MultiplyByVector(xSolution) - s2.xRHS;
  Vec dqApp = (q2-q1) / (2 * eps);

  cout << "Error: " << (dqApp - rhs).inf_norm() << endl;

  Vec p1 = s1.xSolution;
  Vec p2 = s2.xSolution;
  Vec dpApp = (p2-p1) / (2 * eps);

  */ 

  // Compute each atom
  for(i = 0; i < n; i++) 
    {
    MedialAtom &da = dAtoms[i];
    
#if defined(TRANSFER_R2)

    // For each atom, compute F, Fu, Fv
    da.F = - soln[i];
    da.Fu = - xLoopScheme.GetPartialDerivative(0, i, soln.data_block());
    da.Fv = - xLoopScheme.GetPartialDerivative(1, i, soln.data_block());

    xAtoms[i].ComputeBoundaryAtomDerivatives(dAtoms[i], xTempDerivativeTerms[i]);

#elif defined(TRANSFER_LOGR2)

    da.F = - soln[i];
    da.Fu = - xLoopScheme.GetPartialDerivative(0, i, soln.data_block());
    da.Fv = - xLoopScheme.GetPartialDerivative(1, i, soln.data_block());

    da.R = 0.5 * xAtoms[i].R * da.F;
    da.Ru = 0.5 * xAtoms[i].R * da.Fu;
    da.Rv = 0.5 * xAtoms[i].R * da.Fv;

    // Compute the derivatives of the boundary nodes
    xAtoms[i].ComputeBoundaryAtomDerivativesUsingR(dAtoms[i], xTempDerivativeTerms[i]);

#endif

    da.xGradRMagSqrOrig = da.xGradRMagSqr;
    da.xNormalFactorOrig = da.xNormalFactor;
    }


}





#ifdef DOITLATER

int
MeshMedialPDESolver
::TestPartialDerivatives()
{
  bool flagSuccess = true;
  size_t i, j, k;
  
  // Epsilon for central differences
  double xEpsilon = 0.00001;
  double z = 0.5 / xEpsilon;

  // Create a pair of atom arrays
  MedialAtom *atom1 = new MedialAtom[topology->nVertices];
  MedialAtom *atom2 = new MedialAtom[topology->nVertices];
  
  // Create a pair of alternative solutions
  MeshMedialPDESolver m1, m2;
  m1.SetMeshTopology(topology, atom1);
  m2.SetMeshTopology(topology, atom2);

  // Compute the solution in the current model
  this->SolveEquation(NULL, true);

  // Compute the gradient terms
  this->ComputeRHSGradientMatrix();

  // Random generator
  vnl_random rnd;

  // Create a pair of offsets (random variations)
  vector<SMLVec3d> varX(topology->nVertices);
  vector<double> varRho(topology->nVertices);

  // Create the solution vector (current solution)
  vnl_vector<double> phi(topology->nVertices);

  // Populate the variations
  for(i = 0; i < topology->nVertices; i++)
    {
    // Set the offsets
    varX[i][0] = rnd.drand32(-1.0, 1.0);
    varX[i][1] = rnd.drand32(-1.0, 1.0);
    varX[i][2] = rnd.drand32(-1.0, 1.0);
    varRho[i] = rnd.drand32(-1.0, 1.0);
    
    // Set phi
    phi[i] = xAtoms[i].F;

    // Apply to the atoms
    m1.xAtoms[i].X = xAtoms[i].X + xEpsilon * varX[i];
    m1.xAtoms[i].xLapR = xAtoms[i].xLapR + xEpsilon * varRho[i];
    m1.xAtoms[i].F = phi[i];

    m2.xAtoms[i].X = xAtoms[i].X - xEpsilon * varX[i];
    m2.xAtoms[i].xLapR = xAtoms[i].xLapR - xEpsilon * varRho[i];
    m2.xAtoms[i].F = phi[i];
    }

  // Compute the geometry in all three cases
  m1.ComputeMeshGeometry(false);
  m2.ComputeMeshGeometry(false);

  // Compare the gradient with the results
  for(j = 0; j < topology->triangles.size(); j++)
    {
    // Get the triangles
    TriangleGeom &T = xTriangleGeom[j];
    TriangleGeom &T1 = m1.xTriangleGeom[j];
    TriangleGeom &T2 = m2.xTriangleGeom[j];

    // Look at the triangle areas
    double cdArea = z * (T1.xArea - T2.xArea);
    double adArea = 
      dot_product(T.xAreaGrad[0], varX[topology->triangles[j].vertices[0]]) +
      dot_product(T.xAreaGrad[1], varX[topology->triangles[j].vertices[1]]) +
      dot_product(T.xAreaGrad[2], varX[topology->triangles[j].vertices[2]]);
    if(fabs(cdArea - adArea) > 2.0 * xEpsilon)
      {
      flagSuccess = false;
      cout << "j = " << std::setw(5) << j 
        << "; cdArea = " << std::setw(12) << cdArea 
        << "; adArea = " << std::setw(12) << adArea 
        << "; delta = " << std::setw(12) << fabs(cdArea - adArea) << endl;
      }

    // Look at the cotangents of the angles
    for(k = 0; k < 3; k++)
      {
      double cdCot = z * (T1.xCotangent[k] - T2.xCotangent[k]);
      double adCot =
        dot_product(T.xCotGrad[k][0], varX[topology->triangles[j].vertices[0]]) +
        dot_product(T.xCotGrad[k][1], varX[topology->triangles[j].vertices[1]]) +
        dot_product(T.xCotGrad[k][2], varX[topology->triangles[j].vertices[2]]);

      // We have to allow for some slack with the Cot derivatives because for
      // small angles the third derivative can be quite large and the central
      // difference approximation can be inaccurate
      if(fabs(cdCot - adCot) > 2.0 * xEpsilon * std::max(1.0, fabs(T.xCotangent[k])))
        {
        flagSuccess = false;
        cout << "j = " << std::setw(5) << j 
          << "; cdCot = " << std::setw(12) << cdCot 
          << "; adCot = " << std::setw(12) << adCot 
          << "; delta = " << std::setw(12) << fabs(cdCot - adCot) << endl;
        }
      }
    }

  // Access matrix A in the offset models
  size_t *xRowIndex = A.GetRowIndex();
  size_t *xColIndex = A.GetColIndex();
  double *A1 = m1.A.GetSparseData();
  double *A2 = m2.A.GetSparseData();

  // Perform the gradient computation
  m1.FillNewtonMatrix(phi.data_block(), true);
  m2.FillNewtonMatrix(phi.data_block(), true);

  // Now, compute the derivative of the laplacian at every vertex
  for(i = 0; i < topology->nVertices; i++)
    {
    if(topology->IsVertexInternal(i))
      {
      double L1 = 0, L2 = 0, cdLap = 0, adLap = 0;
      for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++)
        {
        L1 += A1[j] * phi[xColIndex[j]];
        L2 += A2[j] * phi[xColIndex[j]];
        adLap += dot_product(W[j], varX[xColIndex[j]]);
        }
      cdLap = z * (L1 - L2);
      
      if(fabs(cdLap - adLap) > 2.0 * xEpsilon)
        {
        flagSuccess = false;
        cout << "i = " << std::setw(5) << i 
          << "; cdLap = " << std::setw(12) << cdLap 
          << "; adLap = " << std::setw(12) << adLap 
          << "; delta = " << std::setw(12) << fabs(cdLap - adLap) << endl;
        }
      }
    else
      {
      double GM1 = 0, GM2 = 0, cdGM = 0, adGM = 0;
      double Fu = 0, Fv = 0;
      for(j = xRowIndex[i]; j < xRowIndex[i+1]; j++) 
        {
        adGM += dot_product(W[j], varX[xColIndex[j]]);
        }

      // Compute gradient magnitude of phi on both subdivision surfaces
      GM1 = m1.xAtoms[i].G.SquaredGradientMagnitude(xAtoms[i].Fu, xAtoms[i].Fv);
      GM2 = m2.xAtoms[i].G.SquaredGradientMagnitude(xAtoms[i].Fu, xAtoms[i].Fv);
      cdGM = z * (GM1 - GM2);

      if(fabs(cdGM - adGM) > 2.0 * xEpsilon)
        {
        flagSuccess = false;
        cout << "i = " << std::setw(5) << i 
          << "; cdGM = " << std::setw(12) << cdGM 
          << "; adGM = " << std::setw(12) << adGM 
          << "; delta = " << std::setw(12) << fabs(cdGM - adGM) << endl;
        }
      }
    }

  // Clean up
  delete atom1; delete atom2;

  return flagSuccess ? 0 : -1;
}

#endif

/**************************************************************************
 * CODE FOR SOLVING THE PDE USING TENSOLVE METHOD
 *************************************************************************/
#ifdef GARBAGE

extern "C" 
{
  void tsneci_(
    int &maxm,int &maxn,int &maxp,double *x0,int &m,int &n,double *typx,double *typf,int &itnlim,
    int &jacflg, double &gradtl, double &steptl, double &ftol,int &method,int &global,
     double &stepmx, double &dlt,int &ipr,double *wrkunc,int &lunc,double *wrknem,int &lnem,
    double *wrknen,int &lnen,int *wrkn,int &lin,
    void (*tensolve_mpde_residual)(double *, double *, int &, int &),
    void (*tensolve_mpde_jacobian)(double *, double *, int &, int &, int &),
    int &msg, double *xp,double *fp,double *gp,int &termcd);

};

// A global variable pointing to the problem to be solved by the tensolve
// algorithm
MeshMedialPDESolver *tensolve_mpde_problem = NULL;

// Method that computes the residual of the tensolve problem

void tensolve_mpde_residual(
  double *x,
  double *f,
  int &m,
  int &n)
{
  cout << "." << flush;
  // Get a handle to the current problem
  MeshMedialPDESolver *mpde = tensolve_mpde_problem;

  // Evaluate the solution at point x (x == phi)
  double F = mpde->FillNewtonRHS(x);

  // Copy the values into the vector f
  for(size_t i = 0; i < m; i++)
    f[i] = -mpde->xRHS[i];
}

// Method that computes the Jacobian of the tensolve problem
void tensolve_mpde_jacobian(
  double *x,
  double *J,
  int &maxm,
  int &m,
  int &n)
{
  cout << "JAC" << endl;
  // Get a handle to the current problem
  MeshMedialPDESolver *mpde = tensolve_mpde_problem;

  // Evaluate the solution at point x (x == phi)
  mpde->FillNewtonRHS(x);
  mpde->FillNewtonMatrix(x, false);

  // Initialize the Jacobian to zeros
  for(size_t q = 0; q < m * n; q++)
    J[q] = 0.0;

  // Copy the values into the output Jacobian matrix
  MeshMedialPDESolver::SparseMat &A = mpde->A;
  for(size_t i = 0; i < A.GetNumberOfRows(); i++)
    {
    for(size_t j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
      {
      // Get the column number and the Jacobian value
      size_t k = A.GetColIndex()[j];
      double jval = A.GetSparseData()[j];

      // Compute the Fortran index. Fortran arrays are column major. Than
      // means, the index is the (column number) * (# rows) + (row number)
      // The Jacobian is M x N (M is the number of equations, N unknowns),
      // so the entry j must be inserted into row i, column k
      size_t idx = k * maxm + i;

      /*
      vnl_vector<double> xx(x, n);
      vnl_vector<double> x1 = xx, x2 = xx;
      double eps = 1.0e-6;
      x1[k] += eps; x2[k] -= eps;
      mpde->FillNewtonRHS(x1.data_block());
      double f1 = -mpde->xRHS[i];
      mpde->FillNewtonRHS(x2.data_block());
      double f2 = -mpde->xRHS[i];

      double jfd = (f1 - f2) / (2 * eps);

      
      if(abs(jfd - jval) > eps)
        printf("Jac[i][k] = %g\t, FDJ = %g\n", jval, jfd);
      */

      
      

      // Set the Jacobian value
      J[idx] = jval;
      }
    }
}

// Method that executes the tensolve algorithm to solve the MPDE
void tensolve_mpde_solve(MeshMedialPDESolver *solver, double *x)
{
  // Store the problem
  tensolve_mpde_problem = solver;

  // Initialize all the parameters of the tensolve call
  int m = solver->A.GetNumberOfRows();
  int n = solver->A.GetNumberOfColumns();
  int maxm = m + n + 2;
  int maxn = n + 2;
  int maxp = (int)ceil(sqrt(n));
  double *x0 = new double[n];
  double *typx  = new double[n];
  double *typf = new double[m];
  int itnlim = 0;
  int jacflg = 1;
  double gradtl = -1;
  double steptl = -1;
  double ftol = -1;
  int method = 1;
  int global = 1;
  double stepmx = -1;
  double dlt = -1.0;
  int ipr = 6;
  int lunc = (int)(2 * ceil(sqrt(n)) + 4);
  double *wrkunc = new double[lunc * maxp];
  int lnem = (int)(n + 2 * ceil(sqrt(n))+11);
  double *wrknem = new double[lnem * maxm];
  int lnen = (int)(2 * ceil(sqrt(n)) + 9);
  double *wrknen = new double[lnen * maxn];
  int lin = 3;
  int *iwrkn = new int[lin * maxn];
  int msg = 0;
  double *xp = new double[n];
  double *fp = new double[m];
  double *gp = new double[n];
  int termcd = 0;

  // Initialize the x0 vector and the scales vector
  for(size_t i = 0; i < n; i++)
    {
    x0[i] = x[i];
    typx[i] = 1.0;
    }

  for(size_t j = 0; j < m; j++)
    {
    typf[j] = 1.0;
    }

  // Call the TENSOLVE routine
  tsneci_(maxm,maxn,maxp,x0,m,n,typx,typf,itnlim,
     jacflg,gradtl,steptl,ftol,method,global,
     stepmx,dlt,ipr,wrkunc,lunc,wrknem,lnem,
     wrknen,lnen,iwrkn,lin,
     &tensolve_mpde_residual,
     &tensolve_mpde_jacobian,
     msg, xp,fp,gp,termcd);

  // Clear memory
  delete wrkunc;
  delete wrknem;
  delete wrknen;
  delete iwrkn;
  delete xp;
  delete fp;
  delete gp;
  delete x0;
  delete typx;
  delete typf;
}

#endif // GARBAGE
