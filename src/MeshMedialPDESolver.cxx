#include "MeshMedialPDESolver.h"
#include "MedialAtomGrid.h"
#include <iomanip>
#include <vector>
#include <vnl/vnl_math.h>
#include <vnl/vnl_random.h>
#include <algorithm>

using namespace std;

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
  SymmetricPositiveDefiniteRealPARDISO pardiso;
  pardiso.SymbolicFactorization(
    n, xPardisoRowIndex.data_block(), xPardisoColIndex.data_block(), A.GetSparseData());

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
  while(!found && k++ < n_iter)
    {
    // Let A = J' * J + mu I
    for(i = 0; i < n; i++)
      A.GetSparseData()[A.GetRowIndex()[i]] += mu;
   
    // Solve the linear system (J' * J + mu I) x = -g
    pardiso.NumericFactorization(A.GetSparseData());
    pardiso.Solve(minus_g.data_block(), h.data_block());

    // Compute the norm of x and h
    double norm_h = h.magnitude();
    double norm_x = x.magnitude();

    // Report the solution
    bool verbose = true;
    if(verbose)
      printf(
        "   -> LM Iter: %3d    F = %8.4e    fmax = %8.4e    fmean = %8.4e"
        "    mu = %8.4e    "
        "|h| = %8.4e    |x| = %8.4e\n",
        k, normsq_fx / 2, fx.inf_norm(), 
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
    { delete x; x = NULL; }
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

    // Now create entries for the sparse matrix
    for(vector<int>::iterator ir = ring.begin(); ir != ring.end(); ++ir)
      {
      // Get the vertex id and index in the mesh's neighbor sparse array
      size_t j = *ir;

      if(topology->IsVertexInternal(i))
        {
        // Add entry in A(phi)
        S[i].push_back(make_pair(j, 0.0));

        // Add entry in A(omega)
        S[n+i].push_back(make_pair(n+j, 0.0));

        // Add entry in U1
        if(i == j)
          S[i].push_back(make_pair(n + i, -1.0));
        }
      else
        {
        // Add entry in N
        S[i+n].push_back(make_pair(j, 0.0));

        // Add entry in U2
        if(i == j)
          S[i].push_back(make_pair(i, 1.0));
        }      
      }
    }

  // Now create the sparse matrix
  M.SetFromSTL(S, 2 * n);

  // Perform PARDISO symbolic factorization
  xPardiso.SymbolicFactorization(M);

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
      else if(i < n)
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
          xIndexAOmega.xNbrIndex[Nt.FindEntryIndex(i, j)] = rM.SparseIndex();
        }
      }
    }

  // Initialize right hand side and solution arrays
  xRHS.set_size(2 * n);
  xPhi.set_size(2 * n);

  // Initialize the triangle area array and other such arrays
  xTriangleGeom = new TriangleGeom[topology->triangles.size()];
  xVertexGeom = new VertexGeom[topology->nVertices];

  // Initialize the tangent weights
  xLoopScheme.SetMesh(topology);

  // Copy the pointer to the atoms. The atoms are passed in to this method
  // and are managed by the parent
  xAtoms = inputAtoms;

  // Initialize the weight array for gradient computation
  W = new SMLVec3d[M.GetNumberOfSparseValues()];
}

MeshMedialPDESolver
::MeshMedialPDESolver()
{
  xTriangleGeom = NULL;
  xVertexGeom = NULL;
  xMapVertexNbrToA = NULL;
  xMapVertexToA = NULL;
  topology = NULL;
  xPardisoColIndex = NULL;
  xPardisoRowIndex = NULL;
  xAtoms = NULL;
  W = NULL;
}

MeshMedialPDESolver
::~MeshMedialPDESolver()
{
  Reset();
}

void
MeshMedialPDESolver
::Reset()
{
  reset_ptr(xTriangleGeom);
  reset_ptr(xVertexGeom);
  reset_ptr(xMapVertexNbrToA);
  reset_ptr(xMapVertexToA);
  reset_ptr(xPardisoColIndex); 
  reset_ptr(xPardisoRowIndex);
  reset_ptr(W);
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

    // Compute the normal vector
    a.ComputeNormalVector();

    // For the boundary atoms, also compute Fu, which is static
    if(!topology->IsVertexInternal(i))
      a.Fu = xLoopScheme.Fu(i, xAtoms);
    } 
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
      M.GetSparseData()[xIndexN.xSelfIndex[i]] = xLoopScheme.GetOwnWeight(1, i);
      for(EdgeWalkAroundVertex it(topology, i); !it.IsAtEnd(); ++it)
        {
        size_t idxMesh = it.GetPositionInMeshSparseArray();
        M.GetSparseData()[xIndexN.xNbrIndex[idxMesh]] = 
          xLoopScheme.GetNeighborWeight(1, it);
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
      // This is the hairy part, where we compute
      // (-fu g12 - Sqrt[fu^2 g12^2 + 4 f0 g22 - fu^2 g11 g22])/g22
      double f = xAtoms[i].F;
      double fu = xAtoms[i].Fu;
      double g11 = xAtoms[i].G.xContravariantTensor[1][1];
      double g12 = xAtoms[i].G.xContravariantTensor[1][2];
      double g22 = xAtoms[i].G.xContravariantTensor[2][2];
      double gInv = xAtoms[i].G.gInv;

      xRHS[i] = xAtoms[i].R * xAtoms[i].R;
      xRHS[n + i] = 
        ( - fu * g12 - sqrt(4 * f * g22 - fu * fu * gInv)) / g22;
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

    // Set the phi in the atom
    a.F = soln[i];

    // Compute the partials of phi in the tangent directions
    a.Fu = xLoopScheme.GetPartialDerivative(0, i, soln);
    a.Fv = xLoopScheme.GetPartialDerivative(1, i, soln);

    // Compute the boundary atoms from this geometric information
    a.ComputeBoundaryAtoms(!topology->IsVertexInternal(i));
    }
}

void
MeshMedialPDESolver
::SetInputData(const SMLVec3d *X, const double *rho, const double *phi)
{
  // Fill the X and rho values of the atoms
  for(size_t i = 0; i < topology->nVertices; ++i)
    {
    xAtoms[i].X = X[i];
    xAtoms[i].xLapR = rho[i];
    if(phi)
      xAtoms[i].F = phi[i];
    }
}

/* 
int
MeshMedialPDESolver
::TestJacobianAndGradient(double *xInitSoln)
{
  size_t i;

  // Success?
  bool flagSuccess = true;

  // Compute the mesh geometry
  ComputeMeshGeometry(false);

  // Set the solution to the current values of phi (or to the initial solution
  // passed in to the method)
  vnl_vector<double> xSoln(topology->nVertices);
  for(i = 0; i < topology->nVertices; i++)
    xSoln[i] = xInitSoln == NULL ? xAtoms[i].F : xInitSoln[i];

  // Compute the Jacobian matrix
  FillNewtonMatrix(xSoln.data_block(), true);

  // Test whether the Jacobian matrix is really Jacobian
  double eps = 0.0001;
  for(i = 0; i < topology->nVertices; i++)
    {
    for(SparseMat::RowIterator r = A.Row(i); !r.IsAtEnd(); ++r)
      {
      size_t j = r.Column(); double xj = xSoln[j];

      // Compute central difference gradient
      xSoln[j] = xj + eps;
      double f1 = ComputeNodeF(i, xSoln.data_block());
      xSoln[j] = xj - eps;
      double f2 = ComputeNodeF(i, xSoln.data_block());
      xSoln[j] = xj;

      double cd = (f1 - f2) * 0.5 / eps;
      double ad = r.Value();
      if(fabs(ad - cd) > eps)
        {
        flagSuccess = false;
        cout << "Jacobian test failed: " << ad << " <> " << cd << endl;
        }
      }
    }

  // Compute finite difference gradient of 0.5 F.F
  Vec cdGrad(topology->nVertices, 0.0);
  for(i = 0; i < topology->nVertices; i++)
    {
    // Compute the central difference
    double xi = xSoln[i];
    xSoln[i] = xi + eps; double f1 = 0.5 * this->FillNewtonRHS(xSoln.data_block());
    xSoln[i] = xi - eps; double f2 = 0.5 * this->FillNewtonRHS(xSoln.data_block());
    xSoln[i] = xi;
    cdGrad[i] = (f1 - f2) * 0.5 / eps;
    }

  // Compute the gradient using Jacobian
  this->FillNewtonRHS(xSoln.data_block());
  Vec adGrad = - A.MultiplyTransposeByVector(xRHS);

  // Report the difference
  double xGradDiff = (adGrad - cdGrad).inf_norm();
  if(xGradDiff > eps)
    {
    flagSuccess = false;
    cout << "Gradient test failed: " << xGradDiff << endl;
    }

  return flagSuccess ? 0 : -1;
}

*/

/*
void
MeshMedialPDESolver
::BacktrackNRC(
  const Vec &xold,          // Starting point for Newton's method
  double fold,              // Value of the function 0.5 F.F at this point
  const Vec &g,             // The gradient direction at xold
  const Vec &p,             // The Newton step direction
  Vec &x,                   // Output: the new point along Newton direction
  double &f)                // Output: the value of the function at the new point
{


}
*/

/*
 NEWTON METHOD


void
MeshMedialPDESolver
::SolveEquation(double *xInitSoln, bool flagGradient)
{
  size_t i;

  cout << "MeshMedialPDESolver::SolveEquation()" << endl;

  // Compute the mesh geometry
  ComputeMeshGeometry(flagGradient);

  // Set the solution to the current values of phi (or to the initial solution
  // passed in to the method)
  vnl_vector<double> xSoln(topology->nVertices);
  for(i = 0; i < topology->nVertices; i++)
    xSoln[i] = xInitSoln == NULL ? xAtoms[i].F : xInitSoln[i];
  vnl_vector<double> xStartingSoln = xSoln;

  // Flag indicating that Newton's method converged (without tricks)
  bool flagConverge = false;

  // Run for up to 50 iterations
  for(i = 0; i < 20; i++)
    {
    // First, fill in the A matrix (only a part of it during most iterations)
    FillNewtonMatrix(xSoln.data_block(), i == 0);

    // Now, fill in the right hand side 
    double xStartingRHS = FillNewtonRHS(xSoln.data_block());

    printf("Iteration: %3d\t Max RHS: %g\t Max PHI: %g\n", 
      i, xRHS.inf_norm(), xSoln.inf_norm());

    // If the right hand side is close enough to zero, we can exit
    if(xRHS.inf_norm() < 1.0e-10) 
      {
      // Set the converge flag
      flagConverge = true;

      // Exit the loop
      break;
      }

    // During the first iteration perform factorization
    if(i == 0)
      xPardiso.SymbolicFactorization(
        topology->nVertices, xPardisoRowIndex, xPardisoColIndex, A.GetSparseData());

    // In the subsequent iterations, only do the numeric factorization and solve
    xPardiso.NumericFactorization(A.GetSparseData());
    xPardiso.Solve(xRHS.data_block(), xEpsilon.data_block());

    // Add epsilon to the current guess
    xSoln += xEpsilon;
    }
  
  // If the system fails to converge, use the SVD-based method, which works
  // with ill-conditioned Jacobians
  if(!flagConverge)
    {
    // Revert to the starting point
    xSoln = xStartingSoln;

    // Do a loop
    for(i = 0; i < 20; i++)
      {
      // Compute the Jacobian matrix and the right hand side
      FillNewtonMatrix(xSoln.data_block(), false);

      // Now, fill in the right hand side 
      FillNewtonRHS(xSoln.data_block());

      // Compute the SVD of the Jacobian matrix

      }



    cout << "{Conv. Fail. " << xRHS.inf_norm() << "}" << endl;

    throw MedialModelException("Convergence Failure");
    }

  // If Newton fails to converge revert to gradient descent
#ifdef GARBAGE  
  if(xPostSolveRHS > 1e-10) 
    {
    cout << "REVERTING TO GRADIENT DESCENT" << endl;

    // Set up the optimizer with target functions
    gsl_multimin_function_fdf my_func;
    my_func.f = &MeshMedialPDESolver::gsl_f;
    my_func.df = &MeshMedialPDESolver::gsl_df;
    my_func.fdf = &MeshMedialPDESolver::gsl_fdf;
    my_func.n = topology->nVertices;
    my_func.params = this;

    // Create the initial solution
    gsl_vector *my_start = gsl_vector_alloc(my_func.n);
    memcpy(my_start->data, xSoln.data_block(), sizeof(double) * my_func.n);

    // Create conjugate gradient minimizer
    gsl_multimin_fdfminimizer *my_min = 
      gsl_multimin_fdfminimizer_alloc( gsl_multimin_fdfminimizer_conjugate_pr, my_func.n);

    // Set up the parameters of the optimizer
    gsl_multimin_fdfminimizer_set(my_min, &my_func, my_start, 0.1, 1e-7);

    // Iterate until we converge
    while(0 == gsl_multimin_fdfminimizer_iterate(my_min))
      {
      cout << gsl_multimin_fdfminimizer_minimum(my_min) << endl;
      }
    cout << endl;

    // Report
    cout << "Conj Grad Final Value: " << gsl_multimin_fdfminimizer_minimum(my_min);
    
    // Get the best solution
    xSoln.copy_in(gsl_multimin_fdfminimizer_x(my_min)->data);
    }
#endif GARBAGE

  // Compute the medial atoms
  ComputeMedialAtoms(xSoln.data_block());
}

*/

void
MeshMedialPDESolver
::SolveEquation(bool flagGradient)
{
  // cout << "MeshMedialPDESolver::SolveEquation()" << endl;

  // Compute the mesh geometry
  ComputeMeshGeometry(flagGradient);

  // Compute the sparse matrix
  FillSparseMatrix(true);

  // Compute the right hand side
  FillRHS();

  // Use pardiso to solve the problem
  xPardiso.NumericFactorization(M);
  xPardiso.Solve(xRHS.data_block(), xPhi.data_block());

  // Compute the medial atoms
  ComputeMedialAtoms(xPhi.data_block());
}

/*

double
MeshMedialPDESolver::gsl_f(const gsl_vector *x, void *params)
{
  // Get the pointer to the solver
  MeshMedialPDESolver *solver = static_cast<MeshMedialPDESolver *>(params);

  // Compute the function for this data
  return solver->FillNewtonRHS(x->data);
}

void
MeshMedialPDESolver::gsl_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df)
{
  // Get the pointer to the solver
  MeshMedialPDESolver *solver = static_cast<MeshMedialPDESolver *>(params);

  // Compute the Jacobian matrix and the function value
  solver->FillNewtonMatrix(v->data, false);
  *f = solver->FillNewtonRHS(v->data);

  // Compute the gradient vector
  Vec grad = -2.0 * solver->A.MultiplyTransposeByVector(solver->xRHS);

  // Copy the gradient in
  std::copy(grad.data_block(), grad.data_block() + grad.size(), df->data);
}

void
MeshMedialPDESolver::gsl_df(const gsl_vector *v, void *params, gsl_vector *df)
{
  double dummy;
  MeshMedialPDESolver::gsl_fdf(v, params, &dummy, df);
}

*/


void
MeshMedialPDESolver
::BeginGradientComputation()
{
}

void
MeshMedialPDESolver
::ComputeAtomVariationalDerivative(MedialAtom *dAtoms)
{
}

#ifdef DOITLATER

void
MeshMedialPDESolver
::ComputeRHSGradientMatrix()
{
  // Fill the elements of matrix W
  size_t i, j;

  // At this point, the structure of the matrix W has been specified. We have
  // to specify the values. This is done one vertex at a time.
  for(i = 0; i < topology->nVertices; i++)
    {
    // Get a reference to the atom, which contains the current solution
    MedialAtom &a = xAtoms[i];

    // Walk around the vertex
    EdgeWalkAroundVertex it(topology, i);
    if(!it.IsOpen())
      {
      // Compute the laplacian of phi (TODO: is it needed?)
      double xLapPhi = 0.0;
      for(j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        xLapPhi += A.GetSparseData()[j] * xAtoms[A.GetColIndex()[j]].F;
      
      // The first component of the partial derivative is the derivative 
      // of the term (xFanArea). It is relatively easy to compute.
      double xRhoOverFanArea = xLapPhi / xVertexGeom[i].xFanArea;

      // Scaling factor for cotangent derivatives
      double xScale = 1.5 / xVertexGeom[i].xFanArea;
      
      // Reference to the weight for the fixed vertex
      SMLVec3d &wii = W[xMapVertexToA[i]];
      wii.fill(0.0);

      // Get the value of phi at the center
      double phiFixed = a.F;

      for( ; !it.IsAtEnd(); ++it)
        {
        // The weight vector for the moving vertex
        SMLVec3d &wij = W[xMapVertexNbrToA[it.GetPositionInMeshSparseArray()]];

        // Get the triangle index in front and behind
        size_t ta = it.TriangleAhead(), tb = it.TriangleBehind();
        size_t va = it.MovingVertexIndexInTriangleAhead();
        size_t vb = it.MovingVertexIndexInTriangleBehind();
        size_t qa = it.OppositeVertexIndexInTriangleAhead();
        size_t qb = it.OppositeVertexIndexInTriangleBehind();
        size_t pa = it.FixedVertexIndexInTriangleAhead();
        size_t pb = it.FixedVertexIndexInTriangleBehind();

        // Get the phi's at the surrounding vertices
        double phiAhead = xAtoms[it.VertexIdAhead()].F;
        double phiBehind = xAtoms[it.VertexIdBehind()].F;
        double phiMoving = xAtoms[it.MovingVertexId()].F;

        // Assign the first component of the weight
        wij = (phiMoving - phiFixed) *
          (xTriangleGeom[ta].xCotGrad[qa][va] + xTriangleGeom[tb].xCotGrad[qb][vb]);
        wij += (phiAhead - phiFixed) * xTriangleGeom[ta].xCotGrad[va][va];
        wij += (phiBehind - phiFixed) * xTriangleGeom[tb].xCotGrad[vb][vb];

        // Scale the vector accumulated so far
        wij *= xScale;

        // Assign the second half of the value to the vector
        wij -= xRhoOverFanArea * 
          (xTriangleGeom[ta].xAreaGrad[va] + xTriangleGeom[tb].xAreaGrad[vb]);

        // The cental value can be computed based on the fact that all weights must
        // add up to zero, so the sum of the derivatives should also be zero
        // wii += xScale * (phiMoving - phiFixed) *
        //  (xTriangleGeom[ta].xCotGrad[qa][pa] + xTriangleGeom[tb].xCotGrad[qb][pb]);
        // wii -= xRhoOverFanArea * xTriangleGeom[ta].xAreaGrad[pa];
        wii -= wij;
        }
      }
    else
      {
      // In this step, we must compute the weights of nu_j (j in Nhd(i)) as
      // they contribute to the right hand side of the variational PDE. To see
      // the derivation, look in the March 27, 2006 notebook (red), page with three
      // stars in the top right corner.
      
      // Get references to the differential geometry
      double (*g)[2] = a.G.xContravariantTensor;

      // Compute the tensor h
      double p1 = g[0][0] * a.Fu + g[1][0] * a.Fv;
      double p2 = g[0][1] * a.Fu + g[1][1] * a.Fv;
      double h11 = - 2 * p1 * p1, h12 = - 2 * p1 * p2, h22 = - 2 * p2 * p2;
      
      // Compute the vectors Y1 and Y2
      SMLVec3d Y1 = h11 * a.Xu + h12 * a.Xv;
      SMLVec3d Y2 = h12 * a.Xu + h22 * a.Xv;

      // Add up to get the weights in the sparse matrix
      for(EdgeWalkAroundVertex it(topology, i); !it.IsAtEnd(); ++it)
        W[xMapVertexNbrToA[it.GetPositionInMeshSparseArray()]] = 
          xLoopScheme.GetNeighborWeight(0, it) * Y1 +
          xLoopScheme.GetNeighborWeight(1, it) * Y2;

      // Finally, add the weight for the point itself (-4 \epsilon)
      W[xMapVertexToA[i]] = 
        xLoopScheme.GetOwnWeight(0, i) * Y1 + 
        xLoopScheme.GetOwnWeight(1, i) * Y2;
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

  // We must initialize pardiso
  xPardiso.SymbolicFactorization(
    topology->nVertices, xPardisoRowIndex, xPardisoColIndex, A.GetSparseData());

  // In the subsequent iterations, only do the numeric factorization and solve
  xPardiso.NumericFactorization(A.GetSparseData());

  // Precompute common terms for atom derivatives
  xTempDerivativeTerms = TempDerivativeTermsArray(topology->nVertices);
  for(i = 0; i < topology->nVertices; i++)
    xAtoms[i].ComputeCommonDerivativeTerms(xTempDerivativeTerms[i]);
}

void
MeshMedialPDESolver
::ComputeAtomVariationalDerivative(MedialAtom *dAtoms)
{
  size_t i;

  // Each variational derivative is a linear system Ax = b. The matrix A
  // should be set correctly from the last call to Solve() and all we have to
  // change are the vectors B.
  vnl_vector<double> rhs(topology->nVertices, 0.0);
  vnl_vector<double> soln(topology->nVertices, 0.0);

  // Repeat for each atom
  for(i = 0; i < topology->nVertices; i++) 
    {
    // Set up the initial right hand side value
    rhs[i] = (topology->IsVertexInternal(i)) ? dAtoms[i].xLapR : 0.0;

    // Compute the rest of the right hand side
    for(size_t j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
      rhs[i] -= dot_product(W[j], dAtoms[A.GetColIndex()[j]].X);
    }

  // Solve the partial differential equation (dPhi/dVar)
  xPardiso.Solve(rhs.data_block(), soln.data_block());

  // Compute each atom
  for(i = 0; i < topology->nVertices; i++) 
    {
    // For each atom, compute F, Fu, Fv, Xu and Xv
    MedialAtom &da = dAtoms[i];
    da.F = soln[i];
    da.Fu = xLoopScheme.GetPartialDerivative(0, i, soln.data_block());
    da.Fv = xLoopScheme.GetPartialDerivative(1, i, soln.data_block());
    da.Xu = xLoopScheme.Xu(i, dAtoms);
    da.Xv = xLoopScheme.Xv(i, dAtoms);
    da.Xuu = da.Xuv = da.Xvv = SMLVec3d(0.0); 

    // Compute the metric tensor derivatives of the atom
    xAtoms[i].ComputeMetricTensorDerivatives(dAtoms[i]);

    // Compute the derivatives of the boundary nodes
    xAtoms[i].ComputeBoundaryAtomDerivatives(dAtoms[i], xTempDerivativeTerms[i]);
    }
}

#endif

/*
void
MeshMedialPDESolver
::ComputeGradient(vector<MedialAtom *> xVariations)
{
  size_t i, q, j;

  // Compute the weight matrix for the right hand side computations
  this->ComputeRHSGradientMatrix();

  // Each variational derivative is a linear system Ax = b. The matrix A
  // should be set correctly from the last call to Solve() and all we have to
  // change are the vectors B.
  vnl_vector<double> rhs(topology->nVertices, 0.0);
  vnl_vector<double> soln(topology->nVertices, 0.0);

  // We must initialize pardiso
  xPardiso.SymbolicFactorization(
    topology->nVertices, xPardisoRowIndex, xPardisoColIndex, A.GetSparseData());

  // In the subsequent iterations, only do the numeric factorization and solve
  xPardiso.NumericFactorization(A.GetSparseData());

  // Precompute common terms for atom derivatives
  MedialAtom::DerivativeTerms *dt = new MedialAtom::DerivativeTerms[topology->nVertices];
  for(i = 0; i < topology->nVertices; i++)
    xAtoms[i].ComputeCommonDerivativeTerms(dt[i]);

  // Compute the derivative for each variation
  for(q = 0; q < xVariations.size(); q++)
    {
    // Get the derivative atoms for this variation
    MedialAtom *dAtoms = xVariations[q];

    // Repeat for each atom
    for(i = 0; i < topology->nVertices; i++) 
      {
      // Set up the initial right hand side value
      rhs[i] = (topology->IsVertexInternal(i)) ? dAtoms[i].xLapR : 0.0;

      // Compute the rest of the right hand side
      for(size_t j = A.GetRowIndex()[i]; j < A.GetRowIndex()[i+1]; j++)
        rhs[i] -= dot_product(W[j], dAtoms[A.GetColIndex()[j]].X);
      }

    // Solve the partial differential equation (dPhi/dVar)
    xPardiso.Solve(rhs.data_block(), soln.data_block());

    // Compute each atom
    for(i = 0; i < topology->nVertices; i++) 
      {
      // For each atom, compute F, Fu, Fv, Xu and Xv
      MedialAtom &da = dAtoms[i];
      da.F = soln[i];
      da.Fu = xLoopScheme.GetPartialDerivative(0, i, soln.data_block());
      da.Fv = xLoopScheme.GetPartialDerivative(1, i, soln.data_block());
      da.Xu = xLoopScheme.Xu(i, dAtoms);
      da.Xv = xLoopScheme.Xv(i, dAtoms);
      da.Xuu = da.Xuv = da.Xvv = SMLVec3d(0.0); 

      // Compute the metric tensor derivatives of the atom
      xAtoms[i].ComputeMetricTensorDerivatives(dAtoms[i]);

      // Compute the derivatives of the boundary nodes
      xAtoms[i].ComputeBoundaryAtomDerivatives(dAtoms[i], dt[i]);
      }
    }

  // Delete junk
  delete dt;
}
*/


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
