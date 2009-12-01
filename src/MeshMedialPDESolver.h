#ifndef __MeshMedialPDESolver_h_
#define __MeshMedialPDESolver_h_

#include "MedialAtom.h"
#include "MedialAtomGrid.h"
#include "SparseMatrix.h"
#include "SubdivisionSurface.h"
#include "GenericMedialModel.h"
#include "SparseSolver.h"
#include <smlmath.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>


/**
 * Mesh Medial PDE solver is similar to the regular medial PDE solver with
 * some substantial differences. The domain of the representation is a
 * triangular mesh of arbitrary shape. Finite difference expressions are
 * computed using mesh-based operators rather than differential geometry
 *
 * Solving an equation involves solving a series of linear Newton step
 * partial differential equations:
 * \Delta \epsilon_i = \rho - \Delta \phi_i, subject to
 * 2 \Nabla \epsilon_i \cdot \Nabla \phi_i - 4 \epsilon_i =
 *            = 4 \phi_i - \| \Nabla \phi_i \|^2 on \partial \Omega
 * To solve this system, we have to construct a matrix that contains the
 * weight of each \epsilon in each of the equations. This is a sparse matrix
 * where each vertex represents a row, and has non-zero values in the
 * columns that correspond to the adjacent vertices
 */
class MeshMedialPDESolver
{
public:
  // Matrix typedefs
  typedef ImmutableSparseMatrix<double> SparseMat;
  typedef SubdivisionSurface::MeshLevel MeshLevel;
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;

  // Constructor
  MeshMedialPDESolver();

  // Destructor
  ~MeshMedialPDESolver();

  // Set the topology of the mesh. This determines the connectivity of the
  // vertices which, in turn, determines the structure of the sparse matrix
  // passed to the sparse solver
  void SetMeshTopology(MeshLevel *topology, MedialAtom *managedAtoms);

  // Our first attempt at a solver method. The flag specifies whether the
  // terms involved in gradient computation should also be computed
  void SolveEquation(bool flagGradient = false, bool flagAllowErrors = false);

  // Compute the common part of the gradient computation
  void BeginGradientComputation();

  // Compute the directional derivative of the solution
  void ComputeAtomVariationalDerivative(MedialAtom * dAtoms);

  // Test the accuracy of partial derivative computations in the gradient code
  // this method should be called after calling solve with some data
  int TestPartialDerivatives();

  // Test whether the Jacobian is computed correctly
  int TestJacobianAndGradient(double *xInitSoln);

  // Get the array of atoms
  MedialAtom *GetAtomArray() const
    { return xAtoms; }

  // This method computes the LBO at each internal node. 
  // vnl_vector<double> ComputeLBO(const double *phi);

private:

  // This method resets all the pointers associated with a mesh
  void Reset();

  // This method is used to compute the sparse matrix A for Newton's method
  void FillSparseMatrix(bool flagInputChange);

  // This method is used to compute the right hand side B for Newton's method
  void FillRHS();

  // This computes the geometry associated with a mesh before running the solver
  void ComputeMeshGeometry(bool flagGradient);

  // This computes the medial atoms once the phi has been solved for
  void ComputeMedialAtoms(const double *soln);

  // Compute the weight matrix used for gradient computations
  void ComputeRHSGradientMatrix();

  // Compute the condition number of the Jacobian
  void ComputeJacobianConditionNumber();


  // An immutable sparse matrix used to represent the PDE. Each row
  // corresponds to a vertex, non-zero values are found between adjacent
  // vertices (those contributing to gradient and LBO computations)
  SparseMat M;

  // Right hand side and the solution of the PDE
  vnl_vector<double> xRHS, xSolution;

  // A pointer to the mesh topology
  MeshLevel *topology;

  // A structure representing triangle geometry in the mesh
  struct TriangleGeom
    {
    // Area of the triangle
    double xArea;

    // Cotangent of the three angles
    double xCotangent[3];

    // The gradient of the area with respect to each vertex
    SMLVec3d xAreaGrad[3];

    // The gradient of the cotangent with respect to each vertex
    SMLVec3d xCotGrad[3][3];
    };

  // A structure representing vertex geometry in the mesh (temp)
  struct VertexGeom
    {
    // Area of the triangle fan around the vertex
    double xFanArea;

    // The two Loop tangent vectors at the vertex
    // SMLVec3d t1, t2;

    // The first fundamental form and its inverse
    // double gCovariant[2][2], gContravariant[2][2];
    };

  // An index object, used to map one sparse array to another
  class MeshMatrixXRef
    {
    public:
    std::vector<size_t> xSelfIndex, xNbrIndex;
    void Initialize(size_t nVertices, size_t nEdges)
      { 
      xSelfIndex.resize(nVertices, 0);
      xNbrIndex.resize(nEdges, 0);
      }
    };

  // There are 3 x-references, for the 3 quadrants of M that change
  MeshMatrixXRef xIndexAPhi, xIndexAOmega, xIndexN;

  // Geometry arrays that store triangle-related and vertex-related info
  TriangleGeom *xTriangleGeom;
  VertexGeom *xVertexGeom;

  // Scheme for computing tangent vectors
  MedialAtomLoopScheme xLoopScheme;

  // For gradient computation, there is an array W, which contains the
  // Jacobian of each finite difference equation with respect to each of the
  // neighbor atoms's positions. This is an array of vectors
  SMLVec3d *WX;

  // For gradient computation of the Neumann boundary condition, we have 
  // some weights for the derivatives of differential geometric properties
  struct NeumannDerivWeights 
    { double wgt_g12, wgt_g22, wgt_g, wgt_f, wgt_fv; } *Wfu;
  

  // Array of medial atoms managed by this solver (should it be managed
  // elsewhere?)
  MedialAtom *xAtoms;

  // Pardiso solver
  SparseSolver *xSolver;

  // LM optimizer callbacks
  static void ComputeLMResidual(void *handle, int n, double *x, double *fx);
  static void ComputeLMJacobian(void *handle, int n, double *x, SparseMat &J);

  // Array that holds temporary data for gradient computation
  typedef std::vector<MedialAtom::DerivativeTerms> TempDerivativeTermsArray;
  TempDerivativeTermsArray xTempDerivativeTerms;
};



#endif
