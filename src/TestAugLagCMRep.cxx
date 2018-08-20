/**
 * This is a test of the new cm-rep approach that combines geodesic shooting
 * and boundary-based medial constraints
 */

#include "MeshTraversal.h"
#include "MedialException.h"
#include "PointSetOptimalControlSystem.h"
#include "PointSetHamiltonianSystem.h"
#include "SparseMatrix.h"

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtksys/SystemTools.hxx>

#include <vnl_matrix.h>
#include <vnl_vector.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_lbfgs.h>
#include <vnl/algo/vnl_lbfgsb.h>
#include <vnl/algo/vnl_conjugate_gradient.h>
#include <vnl/vnl_random.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_brent_minimizer.h>

class FindMedialTriangleCenterObjective : public vnl_cost_function
{
public:
  typedef vnl_vector<double> Vector;
  FindMedialTriangleCenterObjective(vnl_matrix<double> X) 
    : vnl_cost_function(3)
    { 
    this->X = X;
    }

  // Get initial solution
  vnl_vector<double> get_init()
    {
      vnl_vector<double> Y(X.columns(), 0.0);
      for(int i = 0; i < X.rows(); i++)
        Y += X.get_row(i);

      return Y / X.rows();
    }

  // Compute the variance
  virtual double f (vnl_vector<double> const& Y)
    {
    // Triangle vertice
    Vector A1 = X.get_row(0);
    Vector B1 = X.get_row(1);
    Vector C1 = X.get_row(2);
    Vector A2 = X.get_row(3);
    Vector B2 = X.get_row(4);
    Vector C2 = X.get_row(5);

    // Spoke vectors
    Vector AB1 = (B1 - A1), AC1 = (C1 - A1);
    Vector AB2 = (B2 - A2), AC2 = (C2 - A2);
    Vector S1 = (A1 + B1 + C1) / 3 - Y;
    Vector S2 = (A2 + B2 + C2) / 3 - Y;

    // Compute the total error
    double c0 = dot_product(AB1, S1);
    double c1 = dot_product(AC1, S1);
    double c2 = dot_product(AB2, S2);
    double c3 = dot_product(AC2, S2);
    double c4 = dot_product(S1 - S2, S1 + S2);

    return c0 * c0 + c1 * c1 + c2 * c2 + c3 * c3 + 2 * c4 * c4;
    }

  virtual void gradf(vnl_vector<double> const& x, vnl_vector<double>& g)
    {
    this->fdgradf(x, g);
    }

protected:
  vnl_matrix<double> X;
};

vnl_vector<double> FindMedialTriangleCenter(const vnl_matrix<double> &bnd_vertices)
{
  FindMedialTriangleCenterObjective objective(bnd_vertices);
  vnl_vector<double> Y = objective.get_init();

  vnl_lbfgsb optimizer(objective);
  optimizer.set_f_tolerance(1e-9);
  optimizer.set_x_tolerance(1e-4);
  optimizer.set_g_tolerance(1e-6);
  optimizer.set_trace(true);
  optimizer.set_max_function_evals(100);
  optimizer.minimize(Y);

  return Y;
}

/** 
 * A medial mesh representation
 */
struct CMRep 
{
  typedef vnl_vector<double> Vector;
  typedef vnl_matrix<double> Matrix;

  typedef vnl_vector<unsigned int> IdxVector;
  typedef vnl_matrix<unsigned int> IdxMatrix;

  typedef ImmutableSparseMatrix<double> SparseMat;
  typedef SparseMat::RowIterator SparseRowIter;

  // The boundary mesh (VTK)
  vtkSmartPointer<vtkPolyData> bnd_vtk;

  // Index of boundary medial indices and boundary triangle medial indices
  IdxVector bnd_mi, bnd_mti;

  // Index of boundary triangles corresponding to medial triangles
  IdxMatrix med_bti;

  // The half-edge structure of the boundary mesh
  TriangleMesh bnd_tri;

  // Half-edge structure of the medial mesh
  TriangleMesh med_tri;

  // Boundary and medial vertices
  Matrix bnd_vtx, med_vtx;

  // For each medial vertex, the list of corresponding boundary ones
  std::vector< std::vector<unsigned int> > med_bi;

  // Sparse matrices used to compute tangent vectors Qu, Qv from Q's
  SparseMat wgt_Quv[2];

  // Number of boundary vertices
  unsigned int nv, nmv;

  // Number of boundary and medial triangles
  unsigned int nt, nmt;

  // Read from a VTK file (bcm-rep format)
  void ReadVTK(const char *file);
};

void CMRep::ReadVTK(const char *fn)
{
  // Read the mesh from file
  vtkSmartPointer<vtkPolyDataReader> reader = vtkPolyDataReader::New();
  reader->SetFileName(fn);
  reader->Update();
  this->bnd_vtk = reader->GetOutput();

  // Read the normals
  vtkDataArray *b_norm = this->bnd_vtk->GetPointData()->GetNormals();
  vtkDataArray *b_rad = this->bnd_vtk->GetPointData()->GetArray("Radius");

  // Get the number of vertices and triangles
  nv = this->bnd_vtk->GetNumberOfPoints();
  nt = this->bnd_vtk->GetNumberOfCells();
  nmt = nt / 2;

  // Read the point coordinates into a matrix and the medial index array
  this->bnd_mi.set_size(nv);
  this->bnd_vtx.set_size(nv, 3);
  vtkDataArray *da_mi = this->bnd_vtk->GetPointData()->GetArray("MedialIndex");
  for(int k = 0; k < nv; k++)
    {
    this->bnd_mi[k] = da_mi->GetTuple1(k);
    for(int d = 0; d < 3; d++)
      this->bnd_vtx(k, d) = this->bnd_vtk->GetPoint(k)[d];
    }

  // Parse the boundary triangles
  TriangleMeshGenerator bnd_tmg(&this->bnd_tri, nv);
  for(int i = 0; i < nt; i++)
    {
    vtkCell *cell = this->bnd_vtk->GetCell(i);
    if(cell->GetNumberOfPoints() != 3)
      throw ModelIOException("Non-triangle cell in input mesh");

    bnd_tmg.AddTriangle(cell->GetPointId(0), cell->GetPointId(1), cell->GetPointId(2));
    }
  bnd_tmg.GenerateMesh();

  // Match triangles between the two sides of the mesh. This can be done by hashing on
  // the medial index array. So for each triangle, we compute a unique number representing
  // the medial indices of its vertices, and use this to match the boundary triangles
  this->bnd_mti.set_size(nt);
  this->med_bti.set_size(nmt, 2);
  typedef std::map<long long, unsigned int> TriHashMap;
  TriHashMap tri_hash_to_idx;
  unsigned int mt_idx = 0;
  for(unsigned int i = 0; i < nt; i++)
    {
    // Get the medial indices
    long long mi[3];
    for(int j = 0; j < 3; j++)
      mi[j] = this->bnd_mi[this->bnd_tri.triangles[i].vertices[j]];

    // Sort the medial indices
    std::sort(mi, mi+3);

    // Compute a hash
    long long f = 1000000;
    long long tri_hash = f * f * mi[0] + f * mi[1] + mi[2];

    // Associate the hash with the index
    TriHashMap::const_iterator it = tri_hash_to_idx.find(tri_hash);
    if(it == tri_hash_to_idx.end())
      {
      tri_hash_to_idx[tri_hash] = mt_idx;
      this->bnd_mti[i] = mt_idx;
      this->med_bti(mt_idx, 0) = i;
      mt_idx++;
      }
    else
      {
      // Found the corresponding boundary triangle
      this->bnd_mti[i] = it->second;
      this->med_bti(it->second, 1) = i;
      }
    }

  // Compute medial to boundary mappings
  this->nmv = this->bnd_mi.max_value() + 1;
  this->med_bi.resize(nmv);
  for(int i = 0; i < this->nv; i++)
    this->med_bi[this->bnd_mi[i]].push_back(i);

  // Compute initial medial vertices. These can actually be imported from the model
  this->med_vtx.set_size(nmv, 3);
  for(unsigned int i = 0; i < nmv; i++)
    {
    // Get the first boundary point
    unsigned int j = this->med_bi[i][0];

    // Compute the medial point coordinate
    for(int a = 0; a < 3; a++)
      this->med_vtx(i, a) = this->bnd_vtx(j, a) - b_norm->GetComponent(j, a) * b_rad->GetTuple1(j);
    }


  // Compute the weights used to generate Qu and Qv derivatives
  LoopTangentScheme lts;
  lts.SetMesh(&this->bnd_tri);
  
  for(int a = 0; a < 2; a++)
    {
    vnl_sparse_matrix<double> vnl_wgt_Quv(this->nv, this->nv);
    for(int k = 0; k < this->nv; k++)
      {
      vnl_wgt_Quv.put(k, k, lts.GetOwnWeight(a, k));
      for(EdgeWalkAroundVertex walk(&this->bnd_tri, k); !walk.IsAtEnd(); ++walk)
        vnl_wgt_Quv.put(k, walk.MovingVertexId(), lts.GetNeighborWeight(a, walk));
      }

    // Place into an efficient sparse matrix
    this->wgt_Quv[a].SetFromVNL(vnl_wgt_Quv);
    }
}

struct AugLagMedialFitParameters
{
  // Number of timesteps for geodesic shooting
  unsigned int nt;

  // Weight of the kinetic energy term
  double w_kinetic;

  // Sigma of the Gaussian kernel
  double sigma;
  
  // Schedule for mu
  double mu_init, mu_scale;
  unsigned int gradient_iter;
  unsigned int newton_iter;

  // Interpolation mode or fitting mode?
  bool interp_mode;

  // Do we do derivative checks?
  bool check_deriv;

  // Default initializer
  AugLagMedialFitParameters() 
    : nt(40), w_kinetic(0.05), sigma(4.0), 
      mu_init(0.01), mu_scale(1.6), gradient_iter(6000), newton_iter(100),
      interp_mode(false), check_deriv(true)  {}
};



/**
 * Specific implementation of constraints
 */
class PointBasedMediallyConstrainedFittingTraits
{
public:
  typedef vnl_vector<double> Vector;
  typedef vnl_matrix<double> Matrix;
  typedef vnl_matrix_ref<double> MatrixRef;
  typedef vnl_vector_ref<double> VectorRef;

  static int GetNumberOfActiveVertices(CMRep *model)
    {
    return (model->nv + model->nmv);
    }

  static int GetNumberOfConstraintsPerTimepoint(CMRep *model)
    {
    // The constraint set:
    //   - each spoke is orthogonal to the boundary surface (2 * nv_b)
    //   - spokes are equal length (nm_bitangent + 2 * nm_tritangent)
    int nc_t = 2 * model->nv;
    for(int j = 0; j < model->med_bi.size(); j++)
      nc_t += model->med_bi[j].size() - 1;
    return nc_t;
    }

  static void ComputeInitialControl(CMRep *model, CMRep *target, Vector &u, unsigned int nt)
    {
    MatrixRef u_bnd(model->nv, 3, const_cast<double *>(u.data_block()));
    MatrixRef u_med(model->nmv, 3, const_cast<double *>(u.data_block()) + model->nv * 3);

    u_bnd.update((target->bnd_vtx - model->bnd_vtx) / nt, 0, 0);
    u_med.update((target->med_vtx - model->med_vtx) / nt, 0, 0);
    }

  static void ComputeInitialLandmarks(CMRep *model, Matrix &q0)
    {
    q0.update(model->bnd_vtx, 0, 0);
    q0.update(model->med_vtx, model->nv, 0);
    }

  static double ComputeAugmentedLagrangianObjective(
    CMRep *model, CMRep *target, Matrix &qt, Matrix &d_g__d_qt)
    {
    // Get references to the boundary and medial points
    MatrixRef qt_bnd(model->nv, 3, const_cast<double *>(qt.data_block()));

    // Get the partial derivatives for the current timeslice
    MatrixRef d_g__d_qbnd_t(model->nv, 3, d_g__d_qt.data_block());

    // Total objective
    double m_distsq = 0.0;
    for(int i = 0; i < model->nv; i++)
      {
      for(int d = 0; d < 3; d++)
        {
        double del_id = qt_bnd(i,d) - target->bnd_vtx(i,d);
        m_distsq += del_id * del_id;

        // Compute the contribution to d_g__d_q1
        d_g__d_qbnd_t(i,d) += 2 * del_id; 
        }
      }

    return m_distsq;
    }

  struct QuadraticForm 
    {
    ImmutableSparseMatrix<double> A;
    Vector b;
    double c;

    virtual void Initialize(vnl_sparse_matrix<double> &inA, const Vector &in_b, double in_c)
      {
      A.SetFromVNL(inA);
      b = in_b;
      c = in_c;
      }

    /** Compute 0.5 x^t A x + b^t x + c */
    virtual double Compute(const Vector &x, Vector &gradient)
      {
      Vector Ax = A.MultiplyByVector(x);
      gradient = Ax + b;
      return dot_product(x, Ax * 0.5 + b) + c;
      }

    QuadraticForm()
      {
      b.fill(0.0); c = 0.0;
      }
    };

  // Expressions like (M x + d) ^t (M x + d)
  struct SymmetricQuadraticForm : public QuadraticForm
    {
    ImmutableSparseMatrix<double> M;
    Vector d;

    virtual void Initialize(vnl_sparse_matrix<double> &inM, const Vector &in_d)
      {
      M.SetFromVNL(inM);
      d = in_d;

      // Compute the Hessian matrix
      vnl_sparse_matrix<double> H = 2.0 * (inM.transpose() * inM);
      this->A.SetFromVNL(H);

      // Compute linear terms
      this->b = M.MultiplyByVector(d) + M.MultiplyTransposeByVector(d);
      this->c = dot_product(d, d);
      }

    virtual double Compute(const Vector &x, Vector &gradient)
      {
      Vector Mx_d = M.MultiplyByVector(x) + d;
      gradient = M.MultiplyByVector(Mx_d) * 2.0;
      return dot_product(Mx_d, Mx_d);
      }

    };

  // In this problem, the objective and the constraints have the quadratic form
  //     x^t * A * x + b^t * x + c
  // with very sparse A and b
  //
  // This data structure holds the A, b and c for the main objective and the constraints
  struct HessianData 
    {
    // Quadratic parameters of each of the constraints
    std::vector<QuadraticForm> qf_C;

    // Quadratic parameters of the objective function
    SymmetricQuadraticForm qf_F;

    // Hessian of the augmented lagrangian (cached for computations)
    ImmutableSparseMatrix<double> HL_init, HL;
    };


  static void ComputeAugmentedLagrangianConstraints(
    CMRep *model, const Matrix &qt, Matrix &d_g__d_qt, Vector &C, Vector &lambda, double mu)
    {
    // Get references to the boundary and medial points
    MatrixRef qt_bnd(model->nv, 3, const_cast<double *>(qt.data_block()));
    MatrixRef qt_med(model->nmv, 3, const_cast<double *>(qt.data_block()) + model->nv * 3);

    // Get references to the constraints in the list C
    MatrixRef Ct_orth(model->nv, 2, C.data_block());
    VectorRef Ct_spklen(C.size() - 2 * model->nv, C.data_block() + 2 * model->nv);

    // Get the lambdas for the current timeslice
    MatrixRef Lt_orth(model->nv, 2, lambda.data_block());
    VectorRef Lt_spklen(C.size() - 2 * model->nv, lambda.data_block() + 2 * model->nv);

    // Get the partial derivatives for the current timeslice
    MatrixRef d_g__d_qbnd_t(model->nv, 3, d_g__d_qt.data_block());
    MatrixRef d_g__d_qmed_t(model->nmv, 3, d_g__d_qt.data_block() + 3 * model->nv);

    // For each boundary point, compute its directional derivatives using the loop scheme
    for(int j = 0; j < model->nv; j++)
      {
      // Compute the spoke
      vnl_vector_fixed<double, 3> S(0.0);
      for(int a = 0; a < 3; a++)
        S[a] = qt_bnd(j, a) - qt_med(model->bnd_mi(j), a);

      // Iterate over the two tangent vectors at each vertex
      for(int d = 0; d < 2; d++)
        {
        // Initialize the tangent vector to zero
        vnl_vector_fixed<double, 3> Qd(0.0);

        // Compute the tangent vector as weighted sum of neighbors
        for(CMRep::SparseRowIter it = model->wgt_Quv[d].Row(j); !it.IsAtEnd(); ++it)
          for(int a = 0; a < 3; a++)
            Qd[a] += qt_bnd(it.Column(), a) * it.Value();

        // Compute the dot product with the spoke - this should be zero
        Ct_orth(j, d) = dot_product(S, Qd);

        // Compute the partial of this constraint's contribution to the objective
        // with respect to each q_i
        double d_g__d_C = mu * Ct_orth(j, d) - Lt_orth(j, d);

        // Contribution to the ring of vertices
        for(CMRep::SparseRowIter it = model->wgt_Quv[d].Row(j); !it.IsAtEnd(); ++it)
          {
          double z = d_g__d_C * it.Value();
          for(int a = 0; a < 3; a++)
            d_g__d_qbnd_t(it.Column(), a) += z * S[a];
          }

        // Contribution to the medial and boundary vertices
        for(int a = 0; a < 3; a++)
          {
          double z = d_g__d_C * Qd[a];
          d_g__d_qbnd_t(j, a) += z;
          d_g__d_qmed_t(model->bnd_mi[j], a) -= z;
          }
        }
      }

    // Iterate over the medial points now and compute spoke equality constraints
    //   k is index over medial vertices, m is index over constraints
    // TODO: for now we just consider bitangencies, should consider all cases
    for(int k = 0, m = 0; k < model->med_bi.size(); k++)
      {
      if(model->med_bi[k].size() > 1)
        {
        vnl_vector_fixed<double, 3> del0, delj;

        // Compute the length of the first spoke
        int b0 = model->med_bi[k][0];
        for(int a = 0; a < 3; a++)
          del0[a] = qt_bnd(b0,a) - qt_med(k, a);
        double len2_0 = del0.squared_magnitude();

        // Compute the lengths of all other spokes and define constraints
        for(int j = 1; j < model->med_bi[k].size(); j++, m++)
          {
          // Compute the length
          int bj = model->med_bi[k][j];
          for(int a = 0; a < 3; a++)
            delj[a] = qt_bnd(bj,a) - qt_med(k, a);
          double len2_j = delj.squared_magnitude();

          // Define the constraint
          Ct_spklen[m] = len2_j - len2_0;

          // Compute the partials of the constraint with respect to the q's involved
          double d_g__d_C = mu * Ct_spklen[m] - Lt_spklen[m];
          double z = 2.0 * d_g__d_C;

          for(int a = 0; a < 3; a++)
            {
            // Contribution to boundary nodes
            d_g__d_qbnd_t(bj, a) += z * delj[a];
            d_g__d_qbnd_t(b0, a) -= z * del0[a];
            d_g__d_qmed_t(k, a) += z * (del0[a] - delj[a]);
            }
          }
        }
      }
    }

  /** Compute the Augmented Lagrangian and its derivative using precomputed quantities */
  static double ComputeAugmentedLagrangianAndGradient(
    CMRep *model, const Vector &qt, Vector &d_AL__d_qt, Vector &C, Vector &lambda, double mu, HessianData *H)
    {
    // Initialize the outputs with the f/grad of the function f
    double AL = H->qf_F.Compute(qt, d_AL__d_qt); 

    // Vector containing the gradient of the constraint - to be reused
    Vector grad_Cj(qt.size(), 0.0);

    // Iterate over the constraints
    for(int j = 0; j < H->qf_C.size(); j++)
      {
      // Reference to the current quadratic form
      QuadraticForm &Z = H->qf_C[j];

      // Compute the constraint and its gradient
      C[j] = Z.Compute(qt, grad_Cj);

      // Add to the overall objective
      AL += C[j] * (C[j] * mu / 2 - lambda[j]);

      // Add to the overall gradient
      d_AL__d_qt += grad_Cj * (mu * C[j] - lambda[j]);
      }

    return AL;
    }


  /** Compute the Augmented Lagrangian, its derivative and its Hessian using precomputed quantities */
  static double ComputeAugmentedLagrangianJet(
    CMRep *model, const Vector &qt, Vector &d_AL__d_qt, Vector &C, Vector &lambda, double mu, HessianData *H)
    {
    // Get a reference to the Hessian of the Lagrangian
    ImmutableSparseMatrix<double> &HL = H->HL;

    // Initialize the outputs with the Jet of the function f
    double AL = H->qf_F.Compute(qt, d_AL__d_qt); 
    HL = H->HL_init;

    // Vector containing the gradient of the constraint - to be reused
    Vector grad_Cj(qt.size(), 0.0);

    // Iterate over the constraints
    for(int j = 0; j < H->qf_C.size(); j++)
      {
      // Reference to the current quadratic form
      QuadraticForm &Z = H->qf_C[j];

      // Compute the constraint and its gradient
      C[j] = Z.Compute(qt, grad_Cj);

      // Add to the overall objective
      AL += C[j] * (C[j] * mu / 2 - lambda[j]);

      // Add to the overall gradient
      d_AL__d_qt += grad_Cj * (mu * C[j] - lambda[j]);

      // Update the Hessian of the lagrangian
      HL.AddScaledMatrix(Z.A, mu * C[j] - lambda[j]);

      // Update with a scaled outer product
      HL.AddScaledOuterProduct(grad_Cj, grad_Cj, mu);
      }

    return AL;
    }







  /** 
   * Helper function to compute index into the sparse matrix
   */
  static unsigned int IB(CMRep *model, unsigned int i_vtx, unsigned int i_dim)
    {
    return i_vtx * 3 + i_dim;
    }

  static unsigned int IM(CMRep *model, unsigned int i_vtx, unsigned int i_dim)
    {
    return model->nv * 3 + i_vtx * 3 + i_dim;
    }

  /** 
   * This function precomputes different terms used to compute the Hessian that 
   * do not change between iterations (exploiting wide use of quadratic functions
   * in the objective
   */
  static void PrecomputeHessianData(CMRep *model, CMRep *target, HessianData *data)
    {
    // Dimensions of all the sparse matrices
    unsigned int np = model->nv + model->nmv;
    unsigned int k = np * 3;

    // Number of the constraints
    unsigned int nc = GetNumberOfConstraintsPerTimepoint(model), ic = 0;

    // Initialize the data
    data->qf_C.resize(nc);

    // Handle the boundary point ortho constraints
    for(int j = 0; j < model->nv; j++)
      {
      // Iterate over the two tangent vectors at each vertex
      for(int d = 0; d < 2; d++, ic++)
        {
        // Initialize the Hessian for this constraint
        vnl_sparse_matrix<double> H_Cj(k, k);

        for(CMRep::SparseRowIter it = model->wgt_Quv[d].Row(j); !it.IsAtEnd(); ++it)
          {
          for(int a = 0; a < 3; a++)
            {
            unsigned int ib_mov = IB(model, it.Column(), a);
            unsigned int ib_j = IB(model, j, a);
            unsigned int im_j = IM(model, model->bnd_mi(j), a);
            
            H_Cj(ib_mov, ib_j) -= it.Value(); 
            H_Cj(ib_j, ib_mov) -= it.Value();
            H_Cj(ib_mov, im_j) += it.Value();
            H_Cj(im_j, ib_mov) += it.Value();
            }
          }

        // Set the A matrix, the b and c are zeros
        data->qf_C[ic].Initialize(H_Cj, Vector(k, 0.0), 0.0);
        }
      }

    // Handle the spoke equal length constraints
    for(int i = 0; i < model->med_bi.size(); i++)
      {
      int b0 = model->med_bi[i][0];
      if(model->med_bi[i].size() > 1)
        {
        for(int j = 1; j < model->med_bi[i].size(); j++, ic++)
          {
          // Initialize the Hessian for this constraint
          vnl_sparse_matrix<double> H_Cj(k, k);

          int bj = model->med_bi[i][j];
          for(int a = 0; a < 3; a++)
            {
            unsigned int ib_0 = IB(model, b0, a);
            unsigned int ib_j = IB(model, bj, a);
            unsigned int im = IM(model, i, a);

            H_Cj(ib_j, ib_j) = 2.0;
            H_Cj(im, ib_j) = -2.0;
            H_Cj(ib_j, im) = -2.0;

            H_Cj(ib_0, ib_0) = -2.0;
            H_Cj(im, ib_0) = 2.0;
            H_Cj(ib_0, im) = 2.0;
            }

          data->qf_C[ic].Initialize(H_Cj, Vector(k, 0.0), 0.0);
          }
        }
      }

    // Handle the objective function - its Hessian is the identify over the boundary vertices
    vnl_sparse_matrix<double> M_f(k, k);
    Vector d_f(k, 0.0);

    for(int j = 0; j < model->nv; j++)
      {
      for(int a = 0; a < 3; a++)
        {
        unsigned int ib = IB(model, j, a);
        M_f(ib,ib) = 1.0;
        d_f[ib] = - target->bnd_vtx(j, a);
        }
      }

    data->qf_F.Initialize(M_f, d_f);


    // Initialize the Hessian matrix for the whole problem
    vnl_sparse_matrix<double> H_L(k, k);

    // Append the entries from the main objective Hessian
    for(int j = 0; j < k; j++)
      for(ImmutableSparseMatrix<double>::RowIterator it = data->qf_F.A.Row(j); !it.IsAtEnd(); ++it)
        H_L(j, it.Column()) = it.Value();

    // Append the entries for each constraint
    for(int i = 0; i < nc; i++)
      {
      // We need not just the Hessian of Cj but also the entries from H^T H because of the
      // outer product computation. The easiest is to determine all the unique rows of j and
      // all all combinations
      std::vector<int> nzrows;

      for(int j = 0; j < k; j++)
        {
        // This is the iterator
        ImmutableSparseMatrix<double>::RowIterator it = data->qf_C[i].A.Row(j);
        if(it.Size())
          {
          // Add this row to the list of nz rows
          nzrows.push_back(j);

          // Add the non-zero entries from the Hessian
          for(; !it.IsAtEnd(); ++it)
            H_L(j, it.Column()) += 0.0;
          }
        }

      // Add all the cross-product entries
      for(int a = 0; a < nzrows.size(); a++)
        for(int b = 0; b < nzrows.size(); b++)
          H_L(nzrows[a],nzrows[b]) += 0.0;
      }

    data->HL_init.SetFromVNL(H_L);
    data->HL = data->HL_init;
    printf("Hessian has dimension %d x %d with %d non-empty values\n", k, k, (int) data->HL.GetNumberOfSparseValues());
    }

  static void ExportTimepoint(CMRep *model, 
    const Matrix &qt, const Vector &C, const Vector &lambda, const char *fname)
    {
    // Get the boundary and medial vertices
    MatrixRef qt_bnd(model->nv, 3, const_cast<double *>(qt.data_block()));
    MatrixRef qt_med(model->nmv, 3, const_cast<double *>(qt.data_block()) + model->nv * 3);

    // Create a copy of the mesh in the model
    vtkSmartPointer<vtkPolyData> pd = vtkPolyData::New();
    pd->DeepCopy(model->bnd_vtk);

    // Update the points
    for(int i = 0; i < qt_bnd.rows(); i++)
      pd->GetPoints()->SetPoint(i, qt_bnd(i,0), qt_bnd(i,1), qt_bnd(i,2));

    // Get the constraint values for the current timeslice
    MatrixRef Ct_orth(model->nv, 2, const_cast<double *>(C.data_block()));
    VectorRef Ct_spklen(C.size() - 2 * model->nv, const_cast<double *>(C.data_block() + 2 * model->nv));

    // Create a point array for the constraints
    vtkSmartPointer<vtkFloatArray> arr_con = vtkFloatArray::New();
    arr_con->SetNumberOfComponents(3);
    arr_con->SetNumberOfTuples(model->nv);

    // The third component needs to be initialized to zeros b/c we skip over some vertices
    arr_con->FillComponent(2, 0.0f);

    // Assign the ortho constraints
    for(int j = 0; j < model->nv; j++)
      {
      arr_con->SetComponent(j, 0, Ct_orth(j, 0));
      arr_con->SetComponent(j, 1, Ct_orth(j, 1));
      }

    // Assign the spoke constraints
    for(int k = 0, m = 0; k < model->med_bi.size(); k++)
      {
      if(model->med_bi[k].size() > 1)
        {
        int b0 = model->med_bi[k][0];
        double consum = 0.0;
        for(int j = 1; j < model->med_bi[k].size(); j++, m++)
          {
          int bj = model->med_bi[k][j];
          arr_con->SetComponent(bj, 2, Ct_spklen[m]);
          consum += Ct_spklen[m];
          }
        arr_con->SetComponent(b0, 2, -consum);
        }
      }

    arr_con->SetName("Constraints");
    pd->GetPointData()->AddArray(arr_con);

    // Write the model
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
    writer->SetFileName(fname);
    writer->SetInputData(pd);
    writer->Update();
    }
};

template <class TMainObjective>
class BrentObjective : public vnl_cost_function
{
  typedef CMRep::Vector Vector;

public:
  BrentObjective(const Vector &x, const Vector &grad, TMainObjective *obj)
    : vnl_cost_function(1)
    {
    this->grad = grad;
    this->obj = obj;
    this->x = x;
    }

  double f(const vnl_vector<double>& alpha_vec)
     { 
     double alpha = alpha_vec[0];
     Vector xa = x + alpha * grad;
     double f=0.0;
     obj->compute(xa, &f, NULL);

     printf("Alpha: %12.4f,   f: %12.4f\n", alpha, f);
     return f;
     }

  TMainObjective *obj;
  Vector grad;
  Vector x;
};

template <class Traits>
class PointMatchingWithTimeConstraintsAugLagObjective : public vnl_cost_function
{
public:

  typedef vnl_matrix_ref<double> MatrixRef;
  typedef vnl_vector_ref<double> VectorRef;
  typedef CMRep::Matrix Matrix;
  typedef CMRep::Vector Vector;
  typedef PointSetOptimalControlSystem<double, 3> OCSystem;

  PointMatchingWithTimeConstraintsAugLagObjective(
    const AugLagMedialFitParameters &in_param, 
    CMRep *in_model, CMRep *in_target)
    : model(in_model), target(in_target), param(in_param), 
      nvtx(Traits::GetNumberOfActiveVertices(in_model)),
      nvar_t(Traits::GetNumberOfActiveVertices(in_model) * 3), 
      nvar_total(Traits::GetNumberOfActiveVertices(in_model) * 3 * in_param.nt), 
      vnl_cost_function(Traits::GetNumberOfActiveVertices(in_model) * 3 * in_param.nt)
    {
    // Compute the initial u-vector. 
    u_init.set_size(nvar_total);
    for(int t = 0; t < param.nt; t++)
      {
      VectorRef u_init_t(nvar_t, u_init.data_block() + nvar_t * t);
      Traits::ComputeInitialControl(model, target, u_init_t, param.nt);
      }

    // Initialize the initial landmarks
    q0.set_size(nvtx, 3);
    Traits::ComputeInitialLandmarks(model, q0);

    // Initialize the hamiltonian system
    ocsys = new OCSystem(q0, param.sigma, param.nt);

    // There is a control (and constraints) at every time point
    u.resize(param.nt, Matrix(nvtx, 3));
    d_g__d_ut.resize(param.nt, Matrix(nvtx, 3));
    d_g__d_qt.resize(param.nt, Matrix(nvtx, 3));

    // Initialize mu based on parameters
    mu = param.mu_init;

    // The number of constraints per timepoint and total
    nc_t = Traits::GetNumberOfConstraintsPerTimepoint(model);
    nc_total = nc_t * param.nt;

    // Initialize the constraints matrices
    C.set_size(nc_total); C.fill(0.0);

    // Initialize the lambdas
    lambda.set_size(nc_total); lambda.fill(0.0);

    // Verbosity
    verbose = false;
    }

  /** Destructor */
  ~PointMatchingWithTimeConstraintsAugLagObjective() 
    {
    delete ocsys;
    }

  /** Get suitable x for initialization */
  Vector get_xinit() const { return u_init; }

  /** Get number of variables */
  unsigned int get_nvar() const { return u_init.size(); }

  void set_verbose(bool flag) { this->verbose = flag; }

  /** Compute the objective function */
  virtual void compute(const CMRep::Vector &x, double *f, CMRep::Vector *g)
    {
    // TODO: this copy-in is a waste of time and space
    for(int t = 0; t < param.nt; t++)
      u[t].copy_in(x.data_block() + nvar_t * t);

    // Initialize the array of objective function derivatives w.r.t. q_t
    for(int t = 0; t < param.nt; t++)
      d_g__d_qt[t].fill(0.0);

    // Initialize the various components of the objective function
    double m_distsq = 0, m_kinetic = 0, m_barrier = 0, m_lag = 0, m_total = 0;

    // Perform forward flow using the control u
    m_kinetic = ocsys->Flow(u);

    // Compute the goodness of fit metric using the terminal points
    MatrixRef q1(nvtx, 3, const_cast<double *>(ocsys->GetQt(param.nt-1).data_block()));
    MatrixRef d_g__d_qt1(nvtx, 3, d_g__d_qt[param.nt - 1].data_block());
    m_distsq = Traits::ComputeAugmentedLagrangianObjective(model, target, q1, d_g__d_qt1);

    // Compute the constraint violations
    for(int t = 0; t < param.nt; t++)
      {
      // Get the references to C and lambda for the current timepoint
      VectorRef Ct(nc_t, C.data_block() + nc_t * t);
      VectorRef lambda_t(nc_t, lambda.data_block() + nc_t * t);

      // Compute the constraints and gradients
      Traits::ComputeAugmentedLagrangianConstraints(model, ocsys->GetQt(t), d_g__d_qt[t], Ct, lambda_t, mu);
      }

    // Compute the contributions of the constraints to the AL objective
    m_barrier = C.squared_magnitude();
    m_lag = dot_product(C, lambda);

    // Compute total objective
    m_total = m_kinetic * param.w_kinetic 
      + m_distsq 
      + (mu / 2) * m_barrier
      - m_lag;

    if(f)
      *f = m_total;

    if(g)
      {
      if(verbose)
        {
        printf("Mu = %8.6f  |Lam| = %8.6f  DstSq = %8.6f  Kin = %8.6f  Bar = %8.6f  Lag = %8.6f  ETot = %12.8f\n",
          mu, lambda.inf_norm(), m_distsq, m_kinetic * param.w_kinetic, m_barrier * mu / 2, m_lag, m_total);
        }

      // Flow the gradient backward
      ocsys->FlowBackward(u, d_g__d_qt, param.w_kinetic, d_g__d_ut);

      // Compute the final gradient
      unsigned int ig = 0;
      for(int t = 0; t < param.nt; t++)
        for(int i = 0; i < q0.rows(); i++)
          for(int a = 0; a < 3; a++, ig++)
            (*g)[ig] = d_g__d_ut[t](i,a);
      }
    }

  // Update the lambdas
  void update_lambdas()
    {
    lambda -= mu * C;
    }

  void SetMu(double mu) { this->mu = mu; }

  double Mu() const { return mu; }

  void Export(const char *fn_pattern)
    {
    // Iterate over time
    for(int t = 0; t < param.nt; t++)
      {
      const Matrix &qt = ocsys->GetQt(t);
      VectorRef Ct(nc_t, C.data_block() + nc_t * t);
      VectorRef lambda_t(nc_t, lambda.data_block() + nc_t * t);

      char fn[4096];
      sprintf(fn, fn_pattern, t);
      Traits::ExportTimepoint(model, qt, Ct, lambda_t, fn);
      }
    }

protected:

  CMRep *model, *target;
  const AugLagMedialFitParameters &param;

  // Number of variables per time-step and total
  unsigned int nvtx, nvar_t, nvar_total;

  // Number of constraints per time-step and total
  unsigned int nc_t, nc_total;

  // The vector of lagrange multipliers
  CMRep::Vector lambda, u_init;

  // Inputs to the geodesic shooting
  Matrix q0, q1, p1;

  // Values of the constraints
  CMRep::Vector C;

  // The hamiltonian system
  OCSystem *ocsys;

  // Matrix array used to store the control
  OCSystem::MatrixArray u, d_g__d_qt, d_g__d_ut, d_ke__d_ut;

  // Mu
  double mu;

  // Verbosity
  bool verbose;
};




/**
 * This one does geodesic shooting, constraints at the end
 */
template <class Traits>
class PointMatchingWithEndpointConstraintsAugLagObjective : public vnl_cost_function
{
public:

  typedef vnl_matrix_ref<double> MatrixRef;
  typedef vnl_vector_ref<double> VectorRef;
  typedef CMRep::Matrix Matrix;
  typedef CMRep::Vector Vector;
  typedef PointSetHamiltonianSystem<double, 3> HSystem;

  PointMatchingWithEndpointConstraintsAugLagObjective(
    const AugLagMedialFitParameters &in_param, 
    CMRep *in_model, CMRep *in_target)
    : model(in_model), target(in_target), param(in_param), 
      nvtx(Traits::GetNumberOfActiveVertices(in_model)),
      nvar_t(Traits::GetNumberOfActiveVertices(in_model) * 3), 
      vnl_cost_function(Traits::GetNumberOfActiveVertices(in_model) * 3)
    {
    // Compute the initial p0-vector. 
    p0_init.set_size(nvar_t);
    Traits::ComputeInitialControl(model, target, p0_init, param.nt);

    // Initialize the initial landmarks
    q0.set_size(nvtx, 3);
    Traits::ComputeInitialLandmarks(model, q0);

    // Compute constant terms of Hessian
    Traits::PrecomputeHessianData(model, target, &hess_data);

    // Initialize the hamiltonian system
    hsys = new HSystem(q0, param.sigma, param.nt);

    // Initialize the p0 and various derivatives
    p0.set_size(nvtx, 3);
    q1.set_size(nvtx, 3);
    p1.set_size(nvtx, 3);
    d_g__d_p0.set_size(nvtx, 3);
    d_g__d_q1.set_size(nvtx, 3);

    // Initialize mu based on parameters
    mu = param.mu_init;

    // The number of constraints per timepoint and total
    nc_t = Traits::GetNumberOfConstraintsPerTimepoint(model);

    // Initialize the constraints matrices
    C.set_size(nc_t); C.fill(0.0);

    // Initialize the lambdas
    lambda.set_size(nc_t); lambda.fill(0.0);

    // Verbosity
    verbose = false;
    }

  /** Destructor */
  ~PointMatchingWithEndpointConstraintsAugLagObjective() 
    {
    delete hsys;
    }

  /** Get suitable x for initialization */
  Vector get_xinit() const { return p0_init; }

  /** Get number of variables */
  unsigned int get_nvar() const { return p0_init.size(); }

  void set_verbose(bool flag) { this->verbose = flag; }

  /** Compute the objective function */
  virtual void compute(const CMRep::Vector &x, double *f, CMRep::Vector *g)
    {
    // TODO: this copy-in is a waste of time and space
    p0.copy_in(x.data_block());

    // Initialize the array of objective function derivatives w.r.t. q_t
    d_g__d_q1.fill(0.0);

    // Initialize the various components of the objective function
    double m_distsq = 0, m_kinetic = 0, m_barrier = 0, m_lag = 0, m_total = 0;

    // Perform forward flow using the control u
    m_kinetic = hsys->FlowHamiltonian(p0, q1, p1);

    /*

    // Compute the goodness of fit metric using the terminal points
    m_distsq = Traits::ComputeAugmentedLagrangianObjective(model, target, q1, d_g__d_q1);

    // Compute the constraint violations
    Traits::ComputeAugmentedLagrangianConstraints(model, q1, d_g__d_q1, C, lambda, mu);

    // Compute the contributions of the constraints to the AL objective
    m_barrier = C.squared_magnitude();
    m_lag = dot_product(C, lambda);

    */

    // Use the matrix-based computation model
    VectorRef q1_flat(nvar_t, q1.data_block());
    VectorRef d_g__d_q1_flat(nvar_t, d_g__d_q1.data_block());

    // Compute the AL and gradient
    double AL = Traits::ComputeAugmentedLagrangianAndGradient(
      model, q1_flat, d_g__d_q1_flat, C, lambda, mu, &hess_data);

    // Compute the separate contributions of the constraints to the AL objective (for output)
    m_barrier = C.squared_magnitude();
    m_lag = dot_product(C, lambda);
    m_distsq = AL - ( (mu / 2) * m_barrier - m_lag );

    // Compute total objective with the kinetic energy
    m_total = AL + m_kinetic * param.w_kinetic;

    if(f)
      *f = m_total;

    if(g)
      {
      if(verbose)
        {
        printf("Mu = %8.6f  |Lam| = %8.6f  DstSq = %8.6f  Kin = %8.6f  Bar = %8.6f  Lag = %8.6f  ETot = %12.8f\n",
          mu, lambda.inf_norm(), m_distsq, m_kinetic * param.w_kinetic, m_barrier * mu / 2, m_lag, m_total);
        }

      // Flow the gradient backward
      d_g__d_p0.fill(0.0);

      Matrix dummy(nvtx, 3, 0.0);
      hsys->FlowGradientBackward(d_g__d_q1, dummy, d_g__d_p0);

      // Add in the kinetic energy gradient
      hsys->ComputeHamiltonianJet(q0, p0, false);
      for(unsigned int a = 0; a < 3; a++)
        for(unsigned int k = 0; k < nvtx; k++)
          d_g__d_p0(k, a) += param.w_kinetic * hsys->GetHp(a)[k];

      // Compute the final gradient
      unsigned int ig = 0;
      for(int i = 0; i < q0.rows(); i++)
        for(int a = 0; a < 3; a++, ig++)
          (*g)[ig] = d_g__d_p0(i,a);
      }
    }


  // Compute the Allassonniere transversality term G
  void ComputeTransversality(CMRep::Vector &x, Vector &G)
    {
    // TODO: this copy-in is a waste of time and space
    p0.copy_in(x.data_block());

    // Initialize the array of objective function derivatives w.r.t. q_t
    d_g__d_q1.fill(0.0);

    // Perform flow with Jacobian computations
    double m_kinetic = hsys->FlowHamiltonian(p0, q1, p1);

    // Compute the Hessian of the Augmented Lagrangian objective w.r.t. q1
    VectorRef q1_flat(nvar_t, q1.data_block());
    VectorRef d_g__d_q1_flat(nvar_t, d_g__d_q1.data_block());
    double AL = Traits::ComputeAugmentedLagrangianAndGradient(
      model, q1_flat, d_g__d_q1_flat, C, lambda, mu, &hess_data);

    // Multiply the Hessian by the Jacobian of q (a bit ugly)
    G.set_size(nvar_t);

    // Iterate over the rows of the AL Hessian
    for(unsigned int i = 0; i < nvar_t; i++)
      {
      // Get the column in Dq1 corresponding to j
      unsigned int vtx_i = i / 3, dim_i = i % 3;

      // Compute transversality measure G
      G[i] = p1(vtx_i, dim_i) * param.w_kinetic + d_g__d_q1_flat(i);
      }
    }

  // Compute the Allassonniere transversality terms G and DG 
  void ComputeTransversalityAndJacobian(CMRep::Vector &x, Vector &G, Matrix &DG)
    {
    // TODO: this copy-in is a waste of time and space
    p0.copy_in(x.data_block());

    // Initialize the array of objective function derivatives w.r.t. q_t
    d_g__d_q1.fill(0.0);

    // These matrices store the Jacobians of q1 and p1 wrt p0
    Matrix Dq1[3][3], Dp1[3][3];
    for(unsigned int a = 0; a < 3; a++)
      {
      for(unsigned int b = 0; b < 3; b++)
        {
        Dq1[a][b].set_size(nvtx,nvtx);
        Dp1[a][b].set_size(nvtx,nvtx);
        }
      }

    // Perform flow with Jacobian computations
    double m_kinetic = hsys->FlowHamiltonianWithGradient(p0, q1, p1, Dq1, Dp1);

    // Compute the Hessian of the Augmented Lagrangian objective w.r.t. q1
    VectorRef q1_flat(nvar_t, q1.data_block());
    VectorRef d_g__d_q1_flat(nvar_t, d_g__d_q1.data_block());
    double AL = Traits::ComputeAugmentedLagrangianJet(
      model, q1_flat, d_g__d_q1_flat, C, lambda, mu, &hess_data);

    // Multiply the Hessian by the Jacobian of q (a bit ugly)
    DG.set_size(nvar_t, nvar_t);
    G.set_size(nvar_t);
    ImmutableSparseMatrix<double> &HL = hess_data.HL;

    // Iterate over the rows of the AL Hessian
    for(unsigned int i = 0; i < nvar_t; i++)
      {
      // Get the column in Dq1 corresponding to j
      unsigned int vtx_i = i / 3, dim_i = i % 3;

      // Compute transversality measure G
      G[i] = p1(vtx_i, dim_i) * param.w_kinetic + d_g__d_q1_flat(i);

      for(unsigned int j = 0; j < nvar_t; j++)
        {
        // Get the column in Dq1 corresponding to j
        unsigned int vtx_j = j / 3, dim_j = j % 3;

        // We are computing DG(i, j)
        double DG_i_j = Dp1[dim_i][dim_j](vtx_i, vtx_j) * param.w_kinetic;

        for(ImmutableSparseMatrix<double>::RowIterator it = HL.Row(i); !it.IsAtEnd(); ++it)
          {
          // Get the row in Dq1 corresponding to k
          unsigned int k = it.Column();
          unsigned int vtx_k = k / 3, dim_k = k % 3;

          // Multiply and add
          DG_i_j += it.Value() * Dq1[dim_k][dim_j](vtx_k, vtx_j);
          }

        DG(i,j) = DG_i_j;
        }
      }

    // Print the current state
    if(verbose)
      {
      double m_barrier = C.squared_magnitude();
      double m_lag = dot_product(C, lambda);
      double m_distsq = AL - ( (mu / 2) * m_barrier - m_lag );
      double m_total = AL + m_kinetic * param.w_kinetic;

      printf("Mu = %8.6f  |Lam| = %8.6f  DstSq = %8.6f  Kin = %8.6f  Bar = %8.6f  Lag = %8.6f  ETot = %12.8f\n",
        mu, lambda.inf_norm(), m_distsq, m_kinetic * param.w_kinetic, m_barrier * mu / 2, m_lag, m_total);
      }
    }

  // Perform an Allassonniere iteration
  void iter_Allassonniere(CMRep::Vector &x)
    {
    // Transversality terms
    Matrix DG;
    Vector G;


    // Compute transversality at x
    this->ComputeTransversalityAndJacobian(x, G, DG);

    // Apply random perturbation to x

    // DERIVATIVE CHECK
    if(param.check_deriv)
      {
      vnl_random rndy;
      double eps = 1e-6;

      for(int k = 0; k < 16; k++)
        {
        int i = rndy.lrand32(0, x.size());
        Vector xtest = x;
        Vector G1, G2;
        xtest[i] = x[i] - eps;
        this->ComputeTransversality(xtest, G1);

        xtest[i] = x[i] + eps;
        this->ComputeTransversality(xtest, G2);

        Vector dG_i = (G2-G1)/(2.0*eps);

        printf("i = %04d,   |AG| = %12.8f,  |NG| = %12.8f,  |Del| = %12.8f\n", 
          i, DG.get_column(i).inf_norm(), dG_i.inf_norm(), (DG.get_column(i) - dG_i).inf_norm());
        }
      }

    // Perform singular value decomposition on the Hessian matrix, zeroing
    // out the singular values below 1.0 (TODO: use a relative threshold?)
    vnl_svd<double> svd(DG, -0.001);
    int nnz = 0;
    for(int i = 0; i < svd.W().rows(); i++)
      if(svd.W()(i,i) != 0.0)
        nnz++;

    printf("SVD min: %12.8f, max: %12.8f, nnz: %d, rank: %d\n", 
      svd.sigma_min(), svd.sigma_max(), nnz, svd.rank());

    // Compute inv(DG) * G
    Vector del_p0_vec = - svd.solve(G);

    // Perform Brent iterations
    std::cout << "Brent iterations" << std::endl;
    BrentObjective<PointMatchingWithEndpointConstraintsAugLagObjective> brent_obj(x, del_p0_vec, this);
    vnl_brent_minimizer brent(brent_obj);
    double alpha = brent.minimize(0.1);

    // Subtract from x
    x += alpha * del_p0_vec;
    }

  // Update the lambdas
  void update_lambdas()
    {
    lambda -= mu * C;
    }

  void SetMu(double mu) { this->mu = mu; }

  double Mu() const { return mu; }

  void Export(const char *fn_pattern)
    {
    // Iterate over time
    for(int t = 0; t < param.nt; t++)
      {
      const Matrix &qt = hsys->GetQt(t);

      char fn[4096];
      sprintf(fn, fn_pattern, t);
      Traits::ExportTimepoint(model, qt, C, lambda, fn);
      }
    }

protected:

  CMRep *model, *target;
  const AugLagMedialFitParameters &param;

  // Number of variables per time-step and total
  unsigned int nvtx, nvar_t;

  // Number of constraints per time-step and total
  unsigned int nc_t;

  // The vector of lagrange multipliers
  CMRep::Vector lambda, p0_init;

  // Inputs to the geodesic shooting
  Matrix q0, q1, p1;

  // Values of the constraints
  CMRep::Vector C;

  // The hamiltonian system
  HSystem *hsys;

  // Hessian data
  typename Traits::HessianData hess_data;

  // Matrix array used to store the control
  typename HSystem::Matrix p0, d_g__d_q1, d_g__d_p0;

  // Mu
  double mu;

  // Verbosity
  bool verbose;
};


template <class TObjective>
void DerivativeCheck(TObjective &obj, CMRep::Vector &x, int iter)
{
  obj.set_verbose(false);
  printf("******* ANALYTIC GRADIENT TEST (ITER %d) *******\n", iter);
  vnl_random rndy;

  double eps = 1e-6;
  typename CMRep::Vector test_grad(x.size());
  double f_test;
  obj.compute(x, &f_test, &test_grad);
  for(int k = 0; k < 16; k++)
    {
    int i = rndy.lrand32(0, x.size());
    typename CMRep::Vector xtest = x;
    double f1, f2;
    xtest[i] = x[i] - eps;
    obj.compute(xtest, &f1, NULL);

    xtest[i] = x[i] + eps;
    obj.compute(xtest, &f2, NULL);

    printf("i = %04d,   AG = %12.8f,  NG = %12.8f,  Del = %12.8f\n", 
      i, test_grad[i], (f2 - f1) / (2 * eps), 
      fabs(test_grad[i] - (f2-f1)/(2*eps)));
    }

  // Perform optimization
  obj.set_verbose(true);
}






int main(int argc, char *argv[])
{
  // Read the template and target meshes
  if(argc < 3)
    {
    std::cerr << "Usage: " << argv[0] << " template.vtk target.vtk" << std::endl;
    return -1;
    }

  // Read the template mesh
  CMRep m_template;
  m_template.ReadVTK(argv[1]);

  // Reat the target mesh
  CMRep m_target;
  m_target.ReadVTK(argv[2]);

  // Read the target mesh
  /*
  vtkSmartPointer<vtkPolyDataReader> reader = vtkPolyDataReader::New();
  reader->SetFileName(argv[2]);
  reader->Update();
  vtkSmartPointer<vtkPolyData> m_target = reader->GetOutput();

  // Put the coordinates into a matrix
  CMRep::Matrix x_target(m_target->GetNumberOfPoints(), 3);
  for(int i = 0; i < m_target->GetNumberOfPoints(); i++)
    for(int d = 0; d < 3; d++)
      x_target(i,d) = m_target->GetPoint(i)[d];
      */

  // Time to set up the objective function
  AugLagMedialFitParameters param;

  typedef PointBasedMediallyConstrainedFittingTraits TraitsType;
  if(param.interp_mode)
    {
    typedef PointMatchingWithTimeConstraintsAugLagObjective<TraitsType> ObjectiveType;
    ObjectiveType obj(param, &m_template, &m_target);
    obj.set_verbose(true);

    // The current iterate -- start with zero momentum
    CMRep::Vector x_opt = obj.get_xinit();

    // Actually initialize with the distance to the target points

    // Outer iteration loop
    for(int it = 0; it < 10; it++)
      {
      // Inner iteration loop
      if(param.check_deriv)
        DerivativeCheck(obj, x_opt, it);
      
      // Uncomment this code to test derivative computation
      // Create an optimization
      vnl_lbfgsb optimizer(obj);
      // vnl_conjugate_gradient optimizer(obj);
      optimizer.set_f_tolerance(1e-9);
      optimizer.set_x_tolerance(1e-4);
      optimizer.set_g_tolerance(1e-6);
      optimizer.set_trace(true);
      optimizer.set_max_function_evals(param.gradient_iter);

      // vnl_conjugate_gradient optimizer(cost_fn);
      optimizer.minimize(x_opt);

      // Update the lambdas
      obj.update_lambdas();

      // Update the mu
      obj.SetMu(obj.Mu() * param.mu_scale);

      // Export the meshes
      char fn_dir[4096], fn_pattern[4096];
      sprintf(fn_dir, "/tmp/testau_iter_%02d", it);
      sprintf(fn_pattern, "%s/testau_iter_%02d_tp_%s.vtk", fn_dir, it, "%03d");
      vtksys::SystemTools::MakeDirectory(fn_dir);
      obj.Export(fn_pattern);
      }
    }
  else
    {
    typedef PointMatchingWithEndpointConstraintsAugLagObjective<TraitsType> ObjectiveType;
    ObjectiveType obj(param, &m_template, &m_target);
    obj.set_verbose(true);

    // The current iterate -- start with zero momentum
    CMRep::Vector x_opt = obj.get_xinit();

    // Actually initialize with the distance to the target points

    // Outer iteration loop
    for(int it = 0; it < 10; it++)
      {
      // Inner iteration loop
      if(param.check_deriv)
        DerivativeCheck(obj, x_opt, it);
      
      // Uncomment this code to test derivative computation
      // Create an optimization
      vnl_lbfgsb optimizer(obj);
      // vnl_conjugate_gradient optimizer(obj);
      optimizer.set_f_tolerance(1e-9);
      optimizer.set_x_tolerance(1e-4);
      optimizer.set_g_tolerance(1e-6);
      optimizer.set_trace(true);
      optimizer.set_max_function_evals(param.gradient_iter);

      // vnl_conjugate_gradient optimizer(cost_fn);
      optimizer.minimize(x_opt);

      // Now do some iterations of Allassonniere
      for(int i = 0; i < param.newton_iter; i++)
        obj.iter_Allassonniere(x_opt);

      // Update the lambdas
      obj.update_lambdas();

      // Update the mu
      obj.SetMu(obj.Mu() * param.mu_scale);

      // Export the meshes
      char fn_dir[4096], fn_pattern[4096];
      sprintf(fn_dir, "/tmp/testau_iter_%02d", it);
      sprintf(fn_pattern, "%s/testau_iter_%02d_tp_%s.vtk", fn_dir, it, "%03d");
      vtksys::SystemTools::MakeDirectory(fn_dir);
      obj.Export(fn_pattern);
      }
    }


  // Read the two meshes
  // vtkSmartPointer<vtkPolyData> m_targ = ReadMesh(argv[2]);
}
