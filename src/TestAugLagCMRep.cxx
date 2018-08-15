/**
 * This is a test of the new cm-rep approach that combines geodesic shooting
 * and boundary-based medial constraints
 */

#include "MeshTraversal.h"
#include "MedialException.h"
#include "PointSetOptimalControlSystem.h"
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
  unsigned int mu_update_freq;

  // Interpolation mode or fitting mode?
  bool interp_mode;

  // Default initializer
  AugLagMedialFitParameters() 
    : nt(40), w_kinetic(0.05), sigma(4.0), 
      mu_init(0.0001), mu_scale(1.6), mu_update_freq(20),
      interp_mode(true)  {}
};


class AugLagMedialFitObjective : public vnl_cost_function
{
public:

  typedef vnl_matrix_ref<double> MatrixRef;
  typedef vnl_vector_ref<double> VectorRef;
  typedef CMRep::Matrix Matrix;
  typedef CMRep::Vector Vector;
  typedef PointSetOptimalControlSystem<double, 3> OCSystem;

  AugLagMedialFitObjective(
    const AugLagMedialFitParameters &in_param, 
    CMRep *in_model, CMRep *in_target)
    : model(in_model), target(in_target), param(in_param), 
      vnl_cost_function(in_param.nt * (in_model->nv + in_model->nmv) * 3)
    {
    // The constraint set:
    //   - each spoke is orthogonal to the boundary surface (2 * nv_b)
    //   - spokes are equal length (nm_bitangent + 2 * nm_tritangent)
    
    // Compute the number of equal length constraints
    nc_t = 2 * model->nv;
    for(int j = 0; j < model->med_bi.size(); j++)
      nc_t += model->med_bi[j].size() - 1;

    // Initialize the lambdas
    lambda.set_size(nc_t * param.nt);
    lambda.fill(0.0);

    // Generate the initial control vector. This vector is large, number of all 
    // moving points times the number of time steps.
    unsigned int nvar = param.nt * (model->nmv + model->nv) * 3;
    
    // Generate the initial momentum vector. The initial momentum consists of the 
    // momentum on the boundary vertices and the momentum on the medial triangle
    // centers. We initialize with a linear displacement towards the target point
    x_init.set_size(nvar);
    x_init.fill(0.0);

    // Initialize the hamiltonian system
    q0.set_size(model->nv + model->nmv, 3);
    q0.update(model->bnd_vtx, 0, 0);
    q0.update(model->med_vtx, model->nv, 0);
    ocsys = new OCSystem(q0, param.sigma, param.nt);

    // Initialize the u array
    u.resize(param.nt, Matrix(model->nv + model->nmv, 3));
    d_g__d_qt.resize(param.nt, Matrix(model->nv + model->nmv, 3));
    d_obj__d_ut.resize(param.nt, Matrix(model->nv + model->nmv, 3));

    // Initialize x_init to contain a constant velocity field (for now)
    for(int t = 0; t < param.nt; t++)
      {
      MatrixRef x_init_t(model->nv + model->nmv, 3, x_init.data_block() + t * q0.rows());
      x_init_t.update((target->bnd_vtx - model->bnd_vtx) / param.nt, 0, 0);
      x_init_t.update((target->med_vtx - model->med_vtx) / param.nt, model->nv, 0);
      }

    // Initialize mu based on parameters
    mu = param.mu_init;

    // Initialize the constraints matrices
    C.set_size(nc_t * param.nt); C.fill(0.0);

    // Verbosity
    verbose = false;
    }

  /** Destructor */
  ~AugLagMedialFitObjective() 
    {
    delete ocsys;
    }

  /** Get suitable x for initialization */
  Vector get_xinit() const { return x_init; }

  /** Get number of variables */
  unsigned int get_nvar() const { return x_init.size(); }

  void set_verbose(bool flag) { this->verbose = flag; }

  /** Compute the objective function */
  virtual void compute(const CMRep::Vector &x, double *f, CMRep::Vector *g)
    {
    // TODO: this copy-in is a waste of time and space
    for(int t = 0; t < param.nt; t++)
      u[t].copy_in(x.data_block() + (model->nv + model->nmv) * 3 * t);

    // Initialize the various components of the objective function
    double m_distsq = 0, m_kinetic = 0, m_barrier = 0, m_lag = 0, m_total = 0;

    // Initialize the array of objective function derivatives w.r.t. q_t
    for(int t = 0; t < param.nt; t++)
      d_g__d_qt[t].fill(0.0);

    // The reference to the last of these matrices (involved in the final objective)
    MatrixRef d_obj__d_qbnd1(model->nv, 3, d_g__d_qt[param.nt - 1].data_block());

    // Perform forward flow using the control u
    m_kinetic = ocsys->Flow(u);

    // Compute the goodness of fit metric
    MatrixRef q1_bnd(model->nv, 3, const_cast<double *>(ocsys->GetQt(param.nt-1).data_block()));
    for(int i = 0; i < model->nv; i++)
      {
      for(int d = 0; d < 3; d++)
        {
        double del_id = q1_bnd(i,d) - target->bnd_vtx(i,d);
        m_distsq += del_id * del_id;

        // Compute the contribution to d_obj__d_q1
        d_obj__d_qbnd1(i,d) += 2 * del_id; 
        }
      }

    // Compute the constraint violations
    int t_start = param.interp_mode ? 0 : param.nt - 1;
    for(int t = t_start; t < param.nt; t++)
      {
      // Get the q for the current timepoint
      const Matrix &qt = ocsys->GetQt(t);
      MatrixRef qt_bnd(model->nv, 3, const_cast<double *>(qt.data_block()));
      MatrixRef qt_med(model->nmv, 3, const_cast<double *>(qt.data_block()) + model->nv * 3);

      // Get references to the constraints in the list C
      MatrixRef Ct_orth(model->nv, 2, C.data_block() + nc_t * t);
      VectorRef Ct_spklen(nc_t - 2 * model->nv, C.data_block() + nc_t * t + 2 * model->nv);
      
      // Get the lambdas for the current timeslice
      MatrixRef Lt_orth(model->nv, 2, lambda.data_block() + nc_t * t);
      VectorRef Lt_spklen(nc_t - 2 * model->nv, lambda.data_block() + nc_t * t + 2 * model->nv);

      // Get the partial derivatives for the current timeslice
      MatrixRef d_obj__d_qbnd_t(model->nv, 3, d_g__d_qt[t].data_block());
      MatrixRef d_obj__d_qmed_t(model->nmv, 3, d_g__d_qt[t].data_block() + 3 * model->nv);

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
          double d_obj__d_C = mu * Ct_orth(j, d) - Lt_orth(j, d);

          // Contribution to the ring of vertices
          for(CMRep::SparseRowIter it = model->wgt_Quv[d].Row(j); !it.IsAtEnd(); ++it)
            {
            double z = d_obj__d_C * it.Value();
            for(int a = 0; a < 3; a++)
              d_obj__d_qbnd_t(it.Column(), a) += z * S[a];
            }

          // Contribution to the medial and boundary vertices
          for(int a = 0; a < 3; a++)
            {
            double z = d_obj__d_C * Qd[a];
            d_obj__d_qbnd_t(j, a) += z;
            d_obj__d_qmed_t(model->bnd_mi[j], a) -= z;
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
            double d_obj__d_C = mu * Ct_spklen[m] - Lt_spklen[m];
            double z = 2.0 * d_obj__d_C;

            for(int a = 0; a < 3; a++)
              {
              // Contribution to boundary nodes
              d_obj__d_qbnd_t(bj, a) += z * delj[a];
              d_obj__d_qbnd_t(b0, a) -= z * del0[a];
              d_obj__d_qmed_t(k, a) += z * (del0[a] - delj[a]);
              }
            }
          }
        }
      } // loop over time 

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
      ocsys->FlowBackward(u, d_g__d_qt, param.w_kinetic, d_obj__d_ut);

      // Compute the final gradient
      unsigned int ig = 0;
      for(int t = 0; t < param.nt; t++)
        for(int i = 0; i < q0.rows(); i++)
          for(int a = 0; a < 3; a++, ig++)
            (*g)[ig] = d_obj__d_ut[t](i,a);
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
      MatrixRef qt_bnd(model->nv, 3, const_cast<double *>(qt.data_block()));
      MatrixRef qt_med(model->nmv, 3, const_cast<double *>(qt.data_block()) + model->nv * 3);

      // Create a copy of the mesh in the model
      vtkSmartPointer<vtkPolyData> pd = vtkPolyData::New();
      pd->DeepCopy(model->bnd_vtk);

      // Update the points
      for(int i = 0; i < qt_bnd.rows(); i++)
        pd->GetPoints()->SetPoint(i, qt_bnd(i,0), qt_bnd(i,1), qt_bnd(i,2));

      // Get the constraint values for the current timeslice
      MatrixRef Ct_orth(model->nv, 2, C.data_block() + nc_t * t);
      VectorRef Ct_spklen(nc_t - 2 * model->nv, C.data_block() + nc_t * t + 2 * model->nv);

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
      char fn[4096];
      sprintf(fn, fn_pattern, t);
      vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
      writer->SetFileName(fn);
      writer->SetInputData(pd);
      writer->Update();
      }

    }

protected:

  CMRep *model, *target;
  const AugLagMedialFitParameters &param;

  // Number of constraints per time-step
  unsigned int nc_t;

  // The vector of lagrange multipliers
  CMRep::Vector lambda, x_init;

  // Inputs to the geodesic shooting
  Matrix q0, q1, p1;

  // Values of the constraints
  CMRep::Vector C;

  // The hamiltonian system
  OCSystem *ocsys;

  // Matrix array used to store the control
  OCSystem::MatrixArray u, d_g__d_qt, d_obj__d_ut, d_ke__d_ut;

  // Mu
  double mu;

  // Verbosity
  bool verbose;
};


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
  AugLagMedialFitObjective obj(param, &m_template, &m_target);

  // The current iterate -- start with zero momentum
  CMRep::Vector x_opt = obj.get_xinit();

  // Actually initialize with the distance to the target points

  // Outer iteration loop
  for(int it = 0; it < 10; it++)
    {
    // Inner iteration loop
    
    // Uncomment this code to test derivative computation
    obj.set_verbose(false);
    printf("******* ANALYTIC GRADIENT TEST (ITER %d) *******\n", it);
    vnl_random rndy;

    double eps = 1e-6;
    typename CMRep::Vector test_grad(x_opt.size());
    double f_test;
    obj.compute(x_opt, &f_test, &test_grad);
    for(int k = 0; k < 16; k++)
      {
      int i = rndy.lrand32(0, x_opt.size());
      typename CMRep::Vector xtest = x_opt;
      double f1, f2;
      xtest[i] = x_opt[i] - eps;
      obj.compute(xtest, &f1, NULL);

      xtest[i] = x_opt[i] + eps;
      obj.compute(xtest, &f2, NULL);

      int nv = (m_template.nv + m_template.nmv);
      int i_t = i / (nv * 3);
      int i_v = (i - i_t * (nv * 3)) / 3;
      int i_d = i % 3;
      printf("t = %02d, i = %04d, d = %d    AG = %12.8f,  NG = %12.8f,  Del = %12.8f\n", 
        i_t, i_v, i_d, test_grad[i], (f2 - f1) / (2 * eps), 
        fabs(test_grad[i] - (f2-f1)/(2*eps)));
      }

    // Perform optimization
    obj.set_verbose(true);

    // Create an optimization
    vnl_lbfgsb optimizer(obj);
    // vnl_conjugate_gradient optimizer(obj);
    optimizer.set_f_tolerance(1e-9);
    optimizer.set_x_tolerance(1e-4);
    optimizer.set_g_tolerance(1e-6);
    optimizer.set_trace(true);
    optimizer.set_max_function_evals(200);

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


  // Read the two meshes
  // vtkSmartPointer<vtkPolyData> m_targ = ReadMesh(argv[2]);
}
