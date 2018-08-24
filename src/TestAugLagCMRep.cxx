/**
 * This is a test of the new cm-rep approach that combines geodesic shooting
 * and boundary-based medial constraints
 */

#include "MeshTraversal.h"
#include "MedialException.h"
#include "PointSetOptimalControlSystem.h"
#include "PointSetHamiltonianSystem.h"
#include "SparseMatrix.h"
#include "FastLinearInterpolator.h"
#include "MedialAtomGrid.h"
#include "SubdivisionSurface.h"

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

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>

#include <vnl_matrix.h>
#include <vnl_vector.h>
#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_lbfgs.h>
#include <vnl/algo/vnl_lbfgsb.h>
#include <vnl/algo/vnl_conjugate_gradient.h>
#include <vnl/vnl_random.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_brent_minimizer.h>


#include <nlopt.h>

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
  Matrix bnd_vtx, med_vtx, bnd_nrm;

  // Radius function
  Vector med_R;

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
  this->bnd_nrm.set_size(nv, 3);
  vtkDataArray *da_mi = this->bnd_vtk->GetPointData()->GetArray("MedialIndex");
  for(int k = 0; k < nv; k++)
    {
    this->bnd_mi[k] = da_mi->GetTuple1(k);
    for(int d = 0; d < 3; d++)
      {
      this->bnd_vtx(k, d) = this->bnd_vtk->GetPoint(k)[d];
      this->bnd_nrm(k, d) = b_norm->GetComponent(k,d);
      }
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
  this->med_R.set_size(nmv);

  for(unsigned int i = 0; i < nmv; i++)
    {
    // Get the first boundary point
    unsigned int j = this->med_bi[i][0];

    // Compute the medial point coordinate
    this->med_R[i] = b_rad->GetTuple1(j);

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
    : nt(40), w_kinetic(0.05), sigma(2.0), 
      mu_init(1), mu_scale(1.0), gradient_iter(1600), newton_iter(0),
      interp_mode(true), check_deriv(true)  {}
};


struct QuadraticForm 
{
  typedef vnl_vector<double> Vector;
  typedef ImmutableSparseMatrix<double> SparseMat;

  SparseMat A;
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
  typedef vnl_vector<double> Vector;
  typedef ImmutableSparseMatrix<double> SparseMat;

  SparseMat M;
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

// Compute triange area and derivative
double TriangleAreaAndGradient(
  const vnl_vector<double> &A, const vnl_vector<double> &B, const vnl_vector<double> &C,
  vnl_vector<double> &d_a__d_A, vnl_vector<double> &d_a__d_B, vnl_vector<double> &d_a__d_C)
{
  vnl_vector<double> N = 0.25 * vnl_cross_3d(B - A, C - A);

  d_a__d_A = vnl_cross_3d(N, C - B);
  d_a__d_B = vnl_cross_3d(N, A - C);
  d_a__d_C = vnl_cross_3d(N, B - A);

  return N.magnitude() * 2.0;
}



template <class TFunction>
class MeshFunctionVolumeIntegral
{
public:
  typedef vnl_matrix<double> Matrix;
  typedef vnl_vector<double> Vector;

  MeshFunctionVolumeIntegral(CMRep *model, unsigned int n_layers, unsigned int sub_level, TFunction &f)
    : func(f)
    {
    this->model = model;
    this->n_layers = n_layers;

    // Set up the parent level mesh
    *((TriangleMesh *)(&B)) = model->bnd_tri;
    B.SetAsRoot();

    // Set up the child mesh
    SubdivisionSurface::RecursiveSubdivide(&B, &S, sub_level, true);

    vtkSmartPointer<vtkPolyData> test = vtkPolyData::New();
    SubdivisionSurface::ApplySubdivision(model->bnd_vtk, test, S);
    vtkSmartPointer<vtkPolyDataWriter> w = vtkPolyDataWriter::New();
    w->SetInputData(test);
    w->SetFileName("/tmp/subflat.vtk");
    w->Update();



    // Allocate the array of samples
    samples.resize((n_layers + 1) * S.nVertices);

    // Allocate the subdivided boundary and medial vertices
    qb_sub.set_size(S.nVertices, 3);
    qm_sub.set_size(S.nVertices, 3);
    d_obj__d_qb_sub.set_size(S.nVertices, 3);
    d_obj__d_qm_sub.set_size(S.nVertices, 3);

    // Initialize the triangles of the wedges. The quadrilateral faces must be divided into
    // triangles and this is done in such a way that adjacent wedges have the same triangles,
    // so that there is no overlap between wedge triangulations and no holes. 
    for(unsigned int k = 0; k < n_layers; k++)
      {
      for(auto &tri : S.triangles)
        {
        // Create a wedge
        Wedge W;

        // Lambda to get the index of a sample
        auto si = [&](unsigned int l, unsigned int t) { return S.nVertices * l + t; };

        // Indices of the six vertices of this wedge
        unsigned int A = W.V[0] = si(k  , tri.vertices[0]);
        unsigned int B = W.V[1] = si(k  , tri.vertices[1]);
        unsigned int C = W.V[2] = si(k  , tri.vertices[2]);
        unsigned int D = W.V[3] = si(k+1, tri.vertices[0]);
        unsigned int E = W.V[4] = si(k+1, tri.vertices[1]);
        unsigned int F = W.V[5] = si(k+1, tri.vertices[2]);

        // Lambda to set a row in a wedge
        auto set_row = [&](unsigned int row, unsigned int v1, unsigned int v2, unsigned int v3) 
          { W.T(row,0) = v1; W.T(row,1) = v2, W.T(row,2) = v3; };

        // The triangles at the top and bottom of the wedge
        set_row(0, A, B, C);
        set_row(1, D, F, E);

        // The remaining triangles - here we use the ordering of the samples to get consistent 
        // splitting of quad faces
        if(A < C)
          { set_row(2, A, C, F); set_row(3, A, F, D); }
        else
          { set_row(2, C, F, D); set_row(3, C, D, A); }

        if(C < B)
          { set_row(4, C, B, E); set_row(5, C, E, F); }
        else
          { set_row(4, B, E, F); set_row(5, B, F, C); }

        if(B < A)
          { set_row(6, B, A, D); set_row(7, B, D, E); }
        else
          { set_row(6, A, D, E); set_row(7, A, E, B); }

        // Append the wedge
        wedges.push_back(W);
        }
      }
    }

  void Compute(const Matrix &qb, const Matrix &qm, double &v, double &fv)
    {
    // Apply subdivision to outer meshes qb and qm
    qb_sub.fill(0.0); qm_sub.fill(0.0);
    for(unsigned int i = 0; i < S.nVertices; i++)
      {
      for(auto it = S.weights.Row(i); !it.IsAtEnd(); ++it)
        {
        double w = it.Value();
        unsigned int j = it.Column(), j_med = model->bnd_mi[j];

        for(unsigned int a = 0; a < 3; a++)
          {
          qb_sub(i,a) += qb(j, a) * w;
          qm_sub(i,a) += qm(j_med, a) * w;
          }
        }
      }

    // Compute the samples and sample the function at each location
    for(unsigned int k = 0, is = 0; k <= n_layers; k++)
      {
      double l = k * 1.0 / n_layers;
      for(unsigned int i = 0; i < S.nVertices; i++)
        {
        // Get the current sample
        Sample &s = samples[is++];

        // Compute the sample's location
        for(unsigned int a = 0; a < 3; a++)
          s.x[a] = qb_sub(i,a) * (1-l) + qm_sub(i,a) * l;

        // Compute the function and gradient
        s.f = func.Compute(s.x, s.grad_f);

        // Set the volume element of the sample to zero
        s.vol_elt = 0.0;
        }
      }

    // Iterate over the wedges to compute volume
    for(auto &W : wedges)
      {
      // Volume is computed using triple scalar products (divergence formula)
      W.vol = 0.0;
      for(unsigned int t = 0; t < 8; t++)
        {
        const Vec3 &P = samples[W.T[t][0]].x;
        const Vec3 &Q = samples[W.T[t][1]].x;
        const Vec3 &R = samples[W.T[t][2]].x;
        W.vol += ScalarTripleProduct(P, Q - P, R - P) / 6.0;
        }

      // Assign an equal portion of the volume to each of the six vertices
      for(unsigned int v = 0; v < 6; v++)
        samples[W.V[v]].vol_elt += W.vol / 6.0;
      }

    // Integrate the volume and the function value
    v = 0.0; fv = 0.0;
    for(auto &s : samples)
      {
      v += s.vol_elt;
      fv += s.f * s.vol_elt;
      }
    }

  void ComputeAdjoint(
    const Matrix &qb, const Matrix &qm, 
    double d_obj__d_v, double d_obj__d_fv,
    Matrix &d_obj__d_qb, Matrix &d_obj__d_qm)
    {
    for(auto &s : samples)
      {
      // Compute the derivative of the objective wrt each sample's volume element
      s.d_obj__d_vol_elt = d_obj__d_v + s.f * d_obj__d_fv;

      // Initialize the derivative of the objective with respect to the sample
      s.d_obj__d_x = s.vol_elt * d_obj__d_fv * s.grad_f;
      }

    for(auto &W : wedges)
      {
      // For each wedge, compute the derivative of the objective wrt its volume
      double d_obj__d_vol = 0;
      for(unsigned int v = 0; v < 6; v++)
        d_obj__d_vol += samples[W.V[v]].d_obj__d_vol_elt;
      d_obj__d_vol /= 6.0;

      // Now compute the derivative of the volume with respect to each vertex
      for(unsigned int t = 0; t < 8; t++)
        {
        const Vec3 &P = samples[W.T[t][0]].x;
        const Vec3 &Q = samples[W.T[t][1]].x;
        const Vec3 &R = samples[W.T[t][2]].x;
        Vec3 &dP = samples[W.T[t][0]].d_obj__d_x;
        Vec3 &dQ = samples[W.T[t][1]].d_obj__d_x;
        Vec3 &dR = samples[W.T[t][2]].d_obj__d_x;

        // Compute the adjoint of the scalar triple product, which just consists
        // of the cross-products of pairs of the vectors
        dP += vnl_cross_3d(Q, R) * d_obj__d_vol / 6.0;
        dQ += vnl_cross_3d(R, P) * d_obj__d_vol / 6.0;
        dR += vnl_cross_3d(P, Q) * d_obj__d_vol / 6.0;
        }
      }

    // Now compute the derivative of the objective with respect to the subdivided 
    // vertex coordinates
    d_obj__d_qb_sub.fill(0.0); d_obj__d_qm_sub.fill(0.0);

    for(unsigned int k = 0, is = 0; k <= n_layers; k++)
      {
      double l = k * 1.0 / n_layers;
      for(unsigned int i = 0; i < S.nVertices; i++)
        {
        // Get the current sample
        Sample &s = samples[is++];

        // Map its contribution to the subdivision surfaces
        for(unsigned int a = 0; a < 3; a++)
          {
          d_obj__d_qb_sub(i,a) += (1 - l) * s.d_obj__d_x[a];
          d_obj__d_qm_sub(i,a) += l * s.d_obj__d_x[a];
          }
        }
      }

    // Finally, map the contributions back to the original vertices
    d_obj__d_qb.fill(0.0); d_obj__d_qm.fill(0.0);

    for(unsigned int i = 0; i < S.nVertices; i++)
      {
      for(auto it = S.weights.Row(i); !it.IsAtEnd(); ++it)
        {
        double w = it.Value();
        unsigned int j = it.Column(), j_med = model->bnd_mi[j];

        for(unsigned int a = 0; a < 3; a++)
          {
          d_obj__d_qb(j,a)      += d_obj__d_qb_sub(i,a) * w;
          d_obj__d_qm(j_med, a) += d_obj__d_qm_sub(i,a) * w;
          }
        }
      }
    }

protected:
  typedef SubdivisionSurface::MeshLevel MeshLevel;
  typedef vnl_vector_fixed<double, 3> Vec3;

  // The model
  CMRep *model;

  // Number of layers
  unsigned int n_layers;

  // Base layers (input mesh resolution) and subdivided layers
  MeshLevel B, S;

  // Function being integrated
  TFunction &func;

  // Subdivided boundary layer and medial layer
  Matrix qb_sub, qm_sub, d_obj__d_qb_sub, d_obj__d_qm_sub;

  // Array of samples from the objective function
  struct Sample
    {
    Vec3 x, grad_f, d_obj__d_x;
    double f;
    double vol_elt, d_obj__d_vol_elt;
    };

  // A flat list of samples (stored in layer-major order)
  std::vector<Sample> samples;

  // A wedge
  struct Wedge 
    {
    // A wedge constitutes eight triangles, each one references a sample
    vnl_matrix_fixed<unsigned int, 8, 3> T;

    // Also store the vertices of the wedge, each one references a sample
    unsigned int V[6];

    // Volume of the wedge
    double vol;

    // Partial derivatives of the wedge w.r.t. six vertices
    Vec3 d_V__d_x[2][3];
    };

  // A flat list of wedges (stored in layer-major order)
  std::vector<Wedge> wedges;
};


/**
 * A function for testing Dice computations on exact data
 */
struct TestFunction
{
  typedef vnl_vector_fixed<double, 3> Vec3;
  double Compute(const Vec3 &X, Vec3 &grad)
    {
    double x = X[0], y = X[1], z = X[2];
    double f = cos(x * y - z) - sin(y * z - x);
    grad[0] = -sin(x * y - z) * y    - cos(y * z - x) * -1.0;
    grad[1] = -sin(x * y - z) * x    - cos(y * z - x) * z;
    grad[2] = -sin(x * y - z) * -1.0 - cos(y * z - x) * y;
    return f;
    }
};

template <class TFunction>
class DiceOverlapComputation
{
public:
  typedef vnl_matrix<double> Matrix;
  typedef vnl_vector<double> Vector;

  DiceOverlapComputation(
    CMRep *model, unsigned int n_wedges, unsigned int sub_level, 
    double vol_image, TFunction &in_func)
    : integrator(model, n_wedges, sub_level, in_func)
    {
    this->vol_image = vol_image;
    }

  double Compute(const Matrix &qb, const Matrix &qm)
    {
    integrator.Compute(qb, qm, v, fv);
    return 2.0 * fv / (v + vol_image);
    }

  void ComputeAdjoint(
    const Matrix &qb, const Matrix &qm,
    double d_obj__d_dice,
    Matrix &d_obj__d_qb, Matrix &d_obj__d_qm)
    {
    double d_obj__d_fv = d_obj__d_dice * 2.0 / (v + vol_image);
    double d_obj__d_v = - d_obj__d_dice * 2.0 * fv / ((v + vol_image) * (v + vol_image));
    integrator.ComputeAdjoint(qb, qm, d_obj__d_v, d_obj__d_fv, d_obj__d_qb, d_obj__d_qm);
    }

private:
  typedef MeshFunctionVolumeIntegral<TFunction> Integrator;
  Integrator integrator;

  double v, fv, vol_image;
};


class ImageDiceFunction 
{
public:
  typedef itk::Image<double, 3> ImageType;
  typedef FastLinearInterpolator<ImageType, double, 3> InterpolatorType;

  ImageDiceFunction(const char *fname, double sigma)
    {
    // Read input image
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(fname);
    reader->Update();
    m_InputImage = reader->GetOutput();

    // Apply Gaussian smoothing
    typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilter;
    typename GaussianFilter::Pointer fltSmooth = GaussianFilter::New();
    fltSmooth->SetInput(m_InputImage);
    fltSmooth->SetSigma(sigma);
    fltSmooth->Update();
    m_SmoothedImage = fltSmooth->GetOutput();

    // Compute the volume of the image
    m_Volume = 0.0;
    for(size_t i = 0; i < m_SmoothedImage->GetPixelContainer()->Size(); i++)
      m_Volume += m_SmoothedImage->GetBufferPointer()[i];
    for(size_t d = 0; d < 3; d++)
      m_Volume *= m_SmoothedImage->GetSpacing()[d];

    // Create interpolator
    m_Interp = new InterpolatorType(m_SmoothedImage);

    // Compute the transformation from physical space to voxel space
    ComputeVoxelSpaceToNiftiSpaceTransform(m_SmoothedImage);
    }

  ~ImageDiceFunction()
    {
    delete m_Interp;
    }

  typedef vnl_vector_fixed<double, 3> Vec3;

  double Compute(const Vec3 &X, Vec3 &grad)
    {
    double f = 0.0;
    double *g = grad.data_block();

    // Map x to voxel coordinates
    Vec3 x_vox = m_A_RAS_to_IJK * X + m_b_RAS_to_IJK;
    Vec3 grad_vox(0.0);
    double *grad_vox_p = grad_vox.data_block();

    // Interpolate
    InterpolatorType::InOut rc = 
      m_Interp->InterpolateWithGradient(x_vox.data_block(), &f, &grad_vox_p);

    // Map gradient back to physical space
    grad = m_A_RAS_to_IJK * grad_vox;

    if(rc == InterpolatorType::OUTSIDE)
      grad.fill(0.0);

    return f; 
    }

  double GetVolume() const { return m_Volume; }

  // Helper function to map from ITK coordiante space to RAS space
  void ComputeVoxelSpaceToNiftiSpaceTransform(ImageType *image)
    {
    // Generate intermediate terms
    vnl_matrix<double> m_dir, m_ras_matrix;
    vnl_diag_matrix<double> m_scale, m_lps_to_ras;
    vnl_vector<double> v_origin, v_ras_offset;

    // Compute the matrix
    m_dir = image->GetDirection().GetVnlMatrix();
    m_scale.set(image->GetSpacing().GetVnlVector());
    m_lps_to_ras.set(vnl_vector<double>(3, 1.0));
    m_lps_to_ras[0] = -1;
    m_lps_to_ras[1] = -1;
    m_A_IJK_to_RAS = m_lps_to_ras * m_dir * m_scale;

    // Compute the vector
    v_origin = image->GetOrigin().GetVnlVector();
    m_b_IJK_to_RAS = m_lps_to_ras * v_origin;

    // Invert this transform
    m_A_RAS_to_IJK = vnl_matrix_inverse<double>(m_A_IJK_to_RAS);
    m_b_RAS_to_IJK = - m_A_RAS_to_IJK * m_b_IJK_to_RAS;
    }

protected:
  ImageType::Pointer m_InputImage, m_SmoothedImage;
  InterpolatorType *m_Interp;
  double m_Volume;

  vnl_matrix_fixed<double, 3, 3> m_A_RAS_to_IJK, m_A_IJK_to_RAS;
  vnl_vector_fixed<double, 3> m_b_RAS_to_IJK, m_b_IJK_to_RAS;
};


template <class TFunction>
void TestDice(CMRep *model, TFunction &tf, double img_vol, double eps = 1.0e-6)
{
  typedef vnl_matrix<double> Matrix;
  typedef vnl_vector<double> Vector;
  typedef vnl_matrix_ref<double> MatrixRef;
  typedef vnl_vector_ref<double> VectorRef;

  // Storage for variables
  unsigned int k = model->bnd_vtx.size() + model->med_vtx.size();
  Vector y(k, 0.0), d_y(k, 0.0);

  double *py = y.data_block();
  MatrixRef qb(model->nv, 3, py);      py += qb.size();
  MatrixRef qm(model->nmv, 3, py); 
  
  double *pdy = d_y.data_block();
  MatrixRef d_qb(model->nv, 3, pdy);   pdy += qb.size();
  MatrixRef d_qm(model->nmv, 3, pdy);  
  
  // Assign variables
  qb.update(model->bnd_vtx);
  qm.update(model->med_vtx);

  // Dice test
  MeshFunctionVolumeIntegral<TFunction> dicer(model, 5, 2, tf);

  // Values returned by the Dice test
  double fv, v;

  // Lambda for performing dice test
  auto dice_test = [&](double v, double fv) { return 2 * fv / (v + img_vol); };
  // auto dice_test = [&](double v, double fv) { return v + fv; };

  // Compute function and gradient of fv
  double t0 = clock();
  dicer.Compute(qb, qm, v, fv);
  printf("Dice = %12.4f computed in %12.8f ms \n", dice_test(v, fv), (clock() - t0) * 1000 / CLOCKS_PER_SEC);

  // To compute adjoint, we need the adjoint of the Dice function
  double d_dice__d_fv = 2.0 / (v + img_vol);
  double d_dice__d_v  = - 2.0 * fv / ((v + img_vol) * (v + img_vol));
  dicer.ComputeAdjoint(qb,qm,d_dice__d_v,d_dice__d_fv,d_qb,d_qm);
  printf("Dice adjoint computed in %12.8f ms \n", (clock() - t0) * 1000 / CLOCKS_PER_SEC);

  // Numeric differentiation
  vnl_random rndy;
  for(int i = 0; i < 50; i++)
    {
    int j = rndy.lrand32(0, y.size());
    double y_init = y[j];

    y[j] = y_init + eps;
    dicer.Compute(qb,qm,v,fv);
    double f1 = dice_test(v, fv);

    y[j] = y_init - eps;
    dicer.Compute(qb,qm,v,fv);
    double f2 = dice_test(v, fv);

    y[j] = y_init;

    double d_f__d_yj_num = (f1 - f2) / (2 * eps);

    printf("Dice Test: Var: %4d  An: %16.12f  Nu: %16.12f  Del: %16.12f\n", 
      j, d_y[j], d_f__d_yj_num, fabs(d_y[j] - d_f__d_yj_num));
    }
}


/**
 * Specific implementation of constraints with slack variables N and R
 */
template <class TImageFunction>
class PointBasedMediallyConstrainedFittingTraits
{
public:
  typedef vnl_vector<double> Vector;
  typedef vnl_matrix<double> Matrix;
  typedef vnl_vector<unsigned int> IdxVector;
  typedef vnl_matrix<unsigned int> IdxMatrix;
  typedef vnl_matrix_ref<double> MatrixRef;
  typedef vnl_vector_ref<double> VectorRef;

  // Target data represents terms that attract the model
  struct TargetData
    {
    CMRep *target_model;
    TImageFunction *image_function;
    double w_target_model, w_image_function;
    };

  // Components of the variable vector (X)
  struct XComponents
    {
    MatrixRef u_bnd, u_med, N_bnd;
    VectorRef R_med;

    XComponents(CMRep *model, double *p) :
      u_bnd(model->nv,  3,  p),
      u_med(model->nmv, 3,  p + model->nv * 3),
      N_bnd(model->nv,  3,  p + model->nv * 3 + model->nmv * 3),
      R_med(model->nmv,     p + model->nv * 6 + model->nmv * 3) { }
    };

  // Components of the variable + Q vector (Y)
  struct YComponents : public XComponents
    {
    MatrixRef q_bnd, q_med;

    YComponents(CMRep *model, double *p) :
      XComponents(model, p),
      q_bnd(model->nv,  3,  p + model->nv * 6 + model->nmv * 4),
      q_med(model->nmv, 3,  p + model->nv * 9 + model->nmv * 4) {}
    };

  // Get the number of vertices controling diffeomorphic transformation
  static int GetNumberOfActiveVertices(CMRep *model)
    {
    return (model->nv + model->nmv);
    }

  // Get the number of slack variables at a given timepoint. These are variables
  // involved in optimization but not active vertices in geodesic shooting
  static int GetNumberOfSlackVariables(CMRep *model)
    {
    // For each boundary vertex, we have N (3-vector) and for each medial vertex we have R
    return 3 * model->nv + model->nmv;
    }

  // Among all optimization variables, copy the momenta of active vertices into
  // the vector p.
  static void GetActiveVertexMomenta(CMRep *model, const Vector &x, unsigned int t, Vector &p)
    {
    unsigned int nvar_t = GetNumberOfActiveVertices(model) * 3 + GetNumberOfSlackVariables(model);
    p.copy_in(x.data_block() + nvar_t * t);
    }

  static int GetNumberOfConstraintsPerTimepoint(CMRep *model)
    {
    // The constraint set:
    //   - normal orthogonal to the boundary surface (2 * nv_b)
    //   - normal has unit length (nv_b)
    //   - M + normal * R = B (3 * nv_b)
    return 6 * model->nv;
    }

  // Initialize the optimization variable for a timepoint. The input vector x contains the 
  // variables for just the current timepoint
  static void ComputeInitialization(CMRep *model, const TargetData *target, Vector &x, unsigned int nt)
    {
    // Get a pointer to the storage for timepoint t
    XComponents xcmp(model, x.data_block());

    // The initial momentum
    x.fill(0.0);
    // TODO: restore this?
    // xcmp.u_bnd.update((target->target_model->bnd_vtx - model->bnd_vtx) / nt, 0, 0);
    // xcmp.u_med.update((target->target_model->med_vtx - model->med_vtx) / nt, 0, 0);

    // Normals and radii are set to be that of initial timepoint
    for(int i = 0; i < model->nv; i++)
      {
      // N * R is just normalized vector from medial to boundary
      Vector Ni = model->bnd_vtx.get_row(i) - model->med_vtx.get_row(model->bnd_mi[i]);
      xcmp.R_med[model->bnd_mi[i]] = Ni.magnitude();
      Ni.normalize();
      xcmp.N_bnd.set_row(i, Ni);
      }
    }

  static void ComputeInitialLandmarks(CMRep *model, Matrix &q0)
    {
    q0.update(model->bnd_vtx, 0, 0);
    q0.update(model->med_vtx, model->nv, 0);
    }

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


  /** 
   * Compute the Augmented Lagrangian and its derivative using precomputed quantities. 
   * The input vector Y and the output vector d_AL__d_y refer to the concatenation of
   * the optimization variables x and the geodesic landmarks qt for the current time
   * point. The same applies to the Hessian - it is a matrix containing second derivatives
   * in both x and qt
   */
  static double ComputeAugmentedLagrangianAndGradient(
    CMRep *model, TargetData *target, const Vector &Y, Vector &d_AL__d_Y,
    Vector &C, Vector &lambda, double mu, HessianData *H,
    unsigned int t, unsigned int nt)
    {
    // Break the input and output vectors into components
    YComponents Ycmp(model, const_cast<double *>(Y.data_block()));
    YComponents d_Ycmp(model, d_AL__d_Y.data_block());

    // The objective function is only computed at the last timepoint
    double AL = 0.0;
    if(t == nt - 1)
      {
      // Initialize the outputs with the f/grad of the function f
      AL += H->qf_F.Compute(Y, d_AL__d_Y);

      // Compute the Dice overlap
      AL += target->w_image_function * (1.0 - target->image_function->Compute(Ycmp.q_bnd, Ycmp.q_med));

      // Compute the adjoint of the Dice overlap
      target->image_function->ComputeAdjoint(
        Ycmp.q_bnd, Ycmp.q_med, -target->w_image_function, d_Ycmp.q_bnd, d_Ycmp.q_med);
      }

    // Vector containing the gradient of the constraint - to be reused
    Vector grad_Cj(Y.size(), 0.0);

    // Iterate over the constraints
    for(int j = 0; j < H->qf_C.size(); j++)
      {
      // Reference to the current quadratic form
      QuadraticForm &Z = H->qf_C[j];

      // Compute the constraint and its gradient
      C[j] = Z.Compute(Y, grad_Cj);

      // Add to the overall objective
      AL += C[j] * (C[j] * mu / 2 - lambda[j]);

      // Add to the overall gradient
      d_AL__d_Y += grad_Cj * (mu * C[j] - lambda[j]);
      }

    return AL;
    }


  static IdxMatrix MakeIndexMatrix(unsigned int rows, unsigned int cols, unsigned int start_idx)
    {
    IdxMatrix M(rows, cols);
    for(int i = 0; i < M.size(); i++)
      M.data_block()[i] = start_idx + i;
    return M;
    }

  static IdxVector MakeIndexVector(unsigned int rows, unsigned int start_idx)
    {
    IdxVector M(rows);
    for(int i = 0; i < M.size(); i++)
      M.data_block()[i] = start_idx + i;
    return M;
    }

  /** 
   * This function precomputes different terms used to compute the Hessian that 
   * do not change between iterations (exploiting wide use of quadratic functions
   * in the objective
   */
  static void PrecomputeHessianData(CMRep *model, const TargetData *target, HessianData *data)
    {
    // Create index arrays into the input variable Y (includes x's and q's)
    // at the end, k holds total number of input variables
    unsigned int k = 0;
    IdxMatrix i_ub = MakeIndexMatrix(model->nv, 3, k);    k += i_ub.size();
    IdxMatrix i_um = MakeIndexMatrix(model->nmv, 3, k);   k += i_um.size();
    IdxMatrix i_Nb = MakeIndexMatrix(model->nv, 3, k);    k += i_Nb.size();
    IdxVector i_Rm = MakeIndexVector(model->nmv, k);      k += i_Rm.size();
    IdxMatrix i_qb = MakeIndexMatrix(model->nv, 3, k);    k += i_qb.size();
    IdxMatrix i_qm = MakeIndexMatrix(model->nmv, 3, k);   k += i_qm.size();

    // Create indices into the constraints (normal constraints and spoke constraints)
    unsigned int nc = 0;
    IdxMatrix ic_N   = MakeIndexMatrix(model->nv, 3, nc);  nc += ic_N.size();
    IdxMatrix ic_spk = MakeIndexMatrix(model->nv, 3, nc);  nc += ic_spk.size();

    // Initialize the constraint quadratic data
    data->qf_C.resize(nc);

    // Handle the per-boundary point constraints
    for(int j = 0; j < model->nv; j++)
      {
      // Constraint that the normal is orthogonal to the boundary tangent vectors
      // N . Qu = 0, i.e., sum_(j in nhd[i]) w_j N . Q[j] = 0
      for(int d = 0; d < 2; d++)
        {
        // Initialize the Hessian for this constraint
        vnl_sparse_matrix<double> H_Cj(k, k);

        for(CMRep::SparseRowIter it = model->wgt_Quv[d].Row(j); !it.IsAtEnd(); ++it)
          {
          for(int a = 0; a < 3; a++)
            {
            // Index of the neighbor's Q
            unsigned int i_qb_mov = i_qb(it.Column(), a);

            // Index of the boundary normal vector
            unsigned int i_Nb_j = i_Nb(j, a);
            
            // Set the Hessian terms
            H_Cj(i_qb_mov, i_Nb_j) += it.Value(); 
            H_Cj(i_Nb_j, i_qb_mov) += it.Value(); 
            }
          }

        // Set the A matrix, the b and c are zeros
        data->qf_C[ic_N(j, d)].Initialize(H_Cj, Vector(k, 0.0), 0.0);
        }

      // Constraint that the 'normal has unit length' constraint
      vnl_sparse_matrix<double> H_Cj(k, k);
      for(int a = 0; a < 3; a++)
        H_Cj(i_Nb(j,a), i_Nb(j,a)) = 2.0;
      data->qf_C[ic_N(j, 2)].Initialize(H_Cj, Vector(k, 0.0), -1.0);

      // Constraint on the spokes (Qb - N * R - Qm = 0)
      for(int a = 0; a < 3; a++)
        {
        // Indices of the involves quantities
        unsigned int i_qb_j = i_qb(j, a);
        unsigned int i_qm_j = i_qm(model->bnd_mi[j], a);
        unsigned int i_Nb_j = i_Nb(j, a);
        unsigned int i_Rm_j  = i_Rm(model->bnd_mi[j]);

        // The Hessian here only contains terms from N R
        vnl_sparse_matrix<double> H_Cspk(k, k);
        H_Cspk(i_Nb_j, i_Rm_j) = -1;
        H_Cspk(i_Rm_j, i_Nb_j) = -1;

        // There is also an offset vector with +1 in the qb row and -1 in the qm row
        Vector B_Cspk(k, 0.0);
        B_Cspk[i_qb_j] = 1.0;
        B_Cspk[i_qm_j] = -1.0;

        data->qf_C[ic_spk(j, a)].Initialize(H_Cspk, B_Cspk, 0.0);
        }
      }

    // Handle the objective function - it's a symmetric quadratic form in q_b's 
    vnl_sparse_matrix<double> M_f(k, k);
    Vector d_f(k, 0.0);

    double w_target_model_sqrt = sqrt(target->w_target_model);
    for(int j = 0; j < model->nv; j++)
      {
      for(int a = 0; a < 3; a++)
        {
        unsigned int i_qb_j = i_qb(j, a);
        M_f(i_qb_j,i_qb_j) = w_target_model_sqrt;
        d_f[i_qb_j] = - target->target_model->bnd_vtx(j, a) * w_target_model_sqrt;
        }
      }

    data->qf_F.Initialize(M_f, d_f);

    // TODO: initialize the Hessian entries for the Dice objective


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

  typedef std::vector<std::pair<std::string, double> > ConstraintDetail;

  static void UpdateConstraintDetail(ConstraintDetail &con_info, const std::string label, double value)
    {
    for(auto p : con_info)
      {
      if(p.first == label)
        {
        p.second = std::max(value, p.second);
        return;
        }
      }
    con_info.push_back(std::make_pair(label, value));
    }

  /**
   * Get details on different types of constraints and their maximum absolute value
   */
  static void GetConstraintDetails(CMRep *model, const Vector &C, ConstraintDetail &con_info)
    {
    double *pc = const_cast<double *>(C.data_block());
    MatrixRef Ct_N(model->nv, 3, pc);   pc += Ct_N.size();
    MatrixRef Ct_Spk(model->nv, 3, pc); pc += Ct_Spk.size();

    UpdateConstraintDetail(con_info, "C_NrmOrth", Ct_N.get_column(0).inf_norm());
    UpdateConstraintDetail(con_info, "C_NrmOrth", Ct_N.get_column(1).inf_norm());
    UpdateConstraintDetail(con_info, "C_NrmUnit", Ct_N.get_column(2).inf_norm());
    UpdateConstraintDetail(con_info, "C_Spk", Ct_Spk.absolute_value_max());
    }


  static void ExportTimepoint(CMRep *model, const Vector &Y, 
    const Vector &C, const Vector &lambda, const char *fname)
    {
    // Get the references to all the data in vector Y
    YComponents ycmp(model, const_cast<double *>(Y.data_block()));

    // Create a copy of the mesh in the model
    vtkSmartPointer<vtkPolyData> pd = vtkPolyData::New();
    pd->DeepCopy(model->bnd_vtk);

    // Create an array for normals
    vtkSmartPointer<vtkFloatArray> arr_nrm = vtkFloatArray::New();
    arr_nrm->SetNumberOfComponents(3);
    arr_nrm->SetNumberOfTuples(model->nv);

    // Create an array for radius
    vtkSmartPointer<vtkFloatArray> arr_rad = vtkFloatArray::New();
    arr_rad->SetNumberOfComponents(1);
    arr_rad->SetNumberOfTuples(model->nv);
    arr_rad->SetName("Radius");

    // Update the points
    for(int i = 0; i < ycmp.q_bnd.rows(); i++)
      {
      pd->GetPoints()->SetPoint(i, ycmp.q_bnd(i,0), ycmp.q_bnd(i,1), ycmp.q_bnd(i,2));
      arr_nrm->SetTuple3(i, ycmp.N_bnd(i,0), ycmp.N_bnd(i,1), ycmp.N_bnd(i,2));
      arr_rad->SetTuple1(i, ycmp.R_med(model->bnd_mi[i]));
      }

    pd->GetPointData()->SetNormals(arr_nrm);
    pd->GetPointData()->AddArray(arr_rad);

    // Get the constraint values for the current timeslice
    double *pc = const_cast<double *>(C.data_block());
    MatrixRef Ct_N(model->nv, 3, pc);   pc += Ct_N.size();
    MatrixRef Ct_Spk(model->nv, 3, pc); pc += Ct_Spk.size();

    // Create a point array for the constraints
    vtkSmartPointer<vtkFloatArray> arr_con = vtkFloatArray::New();
    arr_con->SetNumberOfComponents(6);
    arr_con->SetNumberOfTuples(model->nv);

    // Assign the ortho constraints
    for(int j = 0; j < model->nv; j++)
      {
      arr_con->SetComponent(j, 0, Ct_N(j, 0));
      arr_con->SetComponent(j, 1, Ct_N(j, 1));
      arr_con->SetComponent(j, 2, Ct_N(j, 2));
      arr_con->SetComponent(j, 3, Ct_Spk(j, 0));
      arr_con->SetComponent(j, 4, Ct_Spk(j, 1));
      arr_con->SetComponent(j, 5, Ct_Spk(j, 2));
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
class PointMatchingWithTimeConstraintsAugLagObjective
{
public:

  typedef vnl_matrix_ref<double> MatrixRef;
  typedef vnl_vector_ref<double> VectorRef;
  typedef CMRep::Matrix Matrix;
  typedef CMRep::Vector Vector;
  typedef PointSetOptimalControlSystem<double, 3> OCSystem;

  typedef typename Traits::TargetData TargetData;

  PointMatchingWithTimeConstraintsAugLagObjective(
    const AugLagMedialFitParameters &in_param, CMRep *in_model, TargetData *in_target)
    : model(in_model), target(in_target), param(in_param)
    {
    // Initialize problem size - include slack variables
    nvtx = Traits::GetNumberOfActiveVertices(model);
    nvar_t = nvtx * 3 + Traits::GetNumberOfSlackVariables(model);
    nvar_total = nvar_t * in_param.nt;

    printf("Objective with time constraints\n");
    printf("  nv = %d, nmv = %d\n", model->nv, model->nmv);
    printf("  nvar_t = %d, nvar_total = %d\n", nvar_t, nvar_total);

    // Compute the initial u-vector. 
    x_init.set_size(nvar_total);
    for(int t = 0; t < param.nt; t++)
      {
      VectorRef x_init_t(nvar_t, x_init.data_block() + nvar_t * t);
      Traits::ComputeInitialization(model, target, x_init_t, param.nt);
      }

    // Initialize the Y vector, which will hold the input vars and the q1's
    Y.set_size(nvar_t + nvtx * 3);
    d_AL__d_Y.set_size(Y.size());

    // Initialize the initial landmarks
    q0.set_size(nvtx, 3);
    Traits::ComputeInitialLandmarks(model, q0);

    // Compute constant terms of Hessian
    Traits::PrecomputeHessianData(model, target, &hess_data);

    // Initialize the hamiltonian system
    ocsys = new OCSystem(q0, param.sigma, param.nt);

    // There is a control (and constraints) at every time point
    u.resize(param.nt, Matrix(nvtx, 3));
    d_g__d_qt.resize(param.nt, Matrix(nvtx, 3));
    d_g__d_ut.resize(param.nt, Matrix(nvtx, 3));

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

    // Iteration counter
    iter_count = 0;
    }

  /** Destructor */
  ~PointMatchingWithTimeConstraintsAugLagObjective() 
    {
    delete ocsys;
    }

  /** Get suitable x for initialization */
  Vector get_xinit() const { return x_init; }

  /** Get the current constraint values */
  const Vector &GetC() const { return C; }

  /** Reset the evaluation counter */
  void ResetCounter() { this->iter_count = 0; }

  /** Get number of variables */
  unsigned int get_nvar() const { return x_init.size(); }

  void set_verbose(bool flag) { this->verbose = flag; }

  /** Print the problem status */
  void iter_print(unsigned int iter, 
    double m_distsq, double m_kinetic, 
    double m_barrier, double m_lag, double m_total)
    {
    // Get the details of the constraints
    typename Traits::ConstraintDetail con_detail;
    Traits::GetConstraintDetails(model, C, con_detail);

    // Combine into a string
    std::string con_text;
    char con_buffer[256]; 
    for(auto ci : con_detail)
      {
      sprintf(con_buffer, " |%s| = %8.4f ", ci.first.c_str(), ci.second);
      con_text += con_buffer;
      }

    printf(
      "Iter %05d  "
      "Mu = %8.4f  |Lam| = %8.4f  DstSq = %8.4f  Kin = %8.4f  "
      "Bar = %8.4f  Lag = %8.4f %s ETot = %12.8f\n",
      iter,
      mu, lambda.inf_norm(), m_distsq, m_kinetic * param.w_kinetic, m_barrier * mu / 2, m_lag, 
      con_text.c_str(), m_total);
    }

  /** Compute the objective function */
  virtual void compute(const CMRep::Vector &x, double *f, CMRep::Vector *g)
    {
    // TODO: this copy-in is a waste of time and space
    for(int t = 0; t < param.nt; t++)
      {
      u[t].copy_in(x.data_block() + nvar_t * t);
      d_g__d_ut[t].fill(0.0);
      }

    // Initialize the various components of the objective function
    double m_distsq = 0, m_kinetic = 0, m_barrier = 0, m_lag = 0, m_total = 0;

    // Perform forward flow using the control u
    m_kinetic = ocsys->Flow(u);

    // Initialize the augmented lagrangian to zero
    double AL = 0.0;

    // References into different parts of Y and d_AL__d_Y
    VectorRef Y_xt(nvar_t,  Y.data_block());
    MatrixRef Y_qt(nvtx, 3, Y.data_block() + nvar_t);
    MatrixRef d_AL__d_Y_qt(nvtx, 3, d_AL__d_Y.data_block() + nvar_t);

    // Compute the constraint violations
    for(int t = 0; t < param.nt; t++)
      {
      // Place the u coefficients and the slack variables into vector Y
      Y_xt.copy_in(x.data_block() + nvar_t * t);

      // Place the coordinates qt at the end of Y
      Y_qt.update(ocsys->GetQt(t));

      // Get the references to C and lambda for the current timepoint
      VectorRef Ct(nc_t, C.data_block() + nc_t * t);
      VectorRef lambda_t(nc_t, lambda.data_block() + nc_t * t);

      // Compute the objective, constraints and gradients for the current timepoint
      d_AL__d_Y.fill(0.0);
      AL += Traits::ComputeAugmentedLagrangianAndGradient(
        model, target, Y, d_AL__d_Y, Ct, lambda_t, mu, &hess_data, t, param.nt);

      // Copy the derivatives with respect to ut and slack variables into output gradient
      if(g)
        {
        // Copy the derivatives with respect to qt into the appropriate array
        d_g__d_qt[t].update(d_AL__d_Y_qt);

        // Copy the derivatives wrt u and slack to output gradient
        VectorRef gt(nvar_t, g->data_block() + nvar_t * t);
        gt.copy_in(d_AL__d_Y.data_block());
        }
      }

    // Compute the contributions of the constraints to the AL objective
    m_barrier = C.squared_magnitude();
    m_lag = dot_product(C, lambda);
    m_distsq = AL - ( (mu / 2) * m_barrier - m_lag );

    // Compute total objective with the kinetic energy
    m_total = AL + m_kinetic * param.w_kinetic;

    if(f)
      *f = m_total;

    if(g)
      {
      // Print status
      if(verbose)
        iter_print(iter_count++, m_distsq, m_kinetic, m_barrier, m_lag, m_total);

      // Flow the gradient backward
      ocsys->FlowBackward(u, d_g__d_qt, param.w_kinetic, d_g__d_ut);

      // Update the gradient based on the backward flow
      for(int t = 0; t < param.nt; t++)
        {
        VectorRef g_ut(nvtx * 3, g->data_block() + nvar_t * t);
        VectorRef del_g_ut(nvtx * 3, d_g__d_ut[t].data_block());
        g_ut += del_g_ut;
        }
      }
    }

  // Update the lambdas
  void update_lambdas()
    {
    lambda -= mu * C;
    }

  void SetMu(double mu) { this->mu = mu; }

  double Mu() const { return mu; }

  void Export(const CMRep::Vector &x, const char *fn_pattern)
    {
    // Iterate over time
    for(int t = 0; t < param.nt; t++)
      {
      // Place the u coefficients and the slack variables into vector Y
      VectorRef xt(nvar_t, const_cast<double *>(x.data_block()) + nvar_t * t);
      Y.update(xt);

      // Place the coordinates qt at the end of Y
      VectorRef qt_flat(nvtx * 3, const_cast<double *>(ocsys->GetQt(t).data_block()));
      Y.update(qt_flat, nvar_t);

      // Get the references to C and lambda for the current timepoint
      VectorRef Ct(nc_t, C.data_block() + nc_t * t);
      VectorRef lambda_t(nc_t, lambda.data_block() + nc_t * t);

      char fn[4096];
      sprintf(fn, fn_pattern, t);
      Traits::ExportTimepoint(model, Y, Ct, lambda_t, fn);
      }
    }

protected:

  // The cm-rep template
  CMRep *model;

  // The data fed into the objective function
  TargetData *target;

  // Parameters of optimization
  const AugLagMedialFitParameters &param;

  // Number of variables per time-step and total
  unsigned int nvtx, nvar_t, nvar_total;

  // Number of constraints per time-step and total
  unsigned int nc_t, nc_total;

  // The vector of lagrange multipliers
  CMRep::Vector C, lambda, x_init;

  // Inputs to the geodesic shooting
  Matrix q0;

  // Vector Y holds ut, slack variables, and qt all in a single storage
  Vector Y, d_AL__d_Y;

  // The hamiltonian system
  OCSystem *ocsys;

  // Hessian data
  typename Traits::HessianData hess_data;

  // Matrix array used to store the control
  OCSystem::MatrixArray u, d_g__d_qt, d_g__d_ut;

  // Mu
  double mu;

  // Verbosity
  bool verbose;
 
  // Iteration counter
  unsigned int iter_count;
};


/**
 * This one does geodesic shooting, constraints at the end
 */
template <class Traits>
class PointMatchingWithEndpointConstraintsAugLagObjective 
{
public:

  typedef vnl_matrix_ref<double> MatrixRef;
  typedef vnl_vector_ref<double> VectorRef;
  typedef CMRep::Matrix Matrix;
  typedef CMRep::Vector Vector;
  typedef PointSetHamiltonianSystem<double, 3> HSystem;

  typedef typename Traits::TargetData TargetData;

  PointMatchingWithEndpointConstraintsAugLagObjective(
    const AugLagMedialFitParameters &in_param, CMRep *in_model, TargetData *in_target)
    : model(in_model), target(in_target), param(in_param)
    {
    // Set the problem size variables
    nvtx = Traits::GetNumberOfActiveVertices(model);
    nvar_t = nvtx * 3 + Traits::GetNumberOfSlackVariables(model);

    // Compute the initial variables
    x_init.set_size(nvar_t);
    Traits::ComputeInitialization(model, target, x_init, param.nt);

    // Initialize the Y vector, which will hold the input vars and the q1's
    Y.set_size(nvar_t + nvtx * 3);
    d_AL__d_Y.set_size(Y.size());

    // Initialize the initial landmarks
    q0.set_size(nvtx, 3);
    Traits::ComputeInitialLandmarks(model, q0);

    // Compute constant terms of Hessian
    Traits::PrecomputeHessianData(model, target, &hess_data);

    // Initialize the hamiltonian system
    hsys = new HSystem(q0, param.sigma, param.nt);

    // Initialize the p0 and various derivatives
    p1.set_size(nvtx, 3);

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

    // Iteration counter
    iter_count = 0;
    }

  /** Destructor */
  ~PointMatchingWithEndpointConstraintsAugLagObjective() 
    {
    delete hsys;
    }

  /** Get suitable x for initialization */
  Vector get_xinit() const { return x_init; }

  /** Get the current constraint values */
  const Vector &GetC() const { return C; }

  /** Reset the evaluation counter */
  void ResetCounter() { this->iter_count = 0; }

  /** Get number of variables */
  unsigned int get_nvar() const { return x_init.size(); }

  void set_verbose(bool flag) { this->verbose = flag; }

  /** Print the problem status */
  void iter_print(unsigned int iter, 
    double m_distsq, double m_kinetic, 
    double m_barrier, double m_lag, double m_total)
    {
    // Get the details of the constraints
    typename Traits::ConstraintDetail con_detail;
    Traits::GetConstraintDetails(model, C, con_detail);

    // Combine into a string
    std::string con_text;
    char con_buffer[256]; 
    for(auto ci : con_detail)
      {
      sprintf(con_buffer, " |%s| = %8.4f ", ci.first.c_str(), ci.second);
      con_text += con_buffer;
      }

    printf(
      "Iter %05d  "
      "Mu = %8.4f  |Lam| = %8.4f  DstSq = %8.4f  Kin = %8.4f  "
      "Bar = %8.4f  Lag = %8.4f %s ETot = %12.8f\n",
      iter,
      mu, lambda.inf_norm(), m_distsq, m_kinetic * param.w_kinetic, m_barrier * mu / 2, m_lag, 
      con_text.c_str(), m_total);
    }


  /** Compute the objective function */
  virtual void compute(const CMRep::Vector &x, double *f, CMRep::Vector *g)
    {
    // The first part of Y contains the optimization variables x, second part has q's
    Y.update(x, 0);

    // Extract the initial momenta and final landmark positions from Y
    MatrixRef p0(nvtx, 3, Y.data_block());
    MatrixRef q1(nvtx, 3, Y.data_block() + nvar_t);

    // Initialize the various components of the objective function
    double m_distsq = 0, m_kinetic = 0, m_barrier = 0, m_lag = 0, m_total = 0;

    // Perform forward flow using the control u
    m_kinetic = hsys->FlowHamiltonian(p0, q1, p1);

    // Set the derivative to zeros
    d_AL__d_Y.fill(0.0);

    // Compute the AL and gradient
    double AL = Traits::ComputeAugmentedLagrangianAndGradient(
      model, target, Y, d_AL__d_Y, C, lambda, mu, &hess_data, param.nt - 1, param.nt);


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
      // Print status
      if(verbose)
        iter_print(iter_count++, m_distsq, m_kinetic, m_barrier, m_lag, m_total);

      // Before
      // std::cout << "d_AL__d_Y[0]: " << std::endl << d_AL__d_Y << std::endl;

      // Expose the pointer to the gradient of AL with respect to p0
      MatrixRef d_AL__d_p0(nvtx, 3, d_AL__d_Y.data_block());

      // Expose the pointer to the gradient of AL with respect to q1
      MatrixRef d_AL__d_q1(nvtx, 3, d_AL__d_Y.data_block() + nvar_t);

      // The beta is a dummy zero vector (because AL does not depend on p1)
      Matrix d_AL__d_p1(nvtx, 3, 0.0);
      hsys->FlowGradientBackward(d_AL__d_q1, d_AL__d_p1, d_AL__d_p0);

      // std::cout << "d_AL__d_Y[1]: " << std::endl << d_AL__d_Y << std::endl;

      // Add in the kinetic energy gradient
      hsys->ComputeHamiltonianJet(q0, p0, false);
      for(unsigned int a = 0; a < 3; a++)
        for(unsigned int k = 0; k < nvtx; k++)
          d_AL__d_p0(k, a) += param.w_kinetic * hsys->GetHp(a)[k];

      // std::cout << "d_AL__d_Y[2]: " << std::endl << d_AL__d_Y << std::endl;

      // Place the final gradient into g
      g->copy_in(d_AL__d_Y.data_block());
      }
    }

/*
 * DON'T KNOW HOW TO DO THIS WITH SLACK VARIABLES...
 *
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
    vnl_svd<double> svd(DG, -0.0001);
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

  */

  // Update the lambdas
  void update_lambdas()
    {
    lambda -= mu * C;
    }

  void SetMu(double mu) { this->mu = mu; }

  double Mu() const { return mu; }

  void Export(const CMRep::Vector &x, const char *fn_pattern)
    {
    Vector Yt = Y;
    MatrixRef qt(nvtx, 3, Yt.data_block() + nvar_t);

    // Iterate over time
    for(int t = 0; t < param.nt; t++)
      {
      // Assign the t'th landmarks to end of Y
      qt.update(hsys->GetQt(t));

      // Create a mesh
      char fn[4096];
      sprintf(fn, fn_pattern, t);
      Traits::ExportTimepoint(model, Yt, C, lambda, fn);
      }
    }

protected:

  // The cm-rep template
  CMRep *model;

  // The data fed into the objective function
  TargetData *target;

  // Parameters of optimization
  const AugLagMedialFitParameters &param;

  // Number of variables per time-step and total
  unsigned int nvtx, nvar_t;

  // Number of constraints per time-step and total
  unsigned int nc_t;

  // Iteration counter
  unsigned int iter_count;

  // The vector of lagrange multipliers
  CMRep::Vector lambda, x_init;

  // Values of the constraints
  CMRep::Vector C;

  // The hamiltonian system
  HSystem *hsys;

  // Hessian data
  typename Traits::HessianData hess_data;

  // Vector Y holds p0, slack variables, and q1 all in a single storage
  Vector Y, d_AL__d_Y;

  // Initial landmarks, endpoint momenta
  Matrix q0, p1;

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


/** ------------- NLOPT SUPPORT -------------------*/

/**
 * Wrapper around a cost function
 */
template <class TInnerFunction>
class vnl_func_wrapper : public vnl_cost_function
{
public:
  vnl_func_wrapper(TInnerFunction &in_func) 
    : func(in_func), vnl_cost_function(in_func.get_xinit().size()) {}

  virtual void compute(const CMRep::Vector &x, double *f, CMRep::Vector *g)
    {
    func.compute(x,f,g);
    }

private:
  TInnerFunction &func;
};

/**
 * This function wraps around a VNL objective function for interfacing with NLOPT
 */
double nlopt_vnl_func(unsigned n, const double *x, double *grad, void *my_func_data)
{
  vnl_cost_function *vnl_cf = static_cast<vnl_cost_function *>(my_func_data);
  vnl_vector_ref<double> x_vec(n, const_cast<double *>(x));
  double f = 0.0;
  
  if(grad)
    {
    vnl_vector_ref<double> grad_vec(n, grad);
    vnl_cf->compute(x_vec, &f, &grad_vec);
    }
  else
    {
    vnl_cf->compute(x_vec, &f, NULL);
    }

  return f;
}


template <class TObjective>
void optimize_auglag(
  TObjective &obj, const AugLagMedialFitParameters &param, double start_mu)
{
  // The current iterate -- start with zero momentum
  CMRep::Vector x_opt = obj.get_xinit();

  // Compute the initial value of mu using the Birgin and Martinez heuristics
  double f_current;
  obj.compute(x_opt, &f_current, NULL);
  double ssq_con = obj.GetC().squared_magnitude();
  double max_con = obj.GetC().inf_norm();
  double mu = std::max(1e-6, std::min(10.0, 2 * fabs(f_current) / ssq_con));
  double ICM = 1e100;

  // Override initial mu...
  mu = start_mu;

  // Set the initial mu
  obj.SetMu(mu);

  // Outer iteration loop
  for(int it = 0; it < 10; it++)
    {
    // Inner iteration loop
    if(param.check_deriv)
      DerivativeCheck(obj, x_opt, it);

    // Few iterations of CGD (because NLOPT goes nuts)
    vnl_func_wrapper<TObjective> obj_vnl(obj);
    vnl_conjugate_gradient optimizer(obj_vnl);
    optimizer.set_f_tolerance(1e-9);
    optimizer.set_x_tolerance(1e-4);
    optimizer.set_g_tolerance(1e-6);
    optimizer.set_trace(true);
    optimizer.set_max_function_evals(5);
    optimizer.minimize(x_opt);

    // Perform the inner optimization
    nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS, x_opt.size());
    nlopt_set_min_objective(opt, nlopt_vnl_func, &obj_vnl);
    nlopt_set_xtol_rel(opt, 1e-5);
    nlopt_set_ftol_rel(opt, 1e-5);
    nlopt_set_maxeval(opt, param.gradient_iter);
    double f_opt;
    int rc = nlopt_optimize(opt, x_opt.data_block(), &f_opt);
    switch(rc)
      {
    case NLOPT_SUCCESS: printf("NLOPT: Success!\n"); break;
    case NLOPT_STOPVAL_REACHED: printf("NLOPT: Reached f_stopval!\n"); break;
    case NLOPT_FTOL_REACHED: printf("NLOPT: Reached f_tol!\n"); break;
    case NLOPT_XTOL_REACHED: printf("NLOPT: Reached x_tol!\n"); break;
    case NLOPT_MAXEVAL_REACHED: printf("NLOPT: Reached max evaluations!\n"); break;
    default: printf("nlopt failed %d!\n",rc);
      }
    nlopt_destroy(opt);

    printf("*** End of inner iteration loop %d ***\n", it);

    // Update the lambdas
    obj.update_lambdas();

    // Update the mu
    obj.compute(x_opt, &f_current, NULL);
    double newICM = obj.GetC().inf_norm();
    printf("Constraint one-norm [before] : %12.4f  [after]: %12.4f\n", ICM, newICM);
    if(newICM > 0.5 * ICM)
      mu *= 10.0;
    obj.SetMu(mu);
    ICM = newICM;

    // Export the meshes
    char fn_dir[4096], fn_pattern[4096];
    sprintf(fn_dir, "/tmp/testau_iter_%02d", it);
    sprintf(fn_pattern, "%s/testau_iter_%02d_tp_%s.vtk", fn_dir, it, "%03d");
    vtksys::SystemTools::MakeDirectory(fn_dir);
    obj.Export(x_opt, fn_pattern);
    }
}


int main(int argc, char *argv[])
{
  // Read the template and target meshes
  if(argc < 3)
    {
    std::cerr << "Usage: " << argv[0] << " template.vtk target.nii.gz target.vtk" << std::endl;
    return -1;
    }

  // Read the template mesh
  CMRep m_template;
  m_template.ReadVTK(argv[1]);

  // Read the target image file
  ImageDiceFunction idf(argv[2], 0.2);

  // Reat the target mesh
  CMRep m_target;
  m_target.ReadVTK(argv[3]);

  // The starting mu
  double start_mu = atof(argv[4]);

  // Time to set up the objective function
  AugLagMedialFitParameters param;

  if(param.check_deriv)
    {
    printf("Testing Dice derivatives with analytic function\n");
    TestFunction test_fun;
    TestDice(&m_template, test_fun, 1.0);

    printf("Testing Dice derivatives with actual image\n");
    TestDice(&m_template, idf, idf.GetVolume(), 1e-4);
    }

  // Initialize the Dice overlap computer
  typedef DiceOverlapComputation<ImageDiceFunction> DiceOverlapType;
  DiceOverlapType dicer(&m_template, 5, 2, idf.GetVolume(), idf);

  // Define the traits for the AL objective
  typedef PointBasedMediallyConstrainedFittingTraits<DiceOverlapType> TraitsType;

  // Set up the target data
  TraitsType::TargetData td;
  td.w_target_model = 0.0;
  td.w_image_function = 100.0;
  td.target_model = &m_target;
  td.image_function = &dicer;

  if(param.interp_mode)
    {
    // Set up the objective function
    typedef PointMatchingWithTimeConstraintsAugLagObjective<TraitsType> ObjectiveType;
    ObjectiveType obj(param, &m_template, &td);
    obj.set_verbose(true);

    // Do the optimization
    optimize_auglag(obj, param, start_mu);
    }
  else
    {
    typedef PointMatchingWithEndpointConstraintsAugLagObjective<TraitsType> ObjectiveType;
    ObjectiveType obj(param, &m_template, &td);
    obj.set_verbose(true);

    // Do the optimization
    optimize_auglag(obj, param, start_mu);
    }
}
