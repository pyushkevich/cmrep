#include "PointSetHamiltonianSystem.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/vnl_cost_function.h"
#include "vnl/algo/vnl_lbfgs.h"

#include <iostream>
#include <algorithm>
#include "util/ReadWriteVTK.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkQuadricClustering.h"
#include "vtkSmartPointer.h"
#include "vtkCell.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"

#include "CommandLineHelper.h"
#include "GreedyAPI.h"

using namespace std;

int usage()
{
  cout << "lmshoot: Geodesic shooting for landmarks" << endl;
  cout << "Usage:" << endl;
  cout << "  lmshoot [options]" << endl;
  cout << "Required Options:" << endl;
  cout << "  -m template.vtk target.vtk : input meshes" << endl;
  cout << "  -o result.vtk              : output mesh (template with initial momentum)" << endl;
  cout << "  -s sigma                   : kernel standard deviation" << endl;
  cout << "  -l lambda                  : weight of landmark distance term" << endl;
  cout << "Additional Options" << endl;
  cout << "  -d dim                     : problem dimension (3)" << endl;
  cout << "  -n N                       : number of time steps (100)" << endl;
  cout << "  -a <L|C|V>                 : data attachment term, L for landmark euclidean distance (default), " << endl;
  cout << "                               C for current metric, V for varifold metric." << endl;
  cout << "  -S sigma                   : kernel standard deviation for current/varifold metric" << endl;
  cout << "  -c mesh.vtk                : optional control point mesh (if different from template.vtk)" << endl;
  cout << "  -p array_name              : read initial momentum from named array in control/template mesh" << endl;
  cout << "  -i iter_grad iter_newt     : max iterations for optimization for gradient descent and newton's" << endl;
  cout << "  -O filepattern             : pattern for saving traced landmark paths (e.g., path%04d.vtk)" << endl;
  cout << "  -f                         : use single-precision float (off by deflt)" << endl;
  cout << "  -C mu0 mu_mult             : test constrained optimization (not for general use)" << endl;
  cout << "  -t n_threads               : limit number of concurrent threads to n_threads" << endl;
  cout << "  -D n                       : perform derivative check (for first n momenta)" << endl;
  cout << "  -L array_name              : use label-restricted data attachment, with label posteriors in given array" << endl;
  return -1;
}

void check(bool condition, const char *format,...)
{
  if(!condition)
    {
    char buffer[256];
    va_list args;
    va_start (args, format);
    vsprintf (buffer,format, args);
    va_end (args);

    cerr << buffer << endl;
    exit(-1);
    }
}

struct ShootingParameters
{
  enum Algorithm { Allassonniere, GradDescent, QuasiAllassonniere };
  enum DataAttachment { Euclidean, Current, Varifold };
  string fnTemplate, fnTarget, fnControlMesh;
  string fnOutput;
  string fnOutputPaths;
  string arrInitialMomentum;
  string arrAttachmentLabelPosteriors;
  double sigma = 0.0;
  double currents_sigma = 0.0;
  double lambda = 0.0;
  unsigned int dim = 3;
  unsigned int N = 100;
  unsigned int iter_grad = 20, iter_newton = 20;
  Algorithm alg = GradDescent;
  DataAttachment attach = Euclidean;
  bool use_float = false;
  unsigned int n_threads = 0;
  unsigned int n_deriv_check = 0;
  bool test_currents_attachment = false;

  // For constrained optimization - just exprimental
  double constrained_mu_init = 0.0, constrained_mu_mult = 0.0;
};


template <class TFloat, unsigned int VDim>
class TriangleCentersAndNormals
{
public:
  typedef PointSetHamiltonianSystem<TFloat, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;
  typedef vnl_matrix<int> Triangulation;

  TriangleCentersAndNormals(const Triangulation &tri, bool normalize);
  void Forward(const Matrix &q);
  void Backward(const Matrix &dE_dC, const Matrix &dE_dN, const Vector &dE_dW_norm, Matrix &dE_dq);
};

template <class TFloat>
class TriangleCentersAndNormals<TFloat, 3>
{
public:
  typedef PointSetHamiltonianSystem<TFloat, 3> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;
  typedef vnl_matrix<int> Triangulation;

  TriangleCentersAndNormals( const Triangulation &tri, bool normalize)
  {
    this->normalize = normalize;
    this->tri = tri;
    this->C.set_size(tri.rows(), 3);
    this->N.set_size(tri.rows(), 3);
    this->U.set_size(tri.rows(), 3);
    this->V.set_size(tri.rows(), 3);
    this->W.set_size(tri.rows(), 3);
    this->W_norm.set_size(tri.rows());
  }

  void Forward(const Matrix &q)
  {
    for(unsigned int i = 0; i < tri.rows(); i++)
      {
      // Get pointer access to the outputs
      TFloat *Ui = U.data_array()[i];
      TFloat *Vi = V.data_array()[i];
      TFloat *Wi = W.data_array()[i];
      TFloat *Ci = C.data_array()[i];
      TFloat *Ni = N.data_array()[i];

      int v0 = tri(i, 0), v1 = tri(i, 1), v2 = tri(i, 2);
      for(unsigned int a = 0; a < 3; a++)
        {
        Ci[a] = (q(v0, a) + q(v1, a) + q(v2, a)) / 3.0;
        Ui[a] = q(v1, a) - q(v0, a);
        Vi[a] = q(v2, a) - q(v0, a);
        }

      if(normalize)
        {
        Wi[0] = 0.5 * (Ui[1] * Vi[2] - Ui[2] * Vi[1]);
        Wi[1] = 0.5 * (Ui[2] * Vi[0] - Ui[0] * Vi[2]);
        Wi[2] = 0.5 * (Ui[0] * Vi[1] - Ui[1] * Vi[0]);

        // Compute the norm of the cross-product
        W_norm[i] = sqrt(Wi[0] * Wi[0] + Wi[1] * Wi[1] + Wi[2] * Wi[2]);
        if(W_norm[i] > 0.0)
          {
          Ni[0] = Wi[0] / W_norm[i];
          Ni[1] = Wi[1] / W_norm[i];
          Ni[2] = Wi[2] / W_norm[i];
          }
        else
          {
          Ni[0] = 0.0;
          Ni[1] = 0.0;
          Ni[2] = 0.0;
          }
        }
      else
        {
        // Compute the cross-product and store in the normal (this is what currents use)
        Ni[0] = 0.5 * (Ui[1] * Vi[2] - Ui[2] * Vi[1]);
        Ni[1] = 0.5 * (Ui[2] * Vi[0] - Ui[0] * Vi[2]);
        Ni[2] = 0.5 * (Ui[0] * Vi[1] - Ui[1] * Vi[0]);
        }
      }
  }

  void Backward(const Matrix &dE_dC, const Matrix &dE_dN, const Vector &dE_dW_norm, Matrix &dE_dq)
  {
    dE_dq.fill(0.0);
    TFloat dU[3], dV[3], dW[3];

    for(unsigned int i = 0; i < tri.rows(); i++)
      {
      // Get pointer access to the outputs
      TFloat *Ui = U.data_array()[i];
      TFloat *Vi = V.data_array()[i];
      TFloat *Wi = W.data_array()[i];
      TFloat *Ni = N.data_array()[i];
      const TFloat *dCi = dE_dC.data_array()[i];
      const TFloat *dNi = dE_dN.data_array()[i];
      TFloat dW_norm = dE_dW_norm[i];

      // Get the vertex indices and the corresponding gradients
      int v0 = tri(i, 0), v1 = tri(i, 1), v2 = tri(i, 2);
      TFloat *dq0 = dE_dq.data_array()[v0];
      TFloat *dq1 = dE_dq.data_array()[v1];
      TFloat *dq2 = dE_dq.data_array()[v2];

      if(normalize)
        {
        // Partial of the norm of W
        if(W_norm[i] > 0.0)
          {
          dW[0] = ((1 - Ni[0]*Ni[0]) * dNi[0] - Ni[0] * Ni[1] * dNi[1] - Ni[0] * Ni[2] * dNi[2] + Wi[0] * dW_norm) / W_norm[i];
          dW[1] = ((1 - Ni[1]*Ni[1]) * dNi[1] - Ni[1] * Ni[0] * dNi[0] - Ni[1] * Ni[2] * dNi[2] + Wi[1] * dW_norm) / W_norm[i];
          dW[2] = ((1 - Ni[2]*Ni[2]) * dNi[2] - Ni[2] * Ni[0] * dNi[0] - Ni[2] * Ni[1] * dNi[1] + Wi[2] * dW_norm) / W_norm[i];
          }
        else
          {
          dW[0] = dNi[0];
          dW[1] = dNi[1];
          dW[2] = dNi[2];
          }

        // Backprop the cross-product
        dU[0] = 0.5 * (Vi[1] * dW[2] - Vi[2] * dW[1]);
        dU[1] = 0.5 * (Vi[2] * dW[0] - Vi[0] * dW[2]);
        dU[2] = 0.5 * (Vi[0] * dW[1] - Vi[1] * dW[0]);

        dV[0] = 0.5 * (Ui[2] * dW[1] - Ui[1] * dW[2]);
        dV[1] = 0.5 * (Ui[0] * dW[2] - Ui[2] * dW[0]);
        dV[2] = 0.5 * (Ui[1] * dW[0] - Ui[0] * dW[1]);
        }
      else
        {
        // Backprop the cross-product
        dU[0] = 0.5 * (Vi[1] * dNi[2] - Vi[2] * dNi[1]);
        dU[1] = 0.5 * (Vi[2] * dNi[0] - Vi[0] * dNi[2]);
        dU[2] = 0.5 * (Vi[0] * dNi[1] - Vi[1] * dNi[0]);

        dV[0] = 0.5 * (Ui[2] * dNi[1] - Ui[1] * dNi[2]);
        dV[1] = 0.5 * (Ui[0] * dNi[2] - Ui[2] * dNi[0]);
        dV[2] = 0.5 * (Ui[1] * dNi[0] - Ui[0] * dNi[1]);
        }

      // Backprop the Ui and Vi
      for(unsigned int a = 0; a < 3; a++)
        {
        // The center contribution
        TFloat dCi_a  = dCi[a] / 3.0;
        dq0[a] += dCi_a - dU[a] - dV[a];
        dq1[a] += dCi_a + dU[a];
        dq2[a] += dCi_a + dV[a];
        }
      }
  }

  // Whether or not we normalize the values, normalization
  // is needed for the Varifold metric but not for Currents
  bool normalize;

  // Intermediate values
  Triangulation tri;
  Matrix U, V, W;
  Vector W_norm;

  // Results
  Matrix C, N;
};


template <class TFloat>
class TriangleCentersAndNormals<TFloat, 2>
{
public:
  typedef PointSetHamiltonianSystem<TFloat, 2> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;
  typedef vnl_matrix<int> Triangulation;

  TriangleCentersAndNormals( const Triangulation &tri, bool normalize)
  {
    this->normalize = normalize;
    this->tri = tri;
    this->C.set_size(tri.rows(), 2);
    this->N.set_size(tri.rows(), 2);
    this->U.set_size(tri.rows(), 2);
    this->W.set_size(tri.rows(), 2);
    this->W_norm.set_size(tri.rows());
  }

  void Forward(const Matrix &q)
  {
    for(unsigned int i = 0; i < tri.rows(); i++)
      {
      // Get pointer access to the outputs
      TFloat *Ui = U.data_array()[i];
      TFloat *Wi = W.data_array()[i];
      TFloat *Ci = C.data_array()[i];
      TFloat *Ni = N.data_array()[i];

      int v0 = tri(i, 0), v1 = tri(i, 1);
      for(unsigned int a = 0; a < 2; a++)
        {
        Ci[a] = (q(v0, a) + q(v1, a)) / 2.0;
        Ui[a] = q(v1, a) - q(v0, a);
        }

      if(normalize)
        {
        // Compute the cross-product
        Wi[0] = Ui[1];
        Wi[1] = -Ui[0];

        // Compute the norm of the cross-product
        W_norm[i] = sqrt(Wi[0] * Wi[0] + Wi[1] * Wi[1]);
        if(W_norm[i] > 0.0)
          {
          Ni[0] = Wi[0] / W_norm[i];
          Ni[1] = Wi[1] / W_norm[i];
          }
        else
          {
          Ni[0] = 0.0;
          Ni[1] = 0.0;
          }
        }
      else
        {
        // Compute the cross-product
        Ni[0] = Ui[1];
        Ni[1] = -Ui[0];
        }
      }
  }

  void Backward(const Matrix &dE_dC, const Matrix &dE_dN, const Vector &dE_dW_norm, Matrix &dE_dq)
  {
    dE_dq.fill(0.0);
    TFloat dU[2], dW[2];

    for(unsigned int i = 0; i < tri.rows(); i++)
      {
      // Get pointer access to the outputs
      TFloat *Ni = N.data_array()[i];
      const TFloat *dCi = dE_dC.data_array()[i];
      const TFloat *dNi = dE_dN.data_array()[i];
      TFloat *Wi = W.data_array()[i];
      TFloat dW_norm = dE_dW_norm[i];

      // Get the vertex indices and the corresponding gradients
      int v0 = tri(i, 0), v1 = tri(i, 1);
      TFloat *dq0 = dE_dq.data_array()[v0];
      TFloat *dq1 = dE_dq.data_array()[v1];

      if(normalize)
      {
        // Partial of the norm of W
        if(W_norm[i] > 0.0)
          {
          dW[0] = ((1 - Ni[0]*Ni[0]) * dNi[0] - Ni[0] * Ni[1] * dNi[1] + Wi[0] * dW_norm) / W_norm[i];
          dW[1] = ((1 - Ni[1]*Ni[1]) * dNi[1] - Ni[1] * Ni[0] * dNi[0] + Wi[0] * dW_norm) / W_norm[i];
          }
        else
          {
          dW[0] = dNi[0];
          dW[1] = dNi[1];
          }

        // Backprop the cross-product
        dU[0] = -dW[1];
        dU[1] = dW[0];
        }
      else
        {
        // Backprop the cross-product
        dU[0] = -dNi[1];
        dU[1] = dNi[0];
        }

      // Backprop the Ui and Vi
      for(unsigned int a = 0; a < 2; a++)
        {
        // The center contribution
        TFloat dCi_a  = dCi[a] / 2.0;
        dq0[a] += dCi_a - dU[a];
        dq1[a] += dCi_a + dU[a];
        }
      }
  }

  // Intermediate values
  bool normalize;
  Triangulation tri;
  Matrix U, W;
  Vector W_norm;

  // Results
  Matrix C, N;
};


template <class TFloat, unsigned int VDim>
class CurrentsAttachmentTerm
{
public:
  typedef PointSetHamiltonianSystem<TFloat, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;

  enum Mode { CURRENTS, VARIFOLD };

  /** Triangulation data type */
  typedef vnl_matrix<int> Triangulation;

  /* Thread-specific data for the computation of the currents norm */
  struct ThreadData {
    Matrix dE_dC, dE_dN;
    Vector dE_dW;
    std::vector<unsigned int> rows;
  };

  /**
   * Intermediate and output data for calls to Currents norm and scalar product
   * computations
   */
  struct CurrentScalarProductData
  {
    // Partials with respect to centers and normals
    Matrix dE_dC, dE_dN;

    // Partials with respect to the weights/areas (varifold only)
    Vector dE_dW;

    // Per-vertex energy component
    Vector z;

    // Per-thread data
    std::vector<ThreadData> td;
  };


  /**
   * Constructor
   *   mode         : Mode (currents or varifolds)
   *   m            : Number of total landmarks (control and template)
   *   qT           : Target vertices
   *   tri_template : Triangles (3D) or edges (2D) of the template
   *   tri_target   : Triangles (3D) or edges (2D) of the target
   *   lab_template : Template label posterior array (per triangle)
   *   lab_target   : Target label posterior array (per triangle)
   */
  CurrentsAttachmentTerm(Mode mode, unsigned int m, const Matrix &qT,
                         const Triangulation &tri_template,
                         const Triangulation &tri_target,
                         const Matrix &lab_template, const Matrix &lab_target,
                         TFloat sigma, unsigned int n_threads)
    : tcan_template(tri_template, mode==VARIFOLD), tcan_target(tri_target, mode==VARIFOLD)
  {
    this->mode = mode;
    this->m = m;
    this->tri_target = tri_target;
    this->tri_template = tri_template;
    this->sigma = sigma;
    this->n_threads = n_threads;
    this->lab_template = lab_template;
    this->lab_target = lab_target;

    // Compute the triangle centers once and for all
    tcan_target.Forward(qT);

    // Setup the data structures for current gradient storage
    SetupCurrentScalarProductData(&this->tri_target, &this->cspd_target);
    SetupCurrentScalarProductData(&this->tri_template, &this->cspd_template);

    // Compute the current norm for the target (fixed quantity)
    cspd_target.z.fill(0.0);
    ComputeCurrentHalfNormSquared(&tcan_target, &cspd_target, lab_target, false);

    // Allocate this norm among template trianges equally (so that the z array
    // after the complete current computation makes sense on a per-triangle
    // basis and can be exported on a mesh surface)
    z0_template.set_size(tri_template.rows());
    z0_template.fill(cspd_target.z.sum() / tri_template.rows());
  }

  /**
   * Compute the 1/2 currents scalar product of a mesh with itself with optional
   * partial derivatives with respect to simplex centers and normals
   */
  double ComputeCurrentHalfNormSquared(
      TriangleCentersAndNormals<TFloat, VDim> *tcan,
      CurrentScalarProductData *cspd,
      const Matrix &label_matrix,
      bool grad)
  {
    // Perform the pairwise computation multi-threaded
    std::vector<std::thread> workers;
    for (unsigned int thr = 0; thr < n_threads; thr++)
      {
      workers.push_back(std::thread([this,thr,tcan,cspd,label_matrix,grad]()
        {
        // Get our thread data
        ThreadData &my_td = cspd->td[thr];

        // Clear the gradients
        if(grad)
          {
          my_td.dE_dC.fill(0.0);
          my_td.dE_dN.fill(0.0);
          my_td.dE_dW.fill(0.0);
          }

        // Get data arrays for fast memory access
        auto c_da = tcan->C.data_array();
        auto n_da = tcan->N.data_array();
        auto w_da = tcan->W_norm.data_block();
        auto dc_da = my_td.dE_dC.data_array();
        auto dn_da = my_td.dE_dN.data_array();
        auto dw_da = my_td.dE_dW.data_block();
        auto z_da = cspd->z.data_block();
        auto l_da = label_matrix.data_array();

        TFloat f = -0.5 / (sigma * sigma);
        TFloat f_times_2 = 2.0 * f;

        int n_labels = label_matrix.columns();

        // Handle triangles for this thread
        for(unsigned int i : my_td.rows)
          {
          TFloat zi = 0.0;
          TFloat *ci = c_da[i], *ni = n_da[i], wi = w_da[i];
          TFloat *d_ci = dc_da[i], *d_ni = dn_da[i], &d_wi = dw_da[i];
          const TFloat *l_i = l_da[i];
          TFloat dq[VDim];

          TFloat label_weight = 0.0;
          for(int l = 0; l < n_labels; l++)
            label_weight += l_i[l]  * l_i[l];

          // Handle the diagonal term
          if(mode == CURRENTS)
            {
            for(unsigned int a = 0; a < VDim; a++)
              zi += 0.5 * label_weight * ni[a] * ni[a];

            // Handle the diagonal term in the gradient
            if(grad)
              for(unsigned int a = 0; a < VDim; a++)
                d_ni[a] += label_weight * ni[a];
            }
          else
            {
            // We know that the dot product of the normal with itself is just one
            // so we don't need to do anything with the normals. We just need to
            // compute the product of the weights
            zi += 0.5 * label_weight * wi * wi;
            if(grad)
              d_wi += label_weight * wi;
            }

          // Handle the off-diagonal terms
          for(unsigned int j = i+1; j < tcan->C.rows(); j++)
            {
            // Compute the kernel
            TFloat *cj = c_da[j], *nj = n_da[j], wj = w_da[j];
            TFloat *d_cj = dc_da[j], *d_nj = dn_da[j], &d_wj = dw_da[j];
            const TFloat *l_j = l_da[j];
            TFloat dist_sq = 0.0;
            TFloat dot_ni_nj = 0.0;
            for(unsigned int a = 0; a < VDim; a++)
              {
              dq[a] = ci[a] - cj[a];
              dist_sq += dq[a] * dq[a];
              dot_ni_nj += ni[a] * nj[a];
              }

            TFloat label_weight = 0.0;
            for(int l = 0; l < n_labels; l++)
              label_weight += l_i[l] * l_j[l];

            TFloat K = label_weight * exp(f * dist_sq);

            if(mode == CURRENTS)
              {
              TFloat zij = K * dot_ni_nj;
              zi += zij;

              if(grad)
                {
                TFloat w = f_times_2 * zij;
                for(unsigned int a = 0; a < VDim; a++)
                  {
                  d_ci[a] += w * dq[a];
                  d_cj[a] -= w * dq[a];
                  d_ni[a] += K * nj[a];
                  d_nj[a] += K * ni[a];
                  }
                }
              }
            else
              {
              TFloat K_wi_wj = K * wi * wj;
              TFloat dot_ni_nj_sq = dot_ni_nj * dot_ni_nj;
              TFloat zij = K_wi_wj * dot_ni_nj_sq;
              zi += zij;

              if(grad)
                {
                TFloat w = f_times_2 * zij;
                for(unsigned int a = 0; a < VDim; a++)
                  {
                  d_ci[a] += w * dq[a];
                  d_cj[a] -= w * dq[a];
                  d_ni[a] += 2 * dot_ni_nj * K_wi_wj * nj[a];
                  d_nj[a] += 2 * dot_ni_nj * K_wi_wj * ni[a];
                  }
                d_wi += K * wj * dot_ni_nj_sq;
                d_wj += K * wi * dot_ni_nj_sq;
                }
              }
            }

          // Store the cost for this template vertex
          z_da[i] += zi;
          }
        }));
      }

    // Run threads and halt until completed
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });

    // Accumulate the thread data
    for(unsigned int i = 0; i < n_threads; i++)
      {
      cspd->dE_dC += cspd->td[i].dE_dC;
      cspd->dE_dN += cspd->td[i].dE_dN;
      cspd->dE_dW += cspd->td[i].dE_dW;
      }
  }

  /**
   * Compute the currents scalar product of two meshes with with optional
   * partial derivatives with respect to simplex centers and normals of the
   * first input
   */
  double ComputeCurrentScalarProduct(
      TriangleCentersAndNormals<TFloat, VDim> *tcan_1,
      TriangleCentersAndNormals<TFloat, VDim> *tcan_2,
      CurrentScalarProductData *cspd_1,
      const Matrix &label_matrix1, const Matrix &label_matrix2,
      bool grad)
  {
    // Perform the pairwise computation multi-threaded
    std::vector<std::thread> workers;
    for (unsigned int thr = 0; thr < n_threads; thr++)
      {
      workers.push_back(std::thread([this,thr,tcan_1,tcan_2,cspd_1,label_matrix1,label_matrix2,grad]()
        {
        // Get data arrays for fast memory access
        auto c1_da = tcan_1->C.data_array(), c2_da = tcan_2->C.data_array();
        auto n1_da = tcan_1->N.data_array(), n2_da = tcan_2->N.data_array();
        auto w1_da = tcan_1->W_norm.data_block(), w2_da = tcan_2->W_norm.data_block();
        auto dc_da = cspd_1->dE_dC.data_array(), dn_da = cspd_1->dE_dN.data_array();
        auto dw_da = cspd_1->dE_dW.data_block();
        auto z1_da = cspd_1->z.data_block();
        auto l1_da = label_matrix1.data_array(), l2_da = label_matrix2.data_array();

        TFloat f = -0.5 / (sigma * sigma);
        TFloat f_times_2 = 2.0 * f;

        int n_labels = label_matrix1.columns();

        // Handle triangles for this thread
        for(unsigned int i = thr; i < tcan_1->C.rows(); i+=n_threads)
          {
          TFloat zi = 0.0;
          TFloat *ci = c1_da[i], *ni = n1_da[i], wi = w1_da[i];
          TFloat *d_ci = dc_da[i], *d_ni = dn_da[i], &d_wi = dw_da[i];
          TFloat dq[VDim];
          const TFloat *l_i = l1_da[i];

          for(unsigned int j = 0; j < tcan_2->C.rows(); j++)
            {
            // Compute the kernel
            TFloat *cj = c2_da[j], *nj = n2_da[j], wj = w2_da[j];
            const TFloat *l_j = l2_da[j];
            TFloat dist_sq = 0.0;
            TFloat dot_ni_nj = 0.0;
            for(unsigned int a = 0; a < VDim; a++)
              {
              dq[a] = ci[a] - cj[a];
              dist_sq += dq[a] * dq[a];
              dot_ni_nj += ni[a] * nj[a];
              }

            TFloat label_weight = 0.0;
            for(int l = 0; l < n_labels; l++)
              label_weight += l_i[l] * l_j[l];

            TFloat K = -label_weight * exp(f * dist_sq);

            if(mode == CURRENTS)
              {
              TFloat zij = K * dot_ni_nj;
              zi += zij;

              if(grad)
                {
                TFloat w = f_times_2 * zij;
                for(unsigned int a = 0; a < VDim; a++)
                  {
                  d_ci[a] += w * dq[a];
                  d_ni[a] += K * nj[a];
                  }
                }
              }
            else
              {
              TFloat K_wi_wj = K * wi * wj;
              TFloat dot_ni_nj_sq = dot_ni_nj * dot_ni_nj;
              TFloat zij = K_wi_wj * dot_ni_nj_sq;
              zi += zij;

              if(grad)
                {
                TFloat w = f_times_2 * zij;
                for(unsigned int a = 0; a < VDim; a++)
                  {
                  d_ci[a] += w * dq[a];
                  d_ni[a] += 2 * dot_ni_nj * K_wi_wj * nj[a];
                  }
                d_wi += K * wj * dot_ni_nj_sq;
                }
              }
            }

          // Store the cost for this template vertex
          z1_da[i] += zi;
          }
        }));
      }

    // Run threads and halt until completed
    std::for_each(workers.begin(), workers.end(), [](std::thread &t) { t.join(); });
  }

  /** Compute the energy and optional gradient of the energy */
  double Compute(const Matrix &q1, Matrix *grad = nullptr)
  {
    // Compute the triangle centers and normals
    tcan_template.Forward(q1);

    // Clear gradient terms if needed
    if(grad)
      {
      cspd_template.dE_dC.fill(0.0);
      cspd_template.dE_dN.fill(0.0);
      cspd_template.dE_dW.fill(0.0);
      }

    // Initialize the z-term
    cspd_template.z = z0_template;
    // double v0 = cspd_template.z.sum();

    // Add the squared norm term
    this->ComputeCurrentHalfNormSquared(&tcan_template, &cspd_template, lab_template, grad != nullptr);
    // double v1 = cspd_template.z.sum();

    // Subtract twice the scalar product term
    this->ComputeCurrentScalarProduct(&tcan_template, &tcan_target, &cspd_template, lab_template, lab_target, grad != nullptr);
    // double v2 = cspd_template.z.sum();

    // printf("0.5*S*S=%f, 0.5*T*T=%f, S*T=%f\n", v0, v1-v0, v2-v1);

    // Backpropagate the gradient to get gradient with respect to q1
    if(grad)
      tcan_template.Backward(cspd_template.dE_dC, cspd_template.dE_dN, cspd_template.dE_dW, *grad);

    // Return the total energy
    return cspd_template.z.sum();
  }

  /** Save the mesh representing current energy term */
  void SaveMesh(const Matrix &q1, const char *fname)
  {
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();

    this->Compute(q1);

    // Assign points
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    for(unsigned int i = 0; i < m; i++)
      pts->InsertNextPoint(q1.data_array()[i]);
    pd->SetPoints(pts);

    // Assign polys
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
    for(unsigned int j = 0; j < tcan_template.tri.rows(); j++)
      {
      vtkIdType ids[VDim];
      for(unsigned int a = 0; a < VDim; a++)
        ids[a] = tcan_template.tri(j, a);
      polys->InsertNextCell(VDim, ids);
      }
    if(VDim == 3)
      pd->SetPolys(polys);
    else
      pd->SetLines(polys);

    // Create the energy array
    vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
    arr->SetName("CurrentEnergy");
    arr->SetNumberOfComponents(1);
    arr->SetNumberOfTuples(tcan_template.tri.rows());
    for(unsigned int j = 0; j < tcan_template.tri.rows(); j++)
      arr->SetComponent(j, 0, cspd_template.z[j]);
    pd->GetCellData()->AddArray(arr);

    WriteVTKData(pd, fname);
  }

  /**
   * Initialize the data for scalar product computation (target or template)
   */
  void SetupCurrentScalarProductData(const Triangulation *tri,
                                     CurrentScalarProductData *cspd)
  {
    // Global objects
    cspd->dE_dC.set_size(tri->rows(), VDim);
    cspd->dE_dN.set_size(tri->rows(), VDim);
    cspd->dE_dW.set_size(tri->rows());
    cspd->z.set_size(tri->rows());
    cspd->td.resize(n_threads);

    // Per-thread objects
    for(unsigned int t = 0; t < n_threads; t++)
      {
      cspd->td[t].dE_dC.set_size(tri->rows(), VDim);
      cspd->td[t].dE_dN.set_size(tri->rows(), VDim);
      cspd->td[t].dE_dW.set_size(tri->rows());
      }

    // Assign lines in pairs, one at the top of the symmetric matrix K and
    // one at the bottom of K. The loop below will not assign the middle
    // line when there is an odd number of points (e.g., line 7 when there are 15)
    for(int i = 0; i < (int) tri->rows()/2; i++)
      {
      int i_thread = i % n_threads;
      cspd->td[i_thread].rows.push_back(i);
      cspd->td[i_thread].rows.push_back((tri->rows()-1) - i);
      }

    // Handle the middle line for odd number of vertices
    if(tri->rows() % 2 == 1)
      cspd->td[(tri->rows() / 2) % n_threads].rows.push_back(tri->rows()/2);
  }

protected:
  unsigned int m;
  Triangulation tri_template, tri_target;

  // Triangle quantity computer
  TriangleCentersAndNormals<TFloat, VDim> tcan_template, tcan_target;

  // Template and target current norm/scalar product data
  CurrentScalarProductData cspd_template, cspd_target;

  // Label posteriors
  Matrix lab_template, lab_target;

  // Current squared norm of the target allocated equally between all
  // the trianges in the template
  Vector z0_template;

  // Kernel sigma
  TFloat sigma;

  // Threading
  unsigned int n_threads;

  // Mode
  Mode mode;
};



template <class TFloat, unsigned int VDim>
class PointSetShootingCostFunction : public vnl_cost_function
{
public:
  typedef PointSetHamiltonianSystem<TFloat, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;
  typedef vnl_matrix<int> Triangulation;

  // Separate type because vnl optimizer is double-only
  typedef vnl_vector<double> DVector;

  PointSetShootingCostFunction(
    const ShootingParameters &param,
      const Matrix &q0, const Matrix &p0, const Matrix &qT,
      Triangulation tri_template,
      Triangulation tri_target,
      const Matrix &lab_template, const Matrix &lab_target)
    : vnl_cost_function(p0.rows() * VDim),
      hsys(q0, param.sigma, param.N, q0.rows() - p0.rows())
    {
    this->p0 = p0;
    this->q0 = q0;
    this->qT = qT;
    this->param = param;
    this->k = p0.rows();
    this->m = q0.rows();
    this->p1.set_size(k,VDim);
    this->q1.set_size(m,VDim);

    for(unsigned int a = 0; a < VDim; a++)
      {
      alpha[a].set_size(m);
      beta[a].set_size(k); beta[a].fill(0.0);
      grad_f[a].set_size(m);
      }

    // Set up the currents attachment
    currents_attachment = nullptr;
    if(param.attach == ShootingParameters::Current || param.attach == ShootingParameters::Varifold)
      {
      // Create the appropriate attachment term (currents or varifolds)
      currents_attachment = new CATerm(
            param.attach == ShootingParameters::Current ? CATerm::CURRENTS : CATerm::VARIFOLD,
            m, qT, tri_template, tri_target, lab_template, lab_target,
            param.currents_sigma, param.n_threads);

      // Allocate the gradient storage
      grad_currents.set_size(m, VDim);
      }
    }

  ~PointSetShootingCostFunction()
  {
    if(currents_attachment)
      delete currents_attachment;
  }

  DVector wide_to_tall(const Vector p[VDim])
    {
    DVector v(p[0].size() * VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        v[pos++] = p[a](i);
    return v;
    }

  DVector wide_to_tall(const Matrix &p)
    {
    DVector v(p.rows() * VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        v[pos++] = p(i,a);
    return v;
    }

  Matrix tall_to_wide(const DVector &v)
    {
    Matrix p(v.size() / VDim, VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        p(i,a) = (TFloat) v[pos++];
    return p;
    }

  virtual double ComputeEuclideanAttachment()
  {
    // Compute the landmark errors
    double E_data = 0.0;
    unsigned int i_temp = (k == m) ? 0 : k;
    for(unsigned int a = 0; a < VDim; a++)
      {
      alpha[a].fill(0.0);
      for(unsigned int i = i_temp; i < m; i++)
        {
        alpha[a](i) = q1(i,a) - qT(i - i_temp, a);
        E_data += 0.5 * alpha[a](i) * alpha[a](i);
        }
      }

    return E_data;
  }

 virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g)
    {
    // Initialize the p0-vector
    p0 = tall_to_wide(x);

    // Perform flow
    double H = hsys.FlowHamiltonian(p0, q1, p1);

    // Compute the data attachment term
    double E_data = 0.0;
    if(param.attach == ShootingParameters::Euclidean)
      E_data = ComputeEuclideanAttachment();
    else if(param.attach == ShootingParameters::Current || param.attach == ShootingParameters::Varifold)
      {
      static int my_iter = 0;
      if(g)
        {
        E_data = currents_attachment->Compute(q1, &grad_currents);
        for(unsigned int i = 0; i < m; i++)
          for(unsigned int a = 0; a < VDim; a++)
            alpha[a][i] = grad_currents(i,a);
        }
      else
        E_data = currents_attachment->Compute(q1);
      }

    if(f)
      *f = H + param.lambda * E_data;

    if(g)
      {
      // Multiply gradient of f. wrt q1 (alpha) by Jacobian of q1 wrt p0
      hsys.FlowGradientBackward(alpha, beta, grad_f);

      // Recompute Hq/Hp at initial timepoint (TODO: why are we doing this?)
      hsys.ComputeHamiltonianJet(q0, p0, false);

      // Complete gradient computation
      for(unsigned int a = 0; a < VDim; a++)
        {
        // Combine the gradient terms
        grad_f[a] = grad_f[a] * param.lambda + hsys.GetHp(a).extract(k);
        }

      // Pack the gradient into the output vector
      *g = wide_to_tall(grad_f);

      // Count this as an iteration
      ++iter;
      }

    // Print current results
    if(verbose && g && f) 
      {
      printf("It = %04d  H = %8.2f  DA = %8.2f  f = %8.2f\n", iter, H, E_data * param.lambda, *f);
      }
    }

 void SetVerbose(bool value) 
   { 
   this->verbose = value;
   }

protected:
  HSystem hsys;
  ShootingParameters param;
  Matrix qT, p0, q0, p1, q1;
  Vector alpha[VDim], beta[VDim], grad_f[VDim];

  // Attachment terms
  typedef CurrentsAttachmentTerm<TFloat, VDim> CATerm;
  CATerm *currents_attachment;

  // For currents, the gradient is supplied as a matrix
  Matrix grad_currents;

  // Number of control points (k) and total points (m)
  unsigned int k, m;

  // Whether to print values at each iteration
  bool verbose = false;
  unsigned int iter = 0;
};



template <class TFloat, unsigned int VDim>
class PointSetShootingLineSearchCostFunction : public vnl_cost_function
{
public:
  typedef PointSetHamiltonianSystem<TFloat, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;

  // Separate type because vnl optimizer is double-only
  typedef vnl_vector<double> DVector;

  PointSetShootingLineSearchCostFunction(
    const ShootingParameters &param, const Matrix &q0, const Matrix &p0, const Matrix &qT, const Matrix &del_p0)
    : vnl_cost_function(q0.rows() * VDim), hsys(q0, param.sigma, param.N)
    {
    this->p0 = p0;
    this->del_p0 = del_p0;
    this->qT = qT;
    this->param = param;
    this->k = q0.rows();
    this->p1.set_size(k,VDim);
    this->q1.set_size(k,VDim);
    }

  virtual double f (vnl_vector<double> const& x)
    {
    TFloat alpha = (TFloat) x[0];

    // Perform flow
    double H = hsys.FlowHamiltonian(p0 + alpha * del_p0, q1, p1);

    // Compute the landmark errors
    double fnorm_sq = 0.0;
    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int i = 0; i < k; i++)
        {
        double d = q1(i,a) - qT(i,a);
        fnorm_sq += d * d;
        }
      }

    // Compute the landmark part of the objective
    double Edist = 0.5 * fnorm_sq;

    // cout << H + param.lambda * Edist << endl;
    return H + param.lambda * Edist;
    }



protected:
  HSystem hsys;
  ShootingParameters param;
  Matrix qT, p0, del_p0, q0, p1, q1;
  unsigned int k;
};



template <class TFloat, unsigned int VDim>
class PointSetShootingTransversalityCostFunction : public vnl_cost_function
{
public:
  typedef PointSetHamiltonianSystem<TFloat, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;

  // Separate type because vnl optimizer is double-only
  typedef vnl_vector<double> DVector;

  PointSetShootingTransversalityCostFunction(
    const ShootingParameters &param, const Matrix &q0, const Matrix &qT)
    : vnl_cost_function(q0.rows() * VDim), hsys(q0, param.sigma, param.N)
    {
    this->p0 = (qT - q0) / param.N;
    this->qT = qT;
    this->param = param;
    this->k = q0.rows();
    this->p1.set_size(k,VDim);
    this->q1.set_size(k,VDim);

    for(unsigned int a = 0; a < VDim; a++)
      {
      alpha[a].set_size(k);
      beta[a].set_size(k); beta[a].fill(0.0);
      G[a].set_size(k);
      grad_f[a].set_size(k);
      }
    }

  DVector wide_to_tall(const Vector p[VDim])
    {
    DVector v(k * VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        v[pos++] = p[a](i);
    return v;
    }

  DVector wide_to_tall(const Matrix &p)
    {
    DVector v(k * VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        v[pos++] = p(i,a);
    return v;
    }

  Matrix tall_to_wide(const DVector &v)
    {
    Matrix p(k,VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        p(i,a) = (TFloat) v[pos++];
    return p;
    }

  virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g)
    {
    // Initialize the p0-vector
    p0 = tall_to_wide(x);

    // Perform flow
    double H = hsys.FlowHamiltonian(p0, q1, p1);

    // Compute G and alpha/beta
    double Gnorm_sq = 0.0, dsq = 0.0;
    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int i = 0; i < k; i++)
        {
        G[a](i) = p1(i, a) + param.lambda * (q1(i, a) - qT(i, a));
        Gnorm_sq += G[a][i] * G[a][i];
        dsq += (q1(i, a) - qT(i, a)) * (q1(i, a) - qT(i, a));

        alpha[a](i) = param.lambda * G[a][i];
        beta[a](i) = G[a][i];
        }
      }

    if(f)
      *f = 0.5 * Gnorm_sq;

    if(g)
      {
      // Multiply gradient of f. wrt q1 (alpha) by Jacobian of q1 wrt p0
      hsys.FlowGradientBackward(alpha, beta, grad_f);

      // Pack the gradient into the output vector
      *g = wide_to_tall(grad_f);
      }

    // Print the current state
    printf("H=%8.6f   Edist=%8.6f   E=%8.6f   |G|=%8.6f\n",
      H, 0.5 * param.lambda * dsq, H + 0.5 * param.lambda * dsq, sqrt(Gnorm_sq));

    }



protected:
  HSystem hsys;
  ShootingParameters param;
  Matrix qT, p0, q0, p1, q1;
  Vector alpha[VDim], beta[VDim], G[VDim], grad_f[VDim];
  unsigned int k;

};


template <class TFloat, unsigned int VDim>
class PointSetShootingProblem
{
public:
  typedef PointSetHamiltonianSystem<TFloat, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;
  typedef vnl_matrix<int> Triangulation;

  // Minimize using the transversality principle
  static void minimize_Allassonniere(const ShootingParameters &param, 
    const Matrix &q0, const Matrix &qT, Matrix &p0);

  static void minimize_QuasiAllassonniere(const ShootingParameters &param,
    const Matrix &q0, const Matrix &qT, Matrix &p0);

  // Minimize using gradient descent
  static void minimize_gradient(
      const ShootingParameters &param,
      const Matrix &q0, const Matrix &qT, Matrix &p0,
      const Triangulation &tri_template, const Triangulation &tri_target,
      const Matrix &lab_template, const Matrix &lab_target);

  static int minimize(const ShootingParameters &param);

private:
  static int TestCurrentsAttachmentTerm(
      const ShootingParameters &param,
      Matrix &q0, Matrix &qT,
      vnl_matrix<int> &tri_template, vnl_matrix<int> &tri_target,
      const Matrix &lab_template, const Matrix &lab_target);
};

#include <vnl/algo/vnl_brent_minimizer.h>

template <class TFloat, unsigned int VDim>
void
PointSetShootingProblem<TFloat, VDim>
::minimize_Allassonniere(const ShootingParameters &param, 
  const Matrix &q0, const Matrix &qT, Matrix &p0)
{
  unsigned int k = q0.rows();

  // Create the hamiltonian system
  HSystem hsys(q0, param.sigma, param.N, 0, param.n_threads);

  // Where to store the results of the flow
  Matrix q1(k,VDim), p1(k,VDim), del_p0(k, VDim), grad_q[VDim][VDim], grad_p[VDim][VDim];
  for(unsigned int a = 0; a < VDim; a++)
    {
    for(unsigned int b = 0; b < VDim; b++)
      {
      grad_p[a][b].set_size(k,k);
      grad_q[a][b].set_size(k,k);
      }
    }

  // Transversality term and Hessian-like term
  Vector G(VDim * k);
  Matrix DG(VDim * k, VDim * k);

  // Perform optimization using the Allassonniere method
  for(unsigned int iter = 0; iter < param.iter_newton; iter++)
    {
    // Perform Hamiltonian flow
    double H = hsys.FlowHamiltonianWithGradient(p0, q1, p1, grad_q, grad_p);

    // Compute the landmark errors
    Matrix lmdiff = q1 - qT;

    // Compute the landmark part of the objective
    double fnorm = lmdiff.frobenius_norm();
    double dsq = fnorm * fnorm;

    // Fill G and Hessian
    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int i = 0; i < k; i++)
        {
        // Compute the transversality vector error vector G
        G(a * k + i) = p1(i,a) + 2 * param.lambda * lmdiff(i,a);

        // Fill the Hessian-like matrix
        for(unsigned int b = 0; b < VDim; b++)
          {
          for(unsigned int j = 0; j < k; j++)
            {
            DG(a * k + i, b * k + j) = grad_p[a][b](i,j) + 2 * param.lambda * grad_q[a][b](i,j);
            }
          }
        }
      }

    // Perform singular value decomposition on the Hessian matrix, zeroing
    // out the singular values below 1.0 (TODO: use a relative threshold?)
    vnl_svd<TFloat> svd(DG, -0.001);

    int nnz = 0;
    for(int i = 0; i < svd.W().rows(); i++)
      if(svd.W()(i,i) != 0.0)
        nnz++;

    printf("SVD min: %12.8f, max: %12.8f, nnz: %d, rank: %d\n", 
      svd.sigma_min(), svd.sigma_max(), nnz, svd.rank());

    // Compute inv(DG) * G
    Vector del_p0_vec = - svd.solve(G);
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        del_p0(i,a) = del_p0_vec(a * k + i);

    // Perform line search - turned out useless
    /*
    typedef PointSetShootingLineSearchCostFunction<TFloat, VDim> CostFn;
    CostFn cost_fn(param, q0, p0, qT, del_p0);

    vnl_brent_minimizer brent(cost_fn);
    brent.set_x_tolerance(0.02);
    // TFloat alpha = brent.minimize_given_bounds(0.0, 0.9, 2.0);
    */

    // Argh! What to do with alpha!
    TFloat alpha = 0.1; // 1.0 - pow( 0.9, iter + 1.0);

    // Print the current state
    printf("Iter %4d   H=%8.6f   Edist=%8.6f   E=%8.6f   |G|=%8.6f   alpha=%8.6f\n",
      iter, H, 0.5 * param.lambda * dsq, H + 0.5 * param.lambda * dsq, G.two_norm(), alpha);

    p0 += alpha * del_p0;
    }
}
/*
template <class TFloat, unsigned int VDim>
void
PointSetShootingProblem<TFloat, VDim>
::minimize_QuasiAllassonniere(const ShootingParameters &param,
  const Matrix &q0, const Matrix &qT, Matrix &p0)
{
  unsigned int k = q0.rows();

  // Create the hamiltonian system
  HSystem hsys(q0, param.sigma, param.N);

  // Where to store the results of the flow
  Matrix q1(k,VDim), p1(k,VDim), p0test(k,VDim);

  // Transversality term
  Vector G(VDim * k), Gnext(VDim * k);
  Vector dp0(VDim * k), yj(VDim * k);

  // Hessian inverse term
  Matrix Hj(VDim * k, VDim * k), Aj(VDim * k, VDim * k), Bj(VDim * k, VDim * k), Cj(VDim * k, VDim * k);
  Hj.set_identity();

  // Compute gradient G for j = 0
  double H = hsys.FlowHamiltonian(p0, q1, p1);
  for(unsigned int a = 0; a < VDim; a++)
    for(unsigned int i = 0; i < k; i++)
      G(a * k + i) = p1(i,a) + 2 * param.lambda * (q1(i,a) - qT(i,a));

  // Iterate
  TFloat alpha = 0.001;
  for(unsigned int iter = 0; iter < param.iter; iter++)
    {
    // Compute the displacement direction
    dp0 = - (Hj * G);

    // Find an alpha that satisfies Wolfe condition
    double targ = fabs(0.9 * dot_product(dp0, G));

    while(alpha > 1e-6)
      {
      // Update the current point
      for(unsigned int a = 0; a < VDim; a++)
        for(unsigned int i = 0; i < k; i++)
          p0test(i,a) = p0(i,a) + alpha * dp0(a * k + i);

      // Compute gradient at updated location
      H = hsys.FlowHamiltonian(p0test, q1, p1);
      for(unsigned int a = 0; a < VDim; a++)
        for(unsigned int i = 0; i < k; i++)
          Gnext(a * k + i) = p1(i,a) + 2 * param.lambda * (q1(i,a) - qT(i,a));

      double test = fabs(dot_product(dp0, Gnext));
      if(test < targ)
        break;
      else
        alpha = 0.5 * alpha;
      }
    if(alpha < 1e-6)
      {
      cerr << "Failed line search" << endl;
      break;
      }

    // Update p0
    p0 = p0test;

    // Compute yj - difference in gradient
    yj = Gnext - G;

    // Update the Hessian inverse matrix

    // Compute the update stuff
    Vector Z = alpha * dp0 - Hj * yj;
    double z = dot_product(Z, yj);

    for(unsigned int i = 0; i < k * VDim; i++)
      {
      for(unsigned int m = 0; m < k * VDim; m++)
        {
        // Hj(i,m) += Z(i) * Z(m) / z;
        }
      }

    // Store Gnext as G
    G = Gnext;

    // Compute the landmark part of the objective
    double fnorm = (q1 - qT).frobenius_norm();
    double dsq = fnorm * fnorm;

    // Print the current state
    printf("Iter %4d   H=%8.6f   l*Dsq=%8.6f   E=%8.6f   |G|=%8.6f\n",
      iter, H, param.lambda * dsq, H + param.lambda * dsq, G.two_norm());
    }
}
*/


#include "vnl/algo/vnl_lbfgsb.h"

template <class TFloat, unsigned int VDim>
void
PointSetShootingProblem<TFloat, VDim>
::minimize_QuasiAllassonniere(const ShootingParameters &param,
  const Matrix &q0, const Matrix &qT, Matrix &p0)
{
  // Create the minimization problem
  typedef PointSetShootingTransversalityCostFunction<TFloat, VDim> CostFn;
  CostFn cost_fn(param, q0, qT);

  // Create initial/final solution
  p0 = (qT - q0) / param.N;
  typename CostFn::DVector x = cost_fn.wide_to_tall(p0);


  // Solve the minimization problem

  vnl_lbfgsb optimizer(cost_fn);
  optimizer.set_f_tolerance(1e-9);
  optimizer.set_x_tolerance(1e-4);
  optimizer.set_g_tolerance(1e-6);
  optimizer.set_trace(true);
  optimizer.set_max_function_evals(param.iter_newton);
  optimizer.minimize(x);

  // Take the optimal solution
  p0 = cost_fn.tall_to_wide(x);
}



template <class TFloat, unsigned int VDim>
void
PointSetShootingProblem<TFloat, VDim>
::minimize_gradient(
    const ShootingParameters &param,
    const Matrix &q0, const Matrix &qT, Matrix &p0,
    const Triangulation &tri_template, const Triangulation &tri_target,
    const Matrix &lab_template, const Matrix &lab_target)
{
  // unsigned int k = q0.rows();

  // Create the minimization problem
  typedef PointSetShootingCostFunction<TFloat, VDim> CostFn;
  CostFn cost_fn(param, q0, p0, qT, tri_template, tri_target, lab_template, lab_target);

  // Create initial/final solution
  typename CostFn::DVector x = cost_fn.wide_to_tall(p0);

  // Uncomment this code to test derivative computation
  if(param.n_deriv_check > 0)
    {
    TFloat eps = 1e-6;
    typename CostFn::DVector test_grad(x.size());
    double f_test;
    cost_fn.compute(x, &f_test, &test_grad);
    for(unsigned int i = 0; i < std::min(p0.size(), param.n_deriv_check); i++)
      {
      typename CostFn::DVector xtest = x;
      double f1, f2;
      xtest[i] = x[i] - eps;
      cost_fn.compute(xtest, &f1, NULL);

      xtest[i] = x[i] + eps;
      cost_fn.compute(xtest, &f2, NULL);

      printf("i = %03d,  AG = %8.4f,  NG = %8.4f\n", i, test_grad[i], (f2 - f1) / (2 * eps));
      }
    }

  // Solve the minimization problem
  cost_fn.SetVerbose(true);

  vnl_lbfgsb optimizer(cost_fn);
  optimizer.set_f_tolerance(1e-9);
  optimizer.set_x_tolerance(1e-4);
  optimizer.set_g_tolerance(1e-6);
  optimizer.set_trace(false);
  optimizer.set_max_function_evals(param.iter_grad);

  // vnl_conjugate_gradient optimizer(cost_fn);
  // optimizer.set_trace(true);
  optimizer.minimize(x);

  // Take the optimal solution
  p0 = cost_fn.tall_to_wide(x);
}



template <class TFloat, unsigned int VDim>
int
PointSetShootingProblem<TFloat, VDim>
::TestCurrentsAttachmentTerm(const ShootingParameters &param,
                             Matrix &q0, Matrix &qT,
                             vnl_matrix<int> &tri_template, vnl_matrix<int> &tri_target,
                             const Matrix &lab_template, const Matrix &lab_target)
{
  int m = q0.rows();
  Matrix grad_currents(m, VDim);

  TriangleCentersAndNormals <TFloat, VDim> tcan(tri_template, true);
  tcan.Forward(q0);
  cout << "TCAN test" << endl;
  cout << tcan.C.get_row(333) << endl;
  cout << tcan.N.get_row(333) << endl;
  cout << tcan.W_norm(333) << endl;

  Matrix eC(tcan.C.rows(), VDim); eC.fill(1.0);
  Matrix eN(tcan.C.rows(), VDim); eN.fill(1.0);
  Vector ea(tcan.C.rows()); ea.fill(1.0);
  Matrix eQ(tcan.C.rows(), VDim); eC.fill(1.0);
  tcan.Backward(eC, eN, ea, eQ);
  cout << eQ.get_row(333) << endl;

  typedef CurrentsAttachmentTerm<TFloat, VDim> CATerm;
  CATerm currents_attachment(
        param.attach == ShootingParameters::Current ? CATerm::CURRENTS : CATerm::VARIFOLD,
        m, qT, tri_template, tri_target, lab_template, lab_target,
        param.currents_sigma, param.n_threads);

  double value = currents_attachment.Compute(q0, &grad_currents);
  printf("Currents Attachment Value: %f\n", value);
}


template <class TFloat, unsigned int VDim>
int
PointSetShootingProblem<TFloat, VDim>
::minimize(const ShootingParameters &param)
{
  // Read the datasets
  vtkPolyData *pTemplate = ReadVTKData(param.fnTemplate);
  vtkPolyData *pTarget = ReadVTKData(param.fnTarget);

  // Read the optional control point dataset
  vtkPolyData *pControl = nullptr;
  if(param.fnControlMesh.length())
    pControl = ReadVTKData(param.fnControlMesh);

  // Get the number of vertices and dimensionality
  if(param.attach == ShootingParameters::Euclidean)
    check(pTemplate->GetNumberOfPoints() == pTarget->GetNumberOfPoints(),
          "Template and target meshes must match for the Landmark attachment term");

  // Get the number of control points
  unsigned int k = pControl ? pControl->GetNumberOfPoints() : pTemplate->GetNumberOfPoints();

  // Get the number of non-control (rider) points
  unsigned int n_riders = pControl ? pTemplate->GetNumberOfPoints() : 0;

  // Total points
  unsigned int m = k + n_riders;

  // Report number of actual vertices being used
  cout << "Performing geodesic shooting with " << k << " control points and " << m << " total landmarks." << endl;

  // Landmarks and initial momentum
  Matrix q0(m,VDim), p0(k,VDim);

  // Initialize the q0 vertex array
  for(unsigned int a = 0; a < VDim; a++)
    {
    // Control point portion of q0
    for(unsigned int i = 0; i < k; i++)
      q0(i,a) = pControl ? pControl->GetPoint(i)[a] : pTemplate->GetPoint(i)[a];

    // Rider point portion of q0
    for(unsigned int i = k; i < m; i++)
      q0(i,a) = pTemplate->GetPoint(i-k)[a];
    }

  // Initialize the target array
  Matrix qT(pTarget->GetNumberOfPoints(), VDim);
  for(unsigned int a = 0; a < VDim; a++)
    {
    for(unsigned int i = 0; i < pTarget->GetNumberOfPoints(); i++)
      {
      qT(i,a) = pTarget->GetPoint(i)[a];

      // In the simple case of Euclidean matching with no riders, initialize the momentum based on the
      // point assignment
      if(m == k && param.attach == ShootingParameters::Euclidean && param.arrInitialMomentum.length() == 0)
        p0(i,a) = (qT(i,a) - q0(i,a)) / param.N;
      }
    }

  // Read the initial momentum array
  if(param.arrInitialMomentum.length())
    {
    vtkDataArray *da_p0 = nullptr;
    if(pControl)
      da_p0 = pControl->GetPointData()->GetArray(param.arrInitialMomentum.c_str());
    else
      da_p0 = pTemplate->GetPointData()->GetArray(param.arrInitialMomentum.c_str());

    check(da_p0 && da_p0->GetNumberOfTuples() == k && da_p0->GetNumberOfComponents() == VDim,
          "Initial momentum array missing or has wrong dimensions");

    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        p0(i,a) = da_p0->GetComponent(i, a);
    }

  // Compute the triangulation of the template for non-landmark metrics
  vnl_matrix<int> tri_template, tri_target;
  Matrix lab_template(pTemplate->GetNumberOfCells(), 1, 1.0);
  Matrix lab_target(pTarget->GetNumberOfCells(), 1, 1.0);

  if(param.attach == ShootingParameters::Current || param.attach == ShootingParameters::Varifold)
    {
    tri_template.set_size(pTemplate->GetNumberOfCells(), VDim);
    for(unsigned int i = 0; i < pTemplate->GetNumberOfCells(); i++)
      {
      if(pTemplate->GetCell(i)->GetNumberOfPoints() != VDim)
        {
        std::cerr << "Wrong number of points in template cell " << i << std::endl;
        return -1;
        }
      for(unsigned int a = 0; a < VDim; a++)
        {
        unsigned int j = pTemplate->GetCell(i)->GetPointId(a);
        tri_template(i,a) = pControl ? j + k : j;
        }
      }

    // Compute the triangulation of the target
    tri_target.set_size(pTarget->GetNumberOfCells(), VDim);
    for(unsigned int i = 0; i < pTarget->GetNumberOfCells(); i++)
      {
      if(pTarget->GetCell(i)->GetNumberOfPoints() != VDim)
        {
        std::cerr << "Wrong number of points in target cell " << i << std::endl;
        return -1;
        }
      for(unsigned int a = 0; a < VDim; a++)
        tri_target(i,a) = pTarget->GetCell(i)->GetPointId(a);
      }

    // Read the label arrays
    if(param.arrAttachmentLabelPosteriors.length())
      {
      vtkDataArray *da_template = pTemplate->GetCellData()->GetArray(param.arrAttachmentLabelPosteriors.c_str());
      vtkDataArray *da_target = pTarget->GetCellData()->GetArray(param.arrAttachmentLabelPosteriors.c_str());
      check(da_template && da_target && da_template->GetNumberOfComponents() == da_target->GetNumberOfComponents(),
            "Label posterior arrays in template and target missing or do not match");
      int n_labels = da_template->GetNumberOfComponents();
      lab_template.set_size(tri_template.rows(), n_labels);
      for(unsigned int i = 0; i < tri_template.rows(); i++)
        for(unsigned int l = 0; l < n_labels; l++)
          lab_template(i,l) = da_template->GetComponent(i,l);
      lab_target.set_size(tri_target.rows(), n_labels);
      for(unsigned int i = 0; i < tri_target.rows(); i++)
        for(unsigned int l = 0; l < n_labels; l++)
          lab_target(i,l) = da_target->GetComponent(i,l);
      }

    if(param.test_currents_attachment)
      {
      TestCurrentsAttachmentTerm(param, q0, qT, tri_template, tri_target, lab_template, lab_target);
      return 0;
      }
    }

  // Run some iterations of gradient descent
  if(param.iter_grad > 0)
    {
    minimize_gradient(param, q0, qT, p0, tri_template, tri_target, lab_template, lab_target);
    }

  if(param.iter_newton > 0)
    {
    minimize_Allassonniere(param, q0, qT, p0);
    }

  // Genererate the output momentum map
  vtkDoubleArray *arr_p = vtkDoubleArray::New();
  arr_p->SetNumberOfComponents(VDim);
  arr_p->SetNumberOfTuples(k);
  arr_p->SetName("InitialMomentum");

  for(unsigned int a = 0; a < VDim; a++)
    arr_p->FillComponent(a, 0);

  for(unsigned int a = 0; a < VDim; a++)
    for(unsigned int i = 0; i < k; i++)
      arr_p->SetComponent(i, a, p0(i,a));

  if(pControl)
    {
    pControl->GetPointData()->AddArray(arr_p);
    WriteVTKData(pControl, param.fnOutput);
    }
  else
    {
    pTemplate->GetPointData()->AddArray(arr_p);
    WriteVTKData(pTemplate, param.fnOutput);
    }

  // If saving paths requested
  if(param.fnOutputPaths.size())
    {
    // Create and flow a system
    HSystem hsys(q0, param.sigma, param.N, m - k, param.n_threads);
    Matrix q1, p1;
    hsys.FlowHamiltonian(p0, q1, p1);

    // Get the number of points in the output mesh
    unsigned int nv = pTemplate->GetNumberOfPoints();

    // Apply the flow to the points in the rest of the mesh
    vtkDoubleArray *arr_v = vtkDoubleArray::New();
    arr_v->SetNumberOfComponents(VDim);
    arr_v->SetNumberOfTuples(nv);
    arr_v->SetName("Velocity");
    pTemplate->GetPointData()->AddArray(arr_v);

    /*
    vtkDoubleArray *arr_p = vtkDoubleArray::New();
    arr_p->SetNumberOfComponents(VDim);
    arr_p->SetNumberOfTuples(np);
    arr_p->SetName("Momentum");
    pTemplate->GetPointData()->AddArray(arr_p);
    */

    // Apply Euler method to the mesh points
    double dt = hsys.GetDeltaT();
    for(unsigned int t = 1; t < param.N; t++)
      {
      for(unsigned int i = 0; i < nv; i++)
        {
        TFloat qi[VDim], vi[VDim];

        for(unsigned int a = 0; a < VDim; a++)
          qi[a] = pTemplate->GetPoint(i)[a];

        // Interpolate the velocity at each mesh point
        hsys.InterpolateVelocity(t-1, qi, vi);

        // Update the position using Euler's method
        for(unsigned int a = 0; a < VDim; a++)
          qi[a] += dt * vi[a];

        for(unsigned int a = 0; a < VDim; a++)
          pTemplate->GetPoints()->SetPoint(i, qi);

        for(unsigned int a = 0; a < VDim; a++)
          {
          arr_v->SetComponent(i, a, vi[a]);
          // arr_p->SetComponent(i, a, hsys.GetPt(t)(i,a));
          }
        }

      // Output the intermediate mesh
      char buffer[1024];
      sprintf(buffer, param.fnOutputPaths.c_str(), t);
      WriteVTKData(pTemplate, buffer);
      }
    }

  return 0;
}

int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();

  ShootingParameters param;

  // Process parameters
  CommandLineHelper cl(argc, argv);
  while(!cl.is_at_end())
    {
    // Read the next command
    std::string arg = cl.read_command();

    if(arg == "-m")
      {
      param.fnTemplate = cl.read_existing_filename();
      param.fnTarget = cl.read_existing_filename();
      }
    else if(arg == "-c")
      {
      param.fnControlMesh = cl.read_existing_filename();
      }
    else if(arg == "-o")
      {
      param.fnOutput = cl.read_output_filename();
      }
    else if(arg == "-O")
      {
      param.fnOutputPaths = cl.read_string();
      }
    else if(arg == "-s")
      {
      param.sigma = cl.read_double();
      }
    else if(arg == "-l")
      {
      param.lambda = cl.read_double();
      }
    else if(arg == "-n")
      {
      param.N = (unsigned int) cl.read_integer();
      }
    else if(arg == "-d")
      {
      param.dim = (unsigned int) cl.read_integer();
      }
    else if(arg == "-i")
      {
      param.iter_grad = (unsigned int) cl.read_integer();
      param.iter_newton = (unsigned int) cl.read_integer();
      }
    else if(arg == "-C")
      {
      param.constrained_mu_init = cl.read_double();
      param.constrained_mu_mult = cl.read_double();
      }
    else if(arg == "-f")
      {
      param.use_float = true;
      }
    else if(arg == "-p")
      {
      param.arrInitialMomentum = cl.read_string();
      }
    else if(arg == "-L")
      {
      param.arrAttachmentLabelPosteriors = cl.read_string();
      }
    else if(arg == "-t")
      {
      param.n_threads = cl.read_integer();
      }
    else if(arg == "-D")
      {
      param.n_deriv_check = cl.read_integer();
      }
    else if(arg == "-a")
      {
      std::string mode = cl.read_string();
      cout << mode << endl;
      if(mode == "L")
        param.attach = ShootingParameters::Euclidean;
      else if(mode == "C")
        param.attach = ShootingParameters::Current;
      else if(mode == "V")
        param.attach = ShootingParameters::Varifold;
      else
        {
        std::cout << "Unknown attachment type " << mode << std::endl;
        return -1;
        }
      }
    else if(arg == "-S")
      {
      param.currents_sigma = cl.read_double();
      }
    else if(arg == "-test-currents")
      {
      param.test_currents_attachment = true;
      }
    else if(arg == "-h")
      {
      return usage();
      }
    else
      {
      cerr << "Unknown option " << arg << endl;
      return -1;
      }
    }

  // Check parameters
  check(param.sigma > 0, "Missing or negative sigma parameter");
  check(param.attach == ShootingParameters::Euclidean || param.currents_sigma > 0,
        "Missing sigma parameter for current/varifold metric");
  check(param.N > 0 && param.N < 10000, "Incorrect N parameter");
  check(param.dim >= 2 && param.dim <= 3, "Incorrect N parameter");
  check(param.fnTemplate.length(), "Missing template filename");
  check(param.fnTarget.length(), "Missing target filename");
  check(param.fnOutput.length(), "Missing output filename");

  // Set the number of threads if not specified
  if(param.n_threads == 0)
    param.n_threads = std::thread::hardware_concurrency();

  // Specialize by dimension
  if(param.use_float)
    {
    if(param.dim == 2)
      return PointSetShootingProblem<float, 2>::minimize(param);
    else
      return PointSetShootingProblem<float, 3>::minimize(param);
    }
  else
    {
    if(param.dim == 2)
      return PointSetShootingProblem<double, 2>::minimize(param);
    else
      return PointSetShootingProblem<double, 3>::minimize(param);
    }


}
