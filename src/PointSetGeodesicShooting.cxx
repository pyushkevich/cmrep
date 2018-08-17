#include "PointSetHamiltonianSystem.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/vnl_cost_function.h"
#include "vnl/algo/vnl_lbfgs.h"

#include <iostream>
#include <algorithm>
#include "util/ReadWriteVTK.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkQuadricClustering.h"
#include "vtkSmartPointer.h"

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
  cout << "  -i iter_grad iter_newt     : max iterations for optimization for gradient descent and newton's" << endl;
  cout << "  -r fraction                : randomly downsample mesh by factor (e.g. 0.01)" << endl;
  cout << "  -q vector                  : apply quadric clustering to mesh with specified" <<  endl;
  cout << "                               number of bins per dimension (e.g. 40x40x30)" << endl;
//  cout << "  -a <A|G>                   : algorithm to use: A: Allassonniere; G: GradDescent (deflt)" << endl;
  cout << "  -O filepattern             : pattern for saving traced landmark paths (e.g., path%04d.vtk)" << endl;
  cout << "  -f                         : use single-precision float (off by deflt)" << endl;
  cout << "  -C mu0 mu_mult             : test constrained optimization (not for general use)" << endl;
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
  string fnTemplate, fnTarget;
  string fnOutput;
  string fnOutputPaths;
  double sigma;
  double lambda;
  double downsample;
  unsigned int dim;
  unsigned int N;
  unsigned int iter_grad, iter_newton;
  Algorithm alg;
  bool use_float;

  // For constrained optimization - just exprimental
  double constrained_mu_init, constrained_mu_mult;

  std::vector<int> qcdiv;

  ShootingParameters() :
    N(100), dim(3), sigma(0.0), lambda(0.0), iter_grad(20), iter_newton(20), downsample(1),
    use_float(false), alg(GradDescent), constrained_mu_init(0.0), constrained_mu_mult(0.0) 
  {}
};

template <class TFloat, unsigned int VDim>
class PointSetShootingCostFunction : public vnl_cost_function
{
public:
  typedef PointSetHamiltonianSystem<TFloat, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;

  // Separate type because vnl optimizer is double-only
  typedef vnl_vector<double> DVector;

  PointSetShootingCostFunction(
    const ShootingParameters &param, const Matrix &q0, const Matrix &qT)
    : vnl_cost_function(q0.rows() * VDim), hsys(q0, param.sigma, param.N)
    {
    this->p0 = (qT - q0) / param.N;
    this->q0 = q0;
    this->qT = qT;
    this->param = param;
    this->k = q0.rows();
    this->p1.set_size(k,VDim);
    this->q1.set_size(k,VDim);

    for(unsigned int a = 0; a < VDim; a++)
      {
      alpha[a].set_size(k);
      beta[a].set_size(k); beta[a].fill(0.0);
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

    // Compute the landmark errors
    double fnorm_sq = 0.0;
    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int i = 0; i < k; i++)
        {
        alpha[a](i) = q1(i,a) - qT(i,a);
        fnorm_sq += alpha[a](i) * alpha[a](i);
        }
      }

    // Compute the landmark part of the objective
    double Edist = 0.5 * fnorm_sq;

    if(f)
      *f = H + param.lambda * Edist;

    if(g)
      {
      // Multiply gradient of f. wrt q1 (alpha) by Jacobian of q1 wrt p0
      hsys.FlowGradientBackward(alpha, beta, grad_f);

      // Recompute Hq/Hp at initial timepoint
      hsys.ComputeHamiltonianJet(q0, p0, false);

      // Complete gradient computation
      for(unsigned int a = 0; a < VDim; a++)
        {
        // Combine the gradient terms
        grad_f[a] = grad_f[a] * param.lambda + hsys.GetHp(a);
        }

      // Pack the gradient into the output vector
      *g = wide_to_tall(grad_f);
      }
    }



protected:
  HSystem hsys;
  ShootingParameters param;
  Matrix qT, p0, q0, p1, q1;
  Vector alpha[VDim], beta[VDim], grad_f[VDim];
  unsigned int k;

};

/**
 * This is just a simple test of a cost function to be used with augmented lagrangian methods
 */
template <class TFloat, unsigned int VDim>
class ExampleConstrainedPointSetShootingCostFunction : public PointSetShootingCostFunction<TFloat, VDim>
{
public:

  typedef PointSetShootingCostFunction<TFloat,VDim> Superclass;
  typedef typename Superclass::HSystem HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;

  // Separate type because vnl optimizer is double-only
  typedef vnl_vector<double> DVector;

  ExampleConstrainedPointSetShootingCostFunction(
    const ShootingParameters &param, const Matrix &q0, const Matrix &qT)
    : Superclass(param, q0, qT)
    {
    // Compute the vectors with which the flows have to be orthogonal
    norm.set_size(q0.rows(), q0.columns());
    for(int i = 0; i < q0.rows(); i++)
      {
      Vector del = qT.get_row(i) - q0.get_row(i);
      double len = del.magnitude();
      norm(i,0) = -del[1] / len;
      norm(i,1) = del[0] / len;
      }
    }

  void set_lambdas(Matrix &m)
    {
    this->lambda = m;
    }

  const Matrix & lambdas() const { return this->lambda; }

  void set_mu(double mu)
    {
    this->mu = mu;
    }

  virtual void compute(vnl_vector<double> const& x, double *f, vnl_vector<double>* g)
    {
    // Initialize the p0-vector
    this->p0 = this->tall_to_wide(x);

    // Perform flow
    double H = this->hsys.FlowHamiltonian(this->p0, this->q1, this->p1);

    // Compute the landmark errors
    double fnorm_sq = 0.0;
    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int i = 0; i < this->k; i++)
        {
        this->alpha[a](i) = this->q1(i,a) - this->qT(i,a);
        fnorm_sq += this->alpha[a](i) * this->alpha[a](i);
        }
      }

    // Compute the landmark part of the objective
    double Edist = 0.5 * fnorm_sq;

    if(f)
      {
      // Compute the objective part
      *f = H + this->param.lambda * Edist;

      // Compute the constraints
      int i_lam = 0;
      for(int t = 0; t < this->hsys.GetN(); t++)
        {
        const Matrix &qt = this->hsys.GetQt(t);
        for(int i = 0; i < this->q0.rows(); i++)
          {
          // C(i) is the dot product of norm vector with current time vector
          Vector del = qt.get_row(i) - this->q0.get_row(i);
          double Cti = dot_product(del, norm.get_row(i));
          *f -= lambda(t, i) * Cti;
          *f += (mu / 2) * Cti * Cti;
          }
        }
      }

    if(g)
      {
      // Create the gradient array to pass to backward flow
      std::vector<Matrix> d_obj__d_qt(this->hsys.GetN());

      // Fill out the gradients due to the constraint parts of the augmented Lagrangian
      for(int t = 0; t < this->hsys.GetN(); t++)
        {
        d_obj__d_qt[t].set_size(this->q0.rows(), VDim);

        const Matrix &qt = this->hsys.GetQt(t);
        for(int i = 0; i < this->q0.rows(); i++)
          {
          Vector del = qt.get_row(i) - this->q0.get_row(i);
          double Cti = dot_product(del, norm.get_row(i));
          double term1 = ( -lambda(t, i) + mu * Cti);
          d_obj__d_qt[t].set_row(i, norm.get_row(i) * term1);
          }
        }

      // Add the parts due to the objective function
      for(int i = 0; i < this->q0.rows(); i++)
        {
        for(int a = 0; a < this->q0.columns(); a++)
          {
          d_obj__d_qt[this->hsys.GetN() - 1](i,a) += this->alpha[a](i) * this->param.lambda;
          }
        }

      // Multiply gradient of f. wrt q1 (alpha) by Jacobian of q1 wrt p0
      this->hsys.FlowTimeVaryingGradientsBackward(d_obj__d_qt, this->grad_f);

      // Recompute Hq/Hp at initial timepoint
      this->hsys.ComputeHamiltonianJet(this->q0, this->p0, false);

      // Complete gradient computation
      for(unsigned int a = 0; a < VDim; a++)
        {
        // Combine the gradient terms
        this->grad_f[a] = this->grad_f[a] + this->hsys.GetHp(a);
        }

      // Pack the gradient into the output vector
      *g = this->wide_to_tall(this->grad_f);
      }
    }

  void update_lambdas(vnl_vector<double> const& x)
    {
    // Initialize the p0-vector
    this->p0 = this->tall_to_wide(x);

    // Perform flow
    double H = this->hsys.FlowHamiltonian(this->p0, this->q1, this->p1);

    for(int t = 0; t < this->hsys.GetN(); t++)
      {
      const Matrix &qt = this->hsys.GetQt(t);
      for(int i = 0; i < this->q0.rows(); i++)
        {
        // C(i) is the dot product of norm vector with current time vector
        Vector del = qt.get_row(i) - this->q0.get_row(i);
        double Cti = dot_product(del, norm.get_row(i));
        lambda(t, i) -= mu * Cti;
        }
      }
    }

protected:

  Matrix lambda;
  Matrix norm;
  double mu;

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

  // Minimize using the transversality principle
  static void minimize_Allassonniere(const ShootingParameters &param, 
    const Matrix &q0, const Matrix &qT, Matrix &p0);

  static void minimize_QuasiAllassonniere(const ShootingParameters &param,
    const Matrix &q0, const Matrix &qT, Matrix &p0);

  // Minimize using gradient descent
  static void minimize_gradient(const ShootingParameters &param, 
    const Matrix &q0, const Matrix &qT, Matrix &p0);

  // Minimize constrained problem using Augmented Lagrangian method
  static void minimize_constrained(const ShootingParameters &param, 
    const Matrix &q0, const Matrix &qT, Matrix &p0);

  static int minimize(const ShootingParameters &param);

private:
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
  HSystem hsys(q0, param.sigma, param.N);

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
::minimize_gradient(const ShootingParameters &param, 
  const Matrix &q0, const Matrix &qT, Matrix &p0)
{
  unsigned int k = q0.rows();

  // Create the minimization problem
  typedef PointSetShootingCostFunction<TFloat, VDim> CostFn;
  CostFn cost_fn(param, q0, qT);

  // Create initial/final solution
  p0 = (qT - q0) / param.N;
  typename CostFn::DVector x = cost_fn.wide_to_tall(p0);

  // Uncomment this code to test derivative computation
  /*
  TFloat eps = 1e-6;
  typename CostFn::DVector test_grad(x.size());
  double f_test;
  cost_fn.compute(x, &f_test, &test_grad);
  for(int i = 0; i < std::min((int) p0.size(), 10); i++)
    {
    typename CostFn::DVector xtest = x;
    double f1, f2;
    xtest[i] = x[i] - eps;
    cost_fn.compute(xtest, &f1, NULL);

    xtest[i] = x[i] + eps;
    cost_fn.compute(xtest, &f2, NULL);

    printf("i = %03d,  AG = %8.4f,  NG = %8.4f\n", i, test_grad[i], (f2 - f1) / (2 * eps));
    }
  */

  // Solve the minimization problem 

  vnl_lbfgsb optimizer(cost_fn);
  optimizer.set_f_tolerance(1e-9);
  optimizer.set_x_tolerance(1e-4);
  optimizer.set_g_tolerance(1e-6);
  optimizer.set_trace(true);
  optimizer.set_max_function_evals(param.iter_grad);

  // vnl_conjugate_gradient optimizer(cost_fn);
  // optimizer.set_trace(true);
  optimizer.minimize(x);

  // Take the optimal solution
  p0 = cost_fn.tall_to_wide(x);
}



template <class TFloat, unsigned int VDim>
void
PointSetShootingProblem<TFloat, VDim>
::minimize_constrained(const ShootingParameters &param, 
  const Matrix &q0, const Matrix &qT, Matrix &p0)
{
  unsigned int k = q0.rows();

  // Create the minimization problem
  typedef ExampleConstrainedPointSetShootingCostFunction<TFloat, VDim> CostFn;
  CostFn cost_fn(param, q0, qT);

  // Create initial/final solution
  p0 = (qT - q0) / param.N;

  // Initialize mu
  double mu = param.constrained_mu_init;
  cost_fn.set_mu(mu);

  // Initialize the lambdas
  typename CostFn::Matrix lambdas(param.N, k);
  lambdas.fill(0.0);
  cost_fn.set_lambdas(lambdas);

  // Do this for some number of iterations
  for(int it = 0; it < 10; it++)
    {
    typename CostFn::DVector x = cost_fn.wide_to_tall(p0);

    // Update the parameters
    printf("Augmented Lagrangian Phase %d (mu = %8.4f, max-lambda = %8.4f)\n",
      it, mu, cost_fn.lambdas().array_inf_norm());

    // Uncomment this code to test derivative computation
    TFloat eps = 1e-6;
    typename CostFn::DVector test_grad(x.size());
    double f_test;
    cost_fn.compute(x, &f_test, &test_grad);
    for(int i = 0; i < std::min((int) p0.size(), 10); i++)
      {
      typename CostFn::DVector xtest = x;
      double f1, f2;
      xtest[i] = x[i] - eps;
      cost_fn.compute(xtest, &f1, NULL);

      xtest[i] = x[i] + eps;
      cost_fn.compute(xtest, &f2, NULL);

      printf("i = %03d,  AG = %12.8f,  NG = %12.8f\n", i, test_grad[i], (f2 - f1) / (2 * eps));
      }

    // Create an optimization
    vnl_lbfgsb optimizer(cost_fn);
    optimizer.set_f_tolerance(1e-9);
    optimizer.set_x_tolerance(1e-4);
    optimizer.set_g_tolerance(1e-6);
    optimizer.set_trace(true);
    optimizer.set_max_function_evals(param.iter_grad / 10);

    // vnl_conjugate_gradient optimizer(cost_fn);
    // optimizer.set_trace(true);
    optimizer.minimize(x);

    // Take the optimal solution from this round
    p0 = cost_fn.tall_to_wide(x);

    // Update the lambdas
    cost_fn.update_lambdas(x);

    // Update the mu
    mu *= param.constrained_mu_mult;
    cost_fn.set_mu(mu);
    }
}





template <class TFloat, unsigned int VDim>
int
PointSetShootingProblem<TFloat, VDim>
::minimize(const ShootingParameters &param)
{
  // Read the datasets
  vtkPolyData *pTemplate = ReadVTKData(param.fnTemplate);
  vtkPolyData *pTarget = ReadVTKData(param.fnTarget);

  // Get the number of vertices and dimensionality
  check(pTemplate->GetNumberOfPoints() == pTarget->GetNumberOfPoints(), "Meshes don't match");

  unsigned int np = pTemplate->GetNumberOfPoints();
  unsigned int k = np;

  // Mapping of matrix index to vertex - start with identity mapping
  std::vector<unsigned int> index(np);
  for(int i = 0; i < np; i++)
    index[i] = i;

  // Downsample meshes if needed
  if(param.downsample < 1.0)
    {
    // Downsample by random shuffle - not very efficient, just lazy
    k = (unsigned int)(param.downsample * np);
    std::random_shuffle(index.begin(), index.end());
    }
  else if(param.qcdiv.size())
    {
    // Store the index of each input point in VTK mesh
    vtkSmartPointer<vtkIntArray> arr_idx = vtkIntArray::New();
    arr_idx->SetNumberOfComponents(1);
    arr_idx->SetNumberOfTuples(np);
    arr_idx->SetName("Index");
    for(int i = 0; i < np; i++)
      arr_idx->SetTuple1(i, i);
    pTemplate->GetPointData()->AddArray(arr_idx);

    // Perform quadric clustering
    vtkSmartPointer<vtkQuadricClustering> qc = vtkQuadricClustering::New();
    qc->SetInputData(pTemplate);
    qc->SetUseInputPoints(1);

    // Set the divisions array - it must be 3
    int div[3];
    for(int a = 0; a < 3; a++)
      div[a] = (a < param.qcdiv.size()) ? param.qcdiv[a] : 1;
    qc->SetNumberOfDivisions(div);

    // Run the filter
    qc->Update();

    // Generate the index
    vtkPolyData *qc_result = qc->GetOutput();
    vtkDataArray *qc_index = qc_result->GetPointData()->GetArray("Index");
    k = qc_result->GetNumberOfPoints();
    for(int i = 0; i < k; i++)
      index[i] = qc_index->GetTuple1(i);
    }

  // Report number of actual vertices being used
  cout << "Performing geodesic shooting with " << k << " landmarks" << endl;

  // Landmarks and initial momentum
  Matrix q0(k,VDim), qT(k,VDim), p0(k,VDim);

  // Initialize the meshes
  for(unsigned int i = 0; i < k; i++)
    {
    for(unsigned int a = 0; a < VDim; a++)
      {
      q0(i,a) = pTemplate->GetPoint(index[i])[a];
      qT(i,a) = pTarget->GetPoint(index[i])[a];
      p0(i,a) = (qT(i,a) - q0(i,a)) / param.N;
      }
    }

  // Run some iterations of gradient descent
  if(param.iter_grad > 0)
    {
    if(param.constrained_mu_init > 0.0)
      minimize_constrained(param, q0, qT, p0);

    else
      minimize_gradient(param, q0, qT, p0);
    }

  if(param.iter_newton > 0)
    {
    minimize_Allassonniere(param, q0, qT, p0);
    }

  // Select the method to run
  /*
  if(param.alg == ShootingParameters::Allassonniere)
    minimize_Allassonniere(param, q0, qT, p0);
  else if(param.alg == ShootingParameters::QuasiAllassonniere)
    minimize_QuasiAllassonniere(param, q0, qT, p0);
  else
    minimize_gradient(param, q0, qT, p0);
    */

  // Genererate the momentum map
  vtkDoubleArray *arr_p = vtkDoubleArray::New();
  arr_p->SetNumberOfComponents(VDim);
  arr_p->SetNumberOfTuples(np);
  arr_p->SetName("InitialMomentum");
  for(unsigned int a = 0; a < VDim; a++)
    arr_p->FillComponent(a, 0);

  for(unsigned int a = 0; a < VDim; a++)
    {
    for(unsigned int i = 0; i < k; i++)
      {
      arr_p->SetComponent(index[i],a,p0(i,a));
      }
    }

  pTemplate->GetPointData()->AddArray(arr_p);
  WriteVTKData(pTemplate, param.fnOutput);

  // If saving paths requested
  if(param.fnOutputPaths.size())
    {
    // Create and flow a system
    HSystem hsys(q0, param.sigma, param.N);
    Matrix q1, p1;
    hsys.FlowHamiltonian(p0, q1, p1);

    // Apply the flow to the points in the rest of the mesh
    vtkDoubleArray *arr_v = vtkDoubleArray::New();
    arr_v->SetNumberOfComponents(VDim);
    arr_v->SetNumberOfTuples(np);
    arr_v->SetName("Velocity");
    pTemplate->GetPointData()->AddArray(arr_v);

    // Apply Euler method to the mesh points
    double dt = hsys.GetDeltaT();
    for(unsigned int t = 1; t < param.N; t++)
      {
      for(unsigned int i = 0; i < np; i++)
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
    else if(arg == "-r")
      {
      param.downsample = cl.read_double();
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
    else if(arg == "-q")
      {
      param.qcdiv = cl.read_int_vector();
      }
    else if (arg == "-a")
      {
      string alg = cl.read_string();
      if(alg == "a" || alg == "A" || alg == "Allassonniere")
        param.alg = ShootingParameters::Allassonniere;
      else if(alg == "q" || alg == "Q" || alg == "QuasiAllassonniere")
        param.alg = ShootingParameters::QuasiAllassonniere;
      else
        param.alg = ShootingParameters::GradDescent;
      }
    else if(arg == "-f")
      {
      param.use_float = true;
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
  check(param.N > 0 && param.N < 10000, "Incorrect N parameter");
  check(param.dim >= 2 && param.dim <= 3, "Incorrect N parameter");
  check(param.fnTemplate.length(), "Missing template filename");
  check(param.fnTarget.length(), "Missing target filename");
  check(param.fnOutput.length(), "Missing output filename");

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
