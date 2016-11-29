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
  cout << "  -i iter                    : max iterations for optimization" << endl;
  cout << "  -r fraction                : randomly downsample mesh by factor (e.g. 0.01)" << endl;
  cout << "  -q vector                  : apply quadric clustering to mesh with specified" <<  endl;
  cout << "                               number of bins per dimension (e.g. 40x40x30)" << endl;
  cout << "  -a <A|G>                   : algorithm to use: A: Allassonniere; G: GradDescent (deflt)" << endl;
  cout << "  -O filepattern             : pattern for saving traced landmark paths (e.g., path%04d.vtk)" << endl;
  return -1;
}

void check(bool condition, char *format,...)
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
  enum Algorithm { Allassonniere, GradDescent };
  string fnTemplate, fnTarget;
  string fnOutput;
  string fnOutputPaths;
  double sigma;
  double lambda;
  double downsample;
  unsigned int dim;
  unsigned int N;
  unsigned int iter;
  Algorithm alg;

  std::vector<int> qcdiv;

  ShootingParameters() :
    N(100), dim(3), sigma(0.0), lambda(0.0), iter(120), downsample(1),
    alg(GradDescent) {}
};

template <unsigned int VDim>
class PointSetShootingCostFunction : public vnl_cost_function
{
public:
  typedef PointSetHamiltonianSystem<double, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;

  PointSetShootingCostFunction(
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
      grad_f[a].set_size(k);
      }
    }

  Vector wide_to_tall(const Vector p[VDim])
    {
    Vector v(k * VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        v[pos++] = p[a](i);
    return v;
    }

  Vector wide_to_tall(const Matrix &p)
    {
    Vector v(k * VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        v[pos++] = p(i,a);
    return v;
    }

  Matrix tall_to_wide(const Vector &v)
    {
    Matrix p(k,VDim);
    int pos = 0;
    for(unsigned int a = 0; a < VDim; a++)
      for(unsigned int i = 0; i < k; i++)
        p(i,a) = v[pos++];
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

      // Nice thing is that the last call to ComputeHamiltonianJet in the 
      // above method is at (q0, p0), so the value of Hp is already valid

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

template <unsigned int VDim>
class PointSetShootingProblem
{
public:
  typedef PointSetHamiltonianSystem<double, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;

  // Minimize using the transversality principle
  static void minimize_Allassonniere(const ShootingParameters &param, 
    const Matrix &q0, const Matrix &qT, Matrix &p0);

  // Minimize using gradient descent 
  static void minimize_gradient(const ShootingParameters &param, 
    const Matrix &q0, const Matrix &qT, Matrix &p0);

  static int minimize(const ShootingParameters &param);

private:
};

template <unsigned int VDim>
void
PointSetShootingProblem<VDim>
::minimize_Allassonniere(const ShootingParameters &param, 
  const Matrix &q0, const Matrix &qT, Matrix &p0)
{
  unsigned int k = q0.rows();

  // Create the hamiltonian system
  HSystem hsys(q0, param.sigma, param.N);

  // Where to store the results of the flow
  Matrix q1(k,VDim), p1(k,VDim), grad_q[VDim][VDim], grad_p[VDim][VDim];
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
  for(unsigned int iter = 0; iter < param.iter; iter++)
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
    vnl_svd<double> svd(DG, 1.0);

    // Compute inv(DG) * G
    Vector del_p0 = - svd.solve(G);
    
    // Print the current state
    printf("Iter %4d   H=%8.6f   l*Dsq=%8.6f   E=%8.6f   |G|=%8.6f\n",
      iter, H, param.lambda * dsq, H + param.lambda * dsq, G.two_norm());

    // Update the solution
    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int i = 0; i < k; i++)
        {
        p0(i,a) += 0.1 * del_p0(a * k + i);
        }
      }
    }
}

template <unsigned int VDim>
void
PointSetShootingProblem<VDim>
::minimize_gradient(const ShootingParameters &param, 
  const Matrix &q0, const Matrix &qT, Matrix &p0)
{
  unsigned int k = q0.rows();

  // Create the minimization problem
  PointSetShootingCostFunction<VDim> cost_fn(param, q0, qT);

  // Create initial/final solution
  p0 = (qT - q0) / param.N;
  Vector x = cost_fn.wide_to_tall(p0);

  // Solve the minimization problem 
  vnl_lbfgs optimizer(cost_fn);
  optimizer.set_f_tolerance(1e-9);
  optimizer.set_x_tolerance(1e-4);
  optimizer.set_g_tolerance(1e-6);
  optimizer.set_trace(true);
  optimizer.set_max_function_evals(param.iter);

  optimizer.minimize(x);

  // Take the optimal solution
  p0 = cost_fn.tall_to_wide(x);
}

template <unsigned int VDim>
int
PointSetShootingProblem<VDim>
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

  // Select the method to run
  if(param.alg == ShootingParameters::Allassonniere)
    minimize_Allassonniere(param, q0, qT, p0);
  else
    minimize_gradient(param, q0, qT, p0);

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
        double qi[VDim], vi[VDim];
        pTemplate->GetPoint(i, qi);

        // Interpolate the velocity at each mesh point
        hsys.InterpolateVelocity(t-1, qi, vi);

        // Update the position using Euler's method
        for(unsigned int a = 0; a < VDim; a++)
          qi[a] += dt * vi[a];

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
      param.iter = (unsigned int) cl.read_integer();
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
      else
        param.alg = ShootingParameters::GradDescent;
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
  if(param.dim == 2)
    return PointSetShootingProblem<2>::minimize(param);
  else
    return PointSetShootingProblem<3>::minimize(param);


}
