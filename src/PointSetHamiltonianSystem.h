#ifndef __PointSetHamiltionianSystem_h_
#define __PointSetHamiltionianSystem_h_

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <vector>

template <class TFloat, unsigned int VDim>
class PointSetHamiltonianSystem
{
public:

  typedef vnl_matrix<TFloat> Matrix;
  typedef vnl_vector<TFloat> Vector;

  typedef vnl_vector_fixed<TFloat, VDim> VecD;

  /**
   * Constructor - set basic parameters of the system and allocate
   * all of the necessary matrices
   * 
   * q0      : N x D vector of template landmark positions
   * sigma   : standard deviation of the Gaussian kernel 
   * N       : number of timesteps for the ODE
   */
  PointSetHamiltonianSystem(
    const Matrix &q0, 
    TFloat sigma,
    unsigned int N);

  /**
   * Compute the Hamiltonian and its derivatives for given p/q. The
   * computation is stored internally in variables Hp, Hq, Hqq, Hqp, Hpp
   */
  TFloat ComputeHamiltonianJet(const Matrix &q, const Matrix &p, bool flag_hessian);

  /**
   * Flow the Hamiltonian system with initial momentum p0 without gradient 
   * computation. Returns the kinetic energy (Hamiltonian value that should 
   * be preserved over the time evolution)
   */
  TFloat FlowHamiltonian(const Matrix &p0, Matrix &q, Matrix &p);

  /**
   * Flow the Hamiltonian system with gradient computation. The gradient is 
   * strored as a VDim x VDim array of k x k matrices. This is a pretty expensive
   * operation because of many matrix multiplications that are involved
   */
  TFloat FlowHamiltonianWithGradient(
    const Matrix &p0, Matrix &q, Matrix &p,
    Matrix grad_q[VDim][VDim], Matrix grad_p[VDim][VDim]);

  /** 
   * Computes the expression alpha' * Q1 + beta' * P1, where alpha and beta are 
   * vectors, and Q1 and P1 are D_p0(q1) and D_p0(p1), respectively.
   *
   * This can be used to efficiently compute the gradient of any function f(q1,p1)
   * with respect to the initial momentum, since 
   *
   *    D_p0(f) = Dq_f' * D_p0(q1) + Dp_f' * D_p0(p1)
   *
   * The method flows alpha and beta back in time, avoiding expensive matrix-matrix
   * multiplication. This method requires FlowHamiltonian to have been run already
   */ 
  void FlowGradientBackward(
    const Vector alpha[VDim], const Vector beta[VDim], Vector result[VDim]);

  const Vector &GetHp(unsigned int d) const { return Hp[d]; }
  const Vector &GetHq(unsigned int d) const { return Hq[d]; }
  const Matrix &GetHqq(unsigned int a, unsigned int b) const { return Hqq[a][b]; }
  const Matrix &GetHqp(unsigned int a, unsigned int b) const { return Hqp[a][b]; }
  const Matrix &GetHpp(unsigned int a, unsigned int b) const { return Hpp[a][b]; }

protected:

  // Initial ladnmark coordinates - fixed for duration
  Matrix q0;

  // Standard deviation of Gaussian kernel; time step
  TFloat sigma, dt;

  // Number of timesteps for integration; number of points
  unsigned int N, k;

  // Gradient of the Hamiltonian components: Hq and Hp
  Vector Hp[VDim], Hq[VDim];

  // Hessian of the Hamiltonian components: Hqq, Hqp, Hpp
  // matrices Hqq and Hpp are symmetric
  Matrix Hqq[VDim][VDim], Hqp[VDim][VDim], Hpp[VDim][VDim];

  // Streamlines - paths of the landmarks over time
  std::vector<Matrix> Qt, Pt;
};



#endif // __PointSetHamiltionianSystem_h_
