#ifndef __PointSetHamiltionianSystem_h_
#define __PointSetHamiltionianSystem_h_

#include <vnl_matrix.h>
#include <vnl_vector.h>

template <class TFloat, class VDim>
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
  double ComputeHamiltonianJet(const Matrix &q, const Matrix &p, bool flag_hessian);

  /**
   * Flow the Hamiltonian system with initial momentum p0 without gradient 
   * computation. Returns the kinetic energy (Hamiltonian value that should 
   * be preserved over the time evolution)
   */
  double FlowHamiltonian(const Matrix &p0, Matrix &q, Matrix &p);

  /**
   * Flow the Hamiltonian system with gradient computation. The gradient is 
   * strored as a VDim x VDim array of k x k matrices
   */
  double FlowHamiltonianWithGradient(
    const Matrix &p0, Matrix &q, Matrix &p,
    Matrix grad_q[VDim][VDim], Matrix grad_p[VDim][VDim]);

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

  // Flow results: Q and P integrated over time
  Matrix Qt[VDim], Q
}



#endif // __PointSetHamiltionianSystem_h_
