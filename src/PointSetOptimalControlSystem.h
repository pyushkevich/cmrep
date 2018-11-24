#ifndef __PointSetOptimalControlSystem_h_
#define __PointSetOptimalControlSystem_h_

#include "PointSetFlowBase.h"

/**
 * The thread data for the flow system
 */
template <class TFloat, unsigned int VDim>
struct PointSetOptimalControlSystemThreadData : public PointSetFlowThreadDataBase
{
  // List of rows handled by this thread
  TFloat KE;
  vnl_vector<TFloat> d_q__d_t[VDim];
  vnl_vector<TFloat> alpha_U[VDim], alpha_Q[VDim];

  // Initialize the thread data
  PointSetOptimalControlSystemThreadData(unsigned int n, unsigned int k)
    {
    for(int a = 0; a < VDim; a++)
      {
      d_q__d_t[a] = vnl_vector<TFloat>(k, 0.0);
      alpha_U[a] = vnl_vector<TFloat>(k, 0.0);
      alpha_Q[a] = vnl_vector<TFloat>(k, 0.0);
      }
    }
};

template <class TFloat, unsigned int VDim>
class PointSetOptimalControlSystem : 
  public PointSetFlowBase<TFloat, VDim, PointSetOptimalControlSystemThreadData<TFloat, VDim> >
{
public:

  typedef PointSetOptimalControlSystemThreadData<TFloat, VDim> ThreadData;
  typedef PointSetFlowBase<TFloat, VDim, ThreadData> Superclass;
  typedef typename Superclass::Matrix Matrix;
  typedef typename Superclass::Vector Vector;
  typedef typename Superclass::VecD VecD;
  typedef typename Superclass::MatrixArray MatrixArray;

  /**
   * Constructor - set basic parameters of the system and allocate
   * all of the necessary matrices
   * 
   * q0      : N x D vector of template landmark positions
   * sigma   : standard deviation of the Gaussian kernel 
   * N       : number of timesteps for the ODE
   */
  PointSetOptimalControlSystem(
    const Matrix &q0, 
    TFloat sigma,
    unsigned int N);

  /**
   * Compute the kinetic energy and the flow velocities at a given timepoint
   */
  TFloat ComputeEnergyAndVelocity(const Matrix &q, const Matrix &u);

  /**
   * Perform forward flow with control u, returning the total kinetic energy
   * of the flow. The endpoint and intermediate timepoints can be queried
   * using GetQt()
   */
  TFloat Flow(const MatrixArray &u);

  /**
   * Perform backward flow in order to compute the gradient of some function
   * f with respect to the control u. The input is the array of partial derivatives
   * of the function f with respect to the path q(t).
   *
   * The function also includes the partial derivatives of the kinetic energy
   * with respect to u, with weight w_kinetic 
   */
  void FlowBackward(const MatrixArray &u, const MatrixArray &d_f__d_qt, 
                    TFloat w_kinetic, MatrixArray &d_f__d_u);

protected:

  // Step of backpropagation
  void PropagateAlphaBackwards(
    const Matrix &q, const Matrix &u, 
    const Vector alpha[], Vector alpha_Q[], Vector alpha_U[]);

  // The current velocity vector
  Vector d_q__d_t[VDim];

  void ComputeEnergyAndVelocityThreadedWorker(
    const Matrix &q, const Matrix &u, ThreadData *tdi);

  void PropagateAlphaBackwardsThreadedWorker(
    const Matrix &q, const Matrix &u, 
    const Vector alpha[], ThreadData *tdi);
};



#endif // __PointSetOptimalControlSystem_h_
