#ifndef __PointSetFlowBase_h_
#define __PointSetFlowBase_h_

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>
#include <vector>

namespace ctpl { class thread_pool; }

struct PointSetFlowThreadDataBase
{
  std::vector<unsigned int> rows;
};

template <class TFloat, unsigned int VDim, class TThreadData>
class PointSetFlowBase
{
public:
  typedef vnl_matrix<TFloat> Matrix;
  typedef vnl_vector<TFloat> Vector;
  typedef vnl_vector_fixed<TFloat, VDim> VecD;
  typedef std::vector<Matrix> MatrixArray;

  
  /**
   * Constructor - set basic parameters of the system and allocate
   * all of the necessary matrices
   * 
   * q0      : N x D vector of template landmark positions
   * sigma   : standard deviation of the Gaussian kernel 
   * N       : number of timesteps for the ODE
   */
  PointSetFlowBase(
    const Matrix &q0, 
    TFloat sigma,
    unsigned int N);

  ~PointSetFlowBase();

  /** Get the number of time steps */
  unsigned int GetN() const { return N; }

  /** Get the point trajectories */
  const Matrix &GetQt(unsigned int t) const { return Qt[t]; }

  /** Get the point velocities */
  const Matrix &GetVt(unsigned int t) const { return Vt[t]; }

  /** Get the step size */
  const TFloat GetDeltaT() const { return dt; }

protected:

  // Multi-threading setup
  void SetupMultiThreaded();

  // Initial ladnmark coordinates - fixed for duration
  Matrix q0;

  // Standard deviation of Gaussian kernel; time step
  TFloat sigma, dt;

  // Number of timesteps for integration; number of points
  unsigned int N, k;

  // Streamlines - paths of the landmarks over time
  std::vector<Matrix> Qt;

  // Streamline velocities
  std::vector<Matrix> Vt;

  // Data associated with each thread
  std::vector<TThreadData> td;

  // A thread pool to handle parallel execution
  ctpl::thread_pool *thread_pool;

};


#endif

