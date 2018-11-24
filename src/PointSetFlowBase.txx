#include "PointSetFlowBase.h"
#include "ctpl_stl.h"

template <class TFloat, unsigned int VDim, class TThreadData>
PointSetFlowBase<TFloat, VDim, TThreadData>
::PointSetFlowBase(
    const Matrix &q0, TFloat sigma, unsigned int N)
{
  // Copy parameters
  this->q0 = q0;
  this->sigma = sigma;
  this->N = N;
  this->k = q0.rows();
  this->dt = 1.0 / (N-1);

  // Set up thread data
  this->SetupMultiThreaded();
}


template <class TFloat, unsigned int VDim, class TThreadData>
PointSetFlowBase<TFloat, VDim, TThreadData>
::~PointSetFlowBase()
{
  delete thread_pool;
}

template <class TFloat, unsigned int VDim, class TThreadData>
void
PointSetFlowBase<TFloat, VDim, TThreadData>
::SetupMultiThreaded()
{
  unsigned int n_threads = std::thread::hardware_concurrency();

  // Initialize the thread data
  td.resize(n_threads, TThreadData(this->N, this->k));

  // Create the thread pool
  thread_pool = new ctpl::thread_pool(n_threads);
 
  // Assign lines in pairs, one at the top of the symmetric matrix K and
  // one at the bottom of K. The loop below will not assign the middle
  // line when there is an odd number of points (e.g., line 7 when there are 15)
  for(int i = 0; i < k/2; i++)
    {
    int i_thread = i % n_threads;
    td[i_thread].rows.push_back(i);
    td[i_thread].rows.push_back((k-1) - i);
    }

  // Handle the middle line for odd number of vertices
  if(k % 2 == 1)
    td[(k / 2) % n_threads].rows.push_back(k/2);
}


