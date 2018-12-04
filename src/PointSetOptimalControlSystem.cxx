#include "PointSetOptimalControlSystem.h"
#include <vnl/vnl_fastops.h>
#include <iostream>

template <class TFloat, unsigned int VDim>
PointSetOptimalControlSystem<TFloat, VDim>
::PointSetOptimalControlSystem(
    const Matrix &q0, TFloat sigma, unsigned int N)
{
  // Copy parameters
  this->q0 = q0;
  this->sigma = sigma;
  this->N = N;
  this->k = q0.rows();
  this->dt = 1.0 / (N-1);

  // Allocate H derivatives
  for(unsigned int a = 0; a < VDim; a++)
    this->d_q__d_t[a].set_size(k);
}


template <class TFloat, unsigned int VDim>
TFloat
PointSetOptimalControlSystem<TFloat, VDim>
::ComputeEnergyAndVelocity(const Matrix &q, const Matrix &u)
{
  // Gaussian factor, i.e., K(z) = exp(f * z)
  TFloat f = -0.5 / (sigma * sigma);

  // Initialize the velocity vector to zero
  for(unsigned int a = 0; a < VDim; a++)
    {
    this->d_q__d_t[a].fill(0.0);
    }

  // Initialize the kinetic energy
  TFloat KE = 0.0;

  // Loop over all points
  for(unsigned int i = 0; i < k; i++)
    {
    // Get a pointer to pi for faster access?
    const TFloat *ui = u.data_array()[i], *qi = q.data_array()[i];

    // The diagonal terms
    for(unsigned int a = 0; a < VDim; a++)
      {
      KE += 0.5 * ui[a] * ui[a];
      d_q__d_t[a](i) += ui[a];
      }

    // Perform symmetric computation
    #pragma omp parallel for
    for(unsigned int j = i+1; j < k; j++)
      {
      const TFloat *uj = u.data_array()[j], *qj = q.data_array()[j];

      // Vector Qi-Qj
      VecD dq;

      // Dot product of Pi and Pj
      TFloat ui_uj = 0.0;

      // Compute above quantities
      for(unsigned int a = 0; a < VDim; a++)
        {
        dq[a] = qi[a] - qj[a];
        ui_uj += ui[a] * uj[a];
        }

      // Compute the Gaussian and its derivatives
      TFloat g = exp(f * dq.squared_magnitude()), g1 = f * g, g2 = f * g1;

      // Accumulate the Hamiltonian
      KE += ui_uj * g;

      // Accumulate the derivatives
      for(unsigned int a = 0; a < VDim; a++)
        {
        // First derivatives
        d_q__d_t[a](i) += g * uj[a];
        d_q__d_t[a](j) += g * ui[a];
        }
      } // loop over j
    } // loop over i

  return KE;
}


template <class TFloat, unsigned int VDim>
void
PointSetOptimalControlSystem<TFloat, VDim>
::FlowBackward(const MatrixArray &u, const MatrixArray &d_g__d_qt,
  TFloat w_kinetic,
  MatrixArray &d_f__d_u)
{
  // This function computes the gradient with respect to u[t] of the 
  // objective function
  //
  //    f(q[t],u[t]) = g(q[t]) + w_kinetic * KE
  //
  // given the gradient of g with respect to q[t]
  //
  // The function f can be written as (' denotes transpose)
  //
  //    f = g(q[t]) + w_kinetic * sum_{t=0}^{T-1} u[t]' (q[t+1]-q[t]) / dt
  //
  // The backflow uses the vector alpha which is updated using the recurrence
  //   alpha[T] = d_f__d_qt[T]
  //   alpha[t-1] = alpha[t] * Q(t,t-1) + d_f__d_qt[t-1]
  // where Q(t,t-1) is the Jacobian of q[t] with respect to q[t-1]
  //
  // Then the partial derivative d_f__d_u[t] is given by
  //   d_f__d_u[t-1] = alpha[t] * U(t-1) + w_kinetic * (q[t]-q[t-1]) / dt
  // where U(t-1) is the Jacobian of q[t] with respect to u[t-1]   
  
  // Allocate and initialize the alpha vector and the products of 
  // alpha with Q(t,t-1) and U(t-1)
  Vector alpha[VDim], alpha_Q[VDim], alpha_U[VDim];

  TFloat wke_factor = w_kinetic * 0.5;

  for(int a = 0; a < VDim; a++)
    {
    alpha[a] = d_g__d_qt[N-1].get_column(a) + wke_factor * u[N-2].get_column(a);
    alpha_Q[a].set_size(k);
    alpha_U[a].set_size(k);
    }

  // Work our way backwards
  for(int t = N-1; t > 0; t--)
    {
    // Propagate gradient backwards
    PropagateAlphaBackwards(Qt[t-1], u[t-1], alpha, alpha_Q, alpha_U);

    // Terms involved in KE computation
    Matrix delta_q = wke_factor * (Qt[t] - Qt[t-1]);
    Matrix delta_u = wke_factor * ((t > 1) ? (u[t-2] - u[t-1]) : -u[t-1]);

    for(int a = 0; a < VDim; a++)
      {
      // Update the gradient of f with respect to u[t-1]
      d_f__d_u[t-1].set_column(a, dt * alpha_U[a] + delta_q.get_column(a));

      // Update the alpha
      alpha[a] += dt * alpha_Q[a] + d_g__d_qt[t-1].get_column(a) + delta_u.get_column(a);
      }
    } 
}

template <class TFloat, unsigned int VDim>
void
PointSetOptimalControlSystem<TFloat, VDim>
::PropagateAlphaBackwards(
    const Matrix &q, const Matrix &u,
    const Vector alpha[], 
    Vector alpha_Q[], Vector alpha_U[])
{
  TFloat f = -0.5 / (sigma * sigma);

  for(int a = 0; a < VDim; a++)
    {
    alpha_Q[a].fill(0.0);

    // We must account for the diagonal term
    alpha_U[a] = alpha[a];
    }

  for(unsigned int i = 0; i < k; i++)
    {
    // Get a pointer to pi for faster access?
    const TFloat *ui = u.data_array()[i], *qi = q.data_array()[i];

    #pragma omp parallel for
    for(unsigned int j = i+1; j < k; j++)
      {
      const TFloat *uj = u.data_array()[j], *qj = q.data_array()[j];

      // Vector Qi-Qj
      VecD dq;

      // Compute above quantities
      for(unsigned int a = 0; a < VDim; a++)
        {
        dq[a] = qi[a] - qj[a];
        }

      // Compute the Gaussian and its derivatives
      TFloat g = exp(f * dq.squared_magnitude()), g1 = f * g;

      // Accumulate the derivatives
      for(unsigned int a = 0; a < VDim; a++)
        {
        TFloat term_2_g1_dqa = 2.0 * g1 * dq[a];
        TFloat alpha_j_ui_plus_alpha_i_uj = 0.0;

        for(unsigned int b = 0; b < VDim; b++)
          {
          alpha_j_ui_plus_alpha_i_uj += alpha[b][j] * ui[b] + alpha[b][i] * uj[b];
          }

        alpha_Q[a][i] += term_2_g1_dqa * alpha_j_ui_plus_alpha_i_uj;
        alpha_Q[a][j] -= term_2_g1_dqa * alpha_j_ui_plus_alpha_i_uj;

        alpha_U[a][i] += g * alpha[a][j];
        alpha_U[a][j] += g * alpha[a][i];
        }
      } // loop over j
    } // loop over i}
}

template <class TFloat, unsigned int VDim>
TFloat
PointSetOptimalControlSystem<TFloat, VDim>
::Flow(const std::vector<Matrix> &u)
{
  // Initialize q
  Matrix q = q0;

  // Allocate the streamline arrays
  Qt.resize(N); Qt[0] = q0;

  // The return value
  TFloat KE;

  // Flow over time
  for(unsigned int t = 1; t < N; t++)
    {
    // Compute the hamiltonian
    KE += dt * ComputeEnergyAndVelocity(q, u[t-1]);

    // Euler update
    for(unsigned int i = 0; i < k; i++)
      for(unsigned int a = 0; a < VDim; a++)
        q(i,a) += dt * d_q__d_t[a](i);

    // Store the flow results
    Qt[t] = q;
    }

  return KE;
}

template class PointSetOptimalControlSystem<double, 2>;
template class PointSetOptimalControlSystem<double, 3>;
template class PointSetOptimalControlSystem<float, 2>;
template class PointSetOptimalControlSystem<float, 3>;
