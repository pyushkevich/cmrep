#include "PointSetHamiltonianSystem.h"

template <class TFloat, class VDim>
PointSetHamiltonianSystem<TFloat, VDim>
::PointSetHamiltonianSystem(
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
    {
    this->Hq.set_size(k, 0.0);
    this->Hp.set_size(k, 0.0);

    for(unsigned int b = 0; b < VDim; b++)
      {
      this->Hqq[a][b].set_size(k,k,0.0);
      this->Hqp[a][b].set_size(k,k,0.0);
      this->Hpp[a][b].set_size(k,k,0.0);
      }
    }
}


template <class TFloat, class VDim>
TFloat
PointSetHamiltonianSystem<TFloat, VDim>
::ComputeHamiltonianJet(const Matrix &q, const Matrix &p, bool flag_hessian)
{
  // Gaussian factor, i.e., K(z) = exp(f * z)
  TFloat f = -0.5 / (sigma * sigma);

  // Initialize the gradient and Hessian terms to zeros
  for(unsigned int a = 0; a < VDim; a++)
    {
    this->Hq[a].fill(0.0);
    this->Hp[a].fill(0.0);

    if(flag_hessian)
      {
      for(unsigned int b = 0; b < VDim; b++)
        {
        this->Hqq[a][b].fill(0.0);
        this->Hqp[a][b].fill(0.0);
        this->Hpp[a][b].fill(0.0);
        }
      }
    }

  // Initialize hamiltonian
  TFloat H = 0.0;

  // Loop over all points
  for(unsigned int i = 0; i < k; i++)
    {
    // TODO: you should be able to do this computation on half the matrix, it's symmetric!
    #pragma omp parallel for
    for(unsigned int j = 0; j < k; j++)
      {
      // Vector Qi-Qj
      VecD dq;

      // Dot product of Pi and Pj
      TFloat pi_pj = 0.0;

      // Compute above quantities
      for(unsigned int a = 0; a < VDim; a++)
        {
        dq[a] = q(i,a) - q(j,a);
        pi_pj += p(i,a) * p(j,a);
        }

      // Compute the Gaussian and its derivatives
      TFloat g = exp(f * dq.squared_magnitude()), g1 = f * g, g2 = f * g1;

      // Accumulate the Hamiltonian
      H += 0.5 * pi_pj * g;

      // Accumulate the derivatives
      for(unsigned int a = 0; a < VDim; a++)
        {
        // First derivatives
        Hq[a](i) += 2 * pi_pj * g1 * dq[a];
        Hp[a](i) += g * p(j,a);

        // Second derivatives
        if(flag_hessian)
          {
          for(unsigned int b = 0; b < VDim; b++)
            {
            TFloat val_qq = 2.0 * pi_pj * (2 * g2 * dq[a] * dq[b] + g1);
            Hqq[a][b](i,j) -= val;
            Hqq[a][b](i,i) += val;

            Hqp[a][b](i,j) += 2.0 * g1 * dq[a] * p(i,a);
            Hqp[a][b](i,i) += 2.0 * g1 * dq[a] * p(j,a);
            }

          Hpp[a][a](i,j) = g;
          }
        }
      } // loop over j
    } // loop over i

  return H;
}

template <class TFloat, class VDim>
TFloat
PointSetHamiltonianSystem<TFloat, VDim>
::FlowHamiltonian(const Matrix &p0, Matrix &q, Matrix &p)
{
  // Initialize q and p
  q = q0; p = p0;

  // The return value
  TFloat H, H0;

  // Flow over time
  for(unsigned int t = 1; t < N; t++)
    {
    // Compute the hamiltonian
    H = ComputeHamiltonianJet(q, p, false);

    // Euler update
    #pragma omp parallel for
    for(unsigned int i = 0; i < k; i++)
      {
      for(unsigned int a = 0; a < VDim; a++)
        {
        q(i,a) += dt * Hp[a](i);
        p(i,a) -= dt * Hq[a](i);
        }
      }

    // store the first hamiltonian value
    if(t == 1)
      H0 = H;
    }

  return H0;
}

template <class TFloat, class VDim>
TFloat
PointSetHamiltonianSystem<TFloat, VDim>
::FlowHamiltonianWithGradient(
  const Matrix &p0, Matrix &q, Matrix &p,
  Matrix grad_q[VDim][VDim], Matrix grad_p[VDim][VDim])
{
  // Initialize q and p
  q = q0; p = p0;

  // We need temporary matrices to store the updates of gradient
  Matrix gupd_q[VDim][VDim], gupd_p[VDim][VDim];

  // Initialize the gradients
  for(unsigned int a = 0; a < VDim; a++)
    {
    for(unsigned int b = 0; b < VDim; b++)
      {
      grad_q[a][b].fill(0.0);
      if(a == b)
        grad_p[a][b].set_identity();
      else
        grad_p[a][b].fill(0.0);

      gupd_p[a][b].set_size(k,k);
      gupd_q[a][b].set_size(k,k);
      }
    }

  // The return value
  TFloat H, H0;

  // Flow over time
  for(unsigned int t = 1; t < N; t++)
    {
    // Compute the hamiltonian
    H = ComputeHamiltonianJet(q, p, true);

    // Euler update
    #pragma omp parallel for
    for(unsigned int i = 0; i < k; i++)
      {
      for(unsigned int a = 0; a < VDim; a++)
        {
        q(i,a) += dt * Hp[a](i);
        p(i,a) -= dt * Hq[a](i);
        }
      }

    // The nastiest part - some matrix multiplications
    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int b = 0; b < VDim; b++)
        {
        gupd_q[a][b].fill(0.0);
        gupd_p[a][b].fill(0.0);
        
        for(unsigned int c = 0; c < VDim; c++)
          {
          gupd_p[a][b] += Hqp[a][c] * grad_p[c][b] + Hqq[a][c] * grad_q[c][b];
          gupd_q[a][b] += Hpp[a][c] * grad_p[c][b] + Hqp[c][a].transpose() * grad_q[c][b];
          }
        }
      }

    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int b = 0; b < VDim; b++)
        {
        grad_q[a][b] += gupd_q[a][b];
        grad_p[a][b] -= gupd_p[a][b];
        }
      }

    // store the first hamiltonian value
    if(t == 1)
      H0 = H;
    }

  return H0;
}

template class PointSetHamiltonianSystem<double, 2>;
template class PointSetHamiltonianSystem<double, 3>;
