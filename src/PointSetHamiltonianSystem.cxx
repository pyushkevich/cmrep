#include "PointSetHamiltonianSystem.h"
#include <vnl/vnl_fastops.h>
#include <iostream>

template <class TFloat, unsigned int VDim>
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
    this->Hq[a].set_size(k);
    this->Hp[a].set_size(k);

    for(unsigned int b = 0; b < VDim; b++)
      {
      this->Hqq[a][b].set_size(k,k);
      this->Hqp[a][b].set_size(k,k);
      this->Hpp[a][b].set_size(k,k);
      }
    }
}


template <class TFloat, unsigned int VDim>
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
            TFloat val_qq = 2.0 * pi_pj * (2 * g2 * dq[a] * dq[b] + ((a == b) ? g1 : 0.0));
            Hqq[a][b](i,j) -= val_qq;
            Hqq[a][b](i,i) += val_qq;

            Hqp[a][b](i,j) += 2.0 * g1 * dq[a] * p(i,b);
            Hqp[a][b](i,i) += 2.0 * g1 * dq[a] * p(j,b);
            }

          Hpp[a][a](i,j) = g;
          }
        }
      } // loop over j
    } // loop over i

  return H;
}

template <class TFloat, unsigned int VDim>
TFloat
PointSetHamiltonianSystem<TFloat, VDim>
::FlowHamiltonian(const Matrix &p0, Matrix &q, Matrix &p)
{
  // Initialize q and p
  q = q0; p = p0;

  // Allocate the streamline arrays
  Qt.resize(N); Qt[0] = q0;
  Pt.resize(N); Pt[0] = p0;

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

    // Store the flow results
    Qt[t] = q; Pt[t] = p;

    // store the first hamiltonian value
    if(t == 1)
      H0 = H;
    }

  return H0;
}

extern "C" {
  int dgemm_(char *, char *, int *, int *, int *, double *, double *, int *, 
    double *, int *, double *, double *, int *);

  int sgemm_(char *, char *, int *, int *, int *, float *, float *, int *,
    float *, int *, float *, float *, int *);
};

/** WARNING - this is only meant for square matrices! */
template <class TFloat> class BlasInterface
{
public:
  typedef vnl_matrix<TFloat> Mat;
  static void add_AB_to_C(const Mat &A, const Mat &B, Mat &C);
  static void add_AtB_to_C(const Mat &A, const Mat &B, Mat &C);

private:
  static void gems(char *opA, char *opB, int *M, int *N, int *K, TFloat *alpha, TFloat *A, int *LDA,
                   TFloat *B, int *LDB, TFloat *beta, TFloat *C, int *LDC);

};

/** WARNING - this is only meant for square matrices! */
template <class TFloat>
void
BlasInterface<TFloat>
::add_AB_to_C(const Mat &A, const Mat &B, Mat &C)
{
  assert(
    A.rows() == B.rows() && A.rows() == C.rows() && A.rows() == A.columns() 
    && A.rows() == B.columns() && A.rows() == C.columns());

  char opA = 'N', opB = 'N';
  int M=A.rows(), N=M, K=M, LDA=K, LDB=N, LDC=M;
  TFloat alpha = 1.0, beta = 1.0;
  BlasInterface<TFloat>::gems(&opA, &opB, &M,&N,&K,&alpha,
    const_cast<TFloat *>(B.data_block()),&LDA,
    const_cast<TFloat *>(A.data_block()),&LDB,
    &beta,
    C.data_block(),&LDC);
}

template <class TFloat>
void
BlasInterface<TFloat>
::add_AtB_to_C(const Mat &A, const Mat &B, Mat &C)
{
  assert(
    A.rows() == B.rows() && A.rows() == C.rows() && A.rows() == A.columns() 
    && A.rows() == B.columns() && A.rows() == C.columns());

  char opA = 'N', opB = 'T';
  int M=A.rows(), N=M, K=M, LDA=K, LDB=N, LDC=M;
  TFloat alpha = 1.0, beta = 1.0;
  BlasInterface<TFloat>::gems(&opA, &opB, &M,&N,&K,&alpha,
    const_cast<TFloat *>(B.data_block()),&LDA,
    const_cast<TFloat *>(A.data_block()),&LDB,
    &beta,
    C.data_block(),&LDC);
}

template <>
void
BlasInterface<double>
::gems(char *opA, char *opB, int *M, int *N, int *K, double *alpha, double *A, int *LDA,
       double *B, int *LDB, double *beta, double *C, int *LDC)
{
  dgemm_(opA, opB, M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC);
}

template <>
void
BlasInterface<float>
::gems(char *opA, char *opB, int *M, int *N, int *K, float *alpha, float *A, int *LDA,
       float *B, int *LDB, float *beta, float *C, int *LDC)
{
  sgemm_(opA, opB, M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC);
}



template <class TFloat, unsigned int VDim>
void
PointSetHamiltonianSystem<TFloat, VDim>
::FlowGradientBackward(
  const Vector alpha1[VDim], 
  const Vector beta1[VDim],
  Vector result[VDim])
{
  // Allocate update vectors for alpha and beta
  Vector alpha[VDim], beta[VDim];
  Vector d_alpha[VDim], d_beta[VDim];
  for(int a = 0; a < VDim; a++)
    {
    alpha[a] = alpha1[a];
    beta[a] = beta1[a];
    d_alpha[a].set_size(k);
    d_beta[a].set_size(k);
    }

  // Work our way backwards
  for(int t = N-1; t > 0; t--)
    {
    // Compute the Hamiltonian on the (i-1)st q/p
    ComputeHamiltonianJet(Qt[t - 1], Pt[t - 1], true);

    // Compute the updates
    for(int a = 0; a < VDim; a++)
      {
      d_alpha[a].fill(0.0);
      d_beta[a].fill(0.0);

      for(int b = 0; b < VDim; b++)
        {
        d_alpha[a] += alpha[b] * (Hqp[a][b]).transpose() - beta[b] * Hqq[b][a];
        d_beta[a] += alpha[b] * Hpp[b][a] - beta[b] * Hqp[b][a];
        }
      }

    // Update the vectors
    for(int a = 0; a < VDim; a++)
      {
      alpha[a] += dt * d_alpha[a];
      beta[a] += dt * d_beta[a];
      }
    } 

  // Finally, what we are really after are the betas
  for(int a = 0; a < VDim; a++)
    {
    result[a] = beta[a];
    }
}

template <class TFloat, unsigned int VDim>
void
PointSetHamiltonianSystem<TFloat, VDim>
::InterpolateVelocity(unsigned int t, const TFloat *x, TFloat *v)
{
  // Gaussian factor, i.e., K(z) = exp(f * z)
  TFloat f = -0.5 / (sigma * sigma);

  // Initialize v to zero
  for(unsigned int a = 0; a < VDim; a++)
    v[a] = 0.0;

  // Compute the velocity for this point
  for(unsigned int i = 0; i < k; i++)
    {
    double dsq = 0.0;
    for(unsigned int a = 0; a < VDim; a++)
      {
      double da = Qt[t](i,a) - x[a];
      dsq += da * da;
      }
    double Kq = exp(dsq * f);
    for(unsigned int a = 0; a < VDim; a++)
      v[a] += Kq * Pt[t](i,a);
    }
}

template <class TFloat, unsigned int VDim>
TFloat
PointSetHamiltonianSystem<TFloat, VDim>
::FlowHamiltonianWithGradient(
  const Matrix &p0, Matrix &q, Matrix &p,
  Matrix grad_q[VDim][VDim], Matrix grad_p[VDim][VDim])
{
  // Initialize q and p
  q = q0; p = p0;

  // Allocate the streamline arrays
  Qt.resize(N); Qt[0] = q0;
  Pt.resize(N); Pt[0] = p0;

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

    Qt[t] = q; Pt[t] = p;

    // The nastiest part - some matrix multiplications
    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int b = 0; b < VDim; b++)
        {
        gupd_q[a][b].fill(0.0);
        gupd_p[a][b].fill(0.0);
        
        for(unsigned int c = 0; c < VDim; c++)
          {
          BlasInterface<TFloat>::add_AB_to_C(Hqp[a][c], grad_p[c][b], gupd_p[a][b]);
          BlasInterface<TFloat>::add_AB_to_C(Hqq[a][c], grad_q[c][b], gupd_p[a][b]);
          BlasInterface<TFloat>::add_AB_to_C(Hpp[a][c], grad_p[c][b], gupd_q[a][b]); 
          BlasInterface<TFloat>::add_AtB_to_C(Hqp[c][a], grad_q[c][b], gupd_q[a][b]); 

          // gupd_p[a][b] += Hqp[a][c] * grad_p[c][b] + Hqq[a][c] * grad_q[c][b];
          // gupd_q[a][b] += Hpp[a][c] * grad_p[c][b] + Hqp[c][a].transpose() * grad_q[c][b];
          }
        }
      }

    for(unsigned int a = 0; a < VDim; a++)
      {
      for(unsigned int b = 0; b < VDim; b++)
        {
        grad_q[a][b] += dt * gupd_q[a][b];
        grad_p[a][b] -= dt * gupd_p[a][b];
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
template class PointSetHamiltonianSystem<float, 2>;
template class PointSetHamiltonianSystem<float, 3>;
