#include "PointSetHamiltonianSystem.h"
#include <iostream>
#include <cstdlib>
#include <vnl/vnl_matlab_filewrite.h>
#include <vnl/vnl_matlab_read.h>
#include <vnl/vnl_random.h>
#include <cstdarg>

void test(double val1, double val2, double tol, const char *format, ...)
{
  if(fabs(val1 - val2) > tol)
    {
    char buffer[256];
    va_list args;
    va_start (args, format);
    vsprintf (buffer,format, args);
    va_end (args);

    std::cerr << "Mismatch in " << buffer << " ";
    std::cerr << val1 << " vs " << val2 << std::endl;
    exit(-1);
    }
}

int main(int argc, char *argv[])
{
  typedef PointSetHamiltonianSystem<double, 2> Ham;

  // Read regression data
  std::ifstream iss(argv[1]);

  // Input data
  vnl_vector<double> in_q0[2], in_p0[2];
  vnl_matlab_read_or_die(iss, in_q0[0], "qx");
  vnl_matlab_read_or_die(iss, in_q0[1], "qy");
  vnl_matlab_read_or_die(iss, in_p0[0], "px");
  vnl_matlab_read_or_die(iss, in_p0[1], "py");
  
  // Initialize the input vectors
  int k = 80;
  Ham::Matrix q0(k, 2), p0(k, 2), q1(k, 2), p1(k, 2);
  for(int a = 0; a < 2; a++)
    {
    q0.set_column(a, in_q0[a]);
    p0.set_column(a, in_p0[a]);
    }

  // Initialize the gradients
  Ham::Matrix grad_q[2][2], grad_p[2][2];
  for(int a = 0; a < 2; a++)
    {
    for(int b = 0; b < 2; b++)
      {
      grad_p[a][b].set_size(k,k);
      grad_q[a][b].set_size(k,k);
      }
    }

  // Create a hamiltonian system
  Ham hsys(q0, 0.08, 100);

  // Compute the Hamiltonian
  double H  = hsys.ComputeHamiltonianJet(q0, p0, true);

  // Read the regression data for the Hamiltonian and derivatives
  vnl_vector<double> r_H, r_Hp[2], r_Hq[2];
  vnl_matrix<double> r_Hqq[2][2], r_Hqp[2][2], r_Hpp[2][2];
  vnl_matlab_read_or_die(iss, r_H, "H");
  vnl_matlab_read_or_die(iss, r_Hp[0], "Hpx");
  vnl_matlab_read_or_die(iss, r_Hp[1], "Hpy");
  vnl_matlab_read_or_die(iss, r_Hq[0], "Hqx");
  vnl_matlab_read_or_die(iss, r_Hq[1], "Hqy");
  vnl_matlab_read_or_die(iss, r_Hqq[0][0], "Hqxqx");
  vnl_matlab_read_or_die(iss, r_Hqq[0][1], "Hqxqy");
  vnl_matlab_read_or_die(iss, r_Hqq[1][1], "Hqyqy");
  vnl_matlab_read_or_die(iss, r_Hqp[0][0], "Hqxpx");
  vnl_matlab_read_or_die(iss, r_Hqp[0][1], "Hqxpy");
  vnl_matlab_read_or_die(iss, r_Hqp[1][0], "Hqypx");
  vnl_matlab_read_or_die(iss, r_Hqp[1][1], "Hqypy");
  vnl_matlab_read_or_die(iss, r_Hpp[0][0], "Hpxpx");
  vnl_matlab_read_or_die(iss, r_Hpp[0][1], "Hpxpy");
  vnl_matlab_read_or_die(iss, r_Hpp[1][1], "Hpypy");

  // Symmetric terms
  r_Hqq[1][0] = r_Hqq[0][1];
  r_Hpp[1][0] = r_Hpp[0][1];

  // Regression testing on the Hamiltonian
  test(H, r_H(0), 1e-8, "H");
  for(int a = 0; a < 2; a++)
    {
    for(int j = 0; j < k; j++)
      {
      test(hsys.GetHp(a)(j), r_Hp[a](j), 1e-8, "Hp[%d][%d]", a, j); 
      test(hsys.GetHq(a)(j), r_Hq[a](j), 1e-8, "Hq[%d][%d]", a, j); 

      for(int b = 0; b < 2; b++)
        {
        for(int l = 0; l < k; l++)
          {
          test(hsys.GetHpp(a,b)(j,l), r_Hpp[a][b](j,l), 1e-8, "Hpp[%d][%d](%d,%d)",a,b,j,l);
          test(hsys.GetHqp(a,b)(j,l), r_Hqp[a][b](j,l), 1e-8, "Hqp[%d][%d](%d,%d)",a,b,j,l);
          test(hsys.GetHqq(a,b)(j,l), r_Hqq[a][b](j,l), 1e-8, "Hqq[%d][%d](%d,%d)",a,b,j,l);
          }
        }
      }
    }

  // Passed this stage of the test
  std::cout << "Passed regression test on the Hamiltonian derivatives" << std::endl;

  double t_start, t_finish;

  // Flow the system without gradient - to see how long it takes
  t_start = clock();
  hsys.FlowHamiltonian(p0, q1, p1);
  t_finish = clock();
  std::cout << "Flow without gradient computed in " 
    << (t_finish - t_start) / CLOCKS_PER_SEC << " sec" << std::endl;

  // Flow the system without gradient
  t_start = clock();
  hsys.FlowHamiltonianWithGradient(p0, q1, p1, grad_q, grad_p);
  t_finish = clock();
  std::cout << "Flow with gradient computed in " 
    << (t_finish - t_start) / CLOCKS_PER_SEC << " sec" << std::endl;

  // Read the regression data for the flow and derivatives
  vnl_vector<double> r_q1[2], r_p1[2];
  vnl_matrix<double> r_gradQ[2][2], r_gradP[2][2];
  vnl_matlab_read_or_die(iss, r_q1[0], "qx_t");
  vnl_matlab_read_or_die(iss, r_q1[1], "qy_t");
  vnl_matlab_read_or_die(iss, r_p1[0], "px_t");
  vnl_matlab_read_or_die(iss, r_p1[1], "py_t");
  vnl_matlab_read_or_die(iss, r_gradQ[0][0], "qx_nx");
  vnl_matlab_read_or_die(iss, r_gradQ[0][1], "qx_ny");
  vnl_matlab_read_or_die(iss, r_gradQ[1][0], "qy_nx");
  vnl_matlab_read_or_die(iss, r_gradQ[1][1], "qy_ny");
  vnl_matlab_read_or_die(iss, r_gradP[0][0], "px_nx");
  vnl_matlab_read_or_die(iss, r_gradP[0][1], "px_ny");
  vnl_matlab_read_or_die(iss, r_gradP[1][0], "py_nx");
  vnl_matlab_read_or_die(iss, r_gradP[1][1], "py_ny");

  for(int i = 0; i < k; i++)
    {
    for(int a = 0; a < 2; a++)
      {
      test(q1(i,a), r_q1[a](i), 1e-8, "q1[%d][%d]", a, i);
      test(p1(i,a), r_p1[a](i), 1e-8, "p1[%d][%d]", a, i);

      for(int j = 0; j < k; j++)
        {
        for(int b = 0; b < 2; b++)
          {
          test(grad_q[a][b](i,j), r_gradQ[a][b](i,j), 1e-8, "grad_q[%d][%d](%d,%d)",a,b,i,j);
          test(grad_p[a][b](i,j), r_gradP[a][b](i,j), 1e-8, "grad_p[%d][%d](%d,%d)",a,b,i,j);
          }
        }
      }
    }

  // Passed this stage of the test
  std::cout << "Passed regression test on the Hamiltonian flow" << std::endl;

  // Test forward-backward gradient computation approach
  Ham::Vector alpha[2], beta[2], dd[2], dd_test[2];
  vnl_random rnd(1234);
  for(int a = 0; a < 2; a++)
    {
    alpha[a].set_size(k);
    beta[a].set_size(k);

    for(int i = 0; i < k; i++)
      {
      alpha[a](k) = rnd.drand32(-1.0, 1.0);
      beta[a](k) = rnd.drand32(-1.0, 1.0);
      }
    }

  hsys.FlowHamiltonian(p0, q1, p1);

  t_start = clock();
  hsys.FlowGradientBackward(alpha, beta, dd);
  t_finish = clock();
  std::cout << "Flow without gradient computed in " 
    << (t_finish - t_start) / CLOCKS_PER_SEC << " sec" << std::endl;

  // Compare the derivative we got with actual computation of the derivative
  for(int a = 0; a < 2; a++)
    {
    dd_test[a].set_size(k);
    dd_test[a].fill(0.0);
    for(int b = 0; b < 2; b++)
      {
      dd_test[a] += alpha[b] * grad_q[b][a] + beta[b] * grad_p[b][a];
      }

    for(int i = 0; i < k; i++)
      {
      test(dd[a](i), dd_test[a](i), 1.e-6, "Backflown Gradient [%d](%d)", a, i);
      }
    }

  std::cout << "Passed check on backward gradient flow" << std::endl;

  // Test the gradient vs. numerical gradient
  double eps = 1e-5;
  Ham::Matrix q1_1 = q1, p1_1 = p1, q1_2 = q1, p1_2 = p1;
  for(int i = 0; i < k; i++)
    {
    for(int a = 0; a < 2; a++)
      {
      Ham::Matrix p0_1 = p0; p0_1(i,a) -= eps;
      Ham::Matrix p0_2 = p0; p0_2(i,a) += eps;

      hsys.FlowHamiltonian(p0_1, q1_1, p1_1);
      hsys.FlowHamiltonian(p0_2, q1_2, p1_2);

      // Numeric derivative of p1 and q1 wrt p0(i,a)
      Ham::Matrix Dp_num = (p1_2 - p1_1) * (0.5/eps);
      Ham::Matrix Dq_num = (q1_2 - q1_1) * (0.5/eps);

      // Compare with the analytic derivative
      for(int j = 0; j < k; j++)
        {
        for(int b = 0; b < 2; b++)
          {
          test(Dp_num(j,b), grad_p[b][a](j,i), 1e-5, "grad_p[%d][%d][%d][%d]",b,a,j,i);
          test(Dq_num(j,b), grad_q[b][a](j,i), 1e-5, "grad_q[%d][%d][%d][%d]",b,a,j,i);
          }
        }
      }
    }

  // Report that we are done
  std::cout << "Passed derivative check on Hamilton flow" << std::endl;

  // Create a larger problem with 800 time points, to test overhead
  int big_k = 800;
  Ham::Matrix big_q0(big_k, 2), big_p0(big_k, 2), big_q1(big_k, 2), big_p1(big_k, 2);
  for(int i = 0; i < big_k; i++)
    {
    double t = 1.0 / (big_k - 1.0);
    double p = (k - 1.0) * t;
    int i0 = floor(p), i1 = ceil(p);
    double v = p - i0;
    for(int a = 0; a < 2; a++)
      {
      big_q0(i,a) = q0(i0, a) * (1.0 - v) + q0(i1, a) * v;
      big_p0(i,a) = p0(i0, a) * (1.0 - v) + p0(i1, a) * v;
      }
    }

  Ham big_hsys(big_q0, 0.08, 100);
  t_start = clock();
  big_hsys.FlowHamiltonian(big_p0, big_q1, big_p1);
  t_finish = clock();
  std::cout << "Big problem: Flow without gradient computed in " 
    << (t_finish - t_start) / CLOCKS_PER_SEC << " sec" << std::endl;
  

  return 0;
};
