// Include the IPOPT header
#include "IpTNLP.hpp"

using namespace Ipopt;

class MedialContrainedOptimizationProblem : public TNLP
{
public:
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

};


bool
MedialContrainedOptimizationProblem::
get_bounds_info(Index n, Number* x_l, Number* x_u,
                Index m, Number* g_l, Number* g_u);
{
  // Set the lower and upper bounds for the variables. The bounds for X are
  // just the imaging cube (or none, in reality), and the bounds on R is that
  // it must be positive
  for(size_t i = 0; i < n; i++)
    {
    if(xModel->IsCoefficientBoundedBelow(i))
      x_l[i] = xModel->GetCoefficientLowerBound(i);
    else
      x_l[i] = nlp_lower_bound_inf;

    if(xModel->IsCoefficientBoundedAbove(i))
      x_u[i] = xModel->GetCoefficientUpperBound(i);
    else
      x_u[i] = nlp_upper_bound_inf;
    }

  // For now, let's assume all the constraints are equality constraints. That
  // may change one day, but not yet
  for(size_t j = 0; j < m; j++)
    {
    g_l[j] = g_u[j] = 0.0;
    }
}



bool
MedialContrainedOptimizationProblem::
get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  // n is the number of unknows
  n = xModel->GetNumberOfCoefficients();

  // m is the number of constraints
  m = xModel->GetNumberOfConstraints();

  // nnz_jac_g is the number of non-zero entries in the jacobian
  nnz_jac_g = 0;
  for(size_t i = 0; i < m; i++)
    {
      nnz_jac_g += xModel->GetConstraintGradientSize(i);
    }
  nnz_h_lag = 0;
  index_style = 0;
}

virtual bool
MedialContrainedOptimizationProblem
::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{

}