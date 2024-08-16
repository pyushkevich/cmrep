#include "IPOptQCQPProblemInterface.h"

#include "itkTimeProbe.h"

using namespace Ipopt;

/******************************************************************
  IPOPT STUFF
  *****************************************************************/

IPOptQCQPProblemInterface::IPOptQCQPProblemInterface(qcqp::Problem &p)
  : m_Problem(p)
{
  m_ConstraintLogFile = NULL;
}

void IPOptQCQPProblemInterface::log_constraints(FILE *flog)
{
  m_ConstraintLogFile = flog;
}

bool IPOptQCQPProblemInterface::get_nlp_info(
    Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag,
    IndexStyleEnum &index_style)
{
  n = m_Problem.GetNumberOfVariables();
  m = m_Problem.GetNumberOfConstraints();
  nnz_jac_g = m_Problem.GetConstraintsJacobian().GetNumberOfSparseValues();
  nnz_h_lag = m_Problem.GetHessianOfLagrangean().GetNumberOfSparseValues();
  index_style = TNLP::C_STYLE;
  return true;
}

bool IPOptQCQPProblemInterface::get_bounds_info(
    Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u)
{
  for(int i = 0; i < n; i++)
    m_Problem.GetVariableBounds(i, x_l[i], x_u[i]);

  for(int j = 0; j < m; j++)
    m_Problem.GetConstraintBounds(j, g_l[j], g_u[j]);

  return true;
}

bool IPOptQCQPProblemInterface::get_starting_point(
    Index n, bool init_x, Number *x, bool init_z,
    Number *z_L, Number *z_U, Index m, bool init_lambda, Number *lambda)
{
  // Set the primal variables
  for(int i = 0; i < n; i++)
    x[i] = m_Problem.GetVariableValue(i);

  // Set the warm start data
  if(init_z)
    {
    if(m_Problem.GetWarmStartZL().size() && m_Problem.GetWarmStartZU().size())
      {
      for(int i = 0; i < n; i++)
        {
        z_L[i] = m_Problem.GetWarmStartZL()[i];
        z_U[i] = m_Problem.GetWarmStartZU()[i];
        }
      }
    else
      return false;
    }

  if(init_lambda)
    {
    if(m_Problem.GetWarmStartLambda().size())
      {
      for(int i = 0; i < m; i++)
        lambda[i] = m_Problem.GetWarmStartLambda()[i];
      }
    else
      return false;
    }

  return true;
}

bool IPOptQCQPProblemInterface::eval_f(
    Index n, const Number *x, bool new_x, Number &obj_value)
{
  // Compute the objective
  obj_value = m_Problem.EvaluateLoss(x);
  return true;
}

bool IPOptQCQPProblemInterface::eval_grad_f(
    Index n, const Number *x, bool new_x, Number *grad_f)
{
  // Compute the gradient
  m_Problem.EvaluateLossGradient(x, grad_f);
  return true;
}

bool IPOptQCQPProblemInterface::eval_g(
    Index n, const Number *x, bool new_x, Index m, Number *g)
{
  if(m_ConstraintLogFile)
    {
    // Keep track of the constraints
    std::map<std::string, double> catmap;
    std::string lastcat;

    for(unsigned int i = 0; i < m; i++)
      {
      // Evaluate the constraint
      g[i] = m_Problem.EvaluateConstraint(x, i);

      // Get the lower and upper bounds
      double lb, ub;
      m_Problem.GetConstraintBounds(i, lb, ub);

      // Compute violation
      double viol = g[i] > ub ? (g[i] - ub) : (g[i] < lb ? lb - g[i] : 0.0);
      double &max_viol = catmap[m_Problem.GetConstraintName(i)];
      max_viol = std::max(max_viol, viol);
      }

    // Print log
    if(m_ConstraintLogIter % 10 == 0)
      {
      for(auto it : catmap)
        fprintf(m_ConstraintLogFile, "%12s ", it.first.c_str());
      fprintf(m_ConstraintLogFile, "\n");
      }
    for(auto it : catmap)
      fprintf(m_ConstraintLogFile, "%12.8f ", it.second);
    fprintf(m_ConstraintLogFile, "\n");
    fflush(m_ConstraintLogFile);
    m_ConstraintLogIter++;
    }
  else
    {
    // If not logging constraints, just one line
    for(unsigned int i = 0; i < m; i++)
      g[i] = m_Problem.EvaluateConstraint(x, i);

    }

  return true;
}

bool IPOptQCQPProblemInterface::eval_jac_g(
    Index n, const Number *x, bool new_x, Index m,
    Index nele_jac, Index *iRow, Index *jCol, Number *values)
{
  // Get the Jacobian sparse matrix
  using SparseTensor = qcqp::Problem::SparseTensor;
  SparseTensor &T = m_Problem.GetConstraintsJacobian();

  if(iRow && jCol)
    {
    // A request for the sparsity structure. We need to iterate the non-NULL
    // entries of the sparse array
    int iEntry = 0;
    for(int row = 0; row < T.GetNumberOfRows(); row++)
      {
      for(auto it = T.Row(row); !it.IsAtEnd(); ++it)
        {
        iRow[iEntry] = row;
        jCol[iEntry] = it.Column();
        iEntry++;
        }
      }
    }
  else
    {
    // A request for the Jacobian values
    int iEntry = 0;
    for(size_t row = 0; row < T.GetNumberOfRows(); row++)
      {
      for(auto it = T.Row(row); !it.IsAtEnd(); ++it)
        {
        auto &wz = it.Value();
        double val = wz.z;
        for(auto it_w : wz.w)
          val += x[std::get<0>(it_w)] * std::get<1>(it_w);
        values[iEntry++] = val;
        }
      }
    }

  return true;
}


bool IPOptQCQPProblemInterface::eval_h(
    Index n, const Number *x, bool new_x, Number obj_factor, Index m,
    const Number *lambda, bool new_lambda, Index nele_hess,
    Index *iRow, Index *jCol, Number *values)
{
  // Get the Hessian sparse matrix
  using SparseTensor = qcqp::Problem::SparseTensor;
  SparseTensor &H = m_Problem.GetHessianOfLagrangean();

  if(iRow && jCol)
    {
    // A request for the sparsity structure. We need to iterate the non-NULL
    // entries of the sparse array
    int iEntry = 0;
    for(int row = 0; row < H.GetNumberOfRows(); row++)
      {
      for(auto it = H.Row(row); !it.IsAtEnd(); ++it)
        {
        iRow[iEntry] = row;
        jCol[iEntry] = it.Column();
        iEntry++;
        }
      }
    }
  else
    {
    // A request for the Jacobian values
    int iEntry = 0;
    for(int row = 0; row < H.GetNumberOfRows(); row++)
      {
      for(auto it = H.Row(row); !it.IsAtEnd(); ++it)
        {
        auto &wz = it.Value();
        double val = wz.z * obj_factor;
        for(auto it_w : wz.w)
          val += lambda[std::get<0>(it_w)] * std::get<1>(it_w);

        values[iEntry++] = val;
        }
      }
    }

  return true;
}

void IPOptQCQPProblemInterface::finalize_solution(
    SolverReturn status, Index n, const Number *x,
    const Number *z_L, const Number *z_U, Index m,
    const Number *g, const Number *lambda, Number obj_value,
    const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq)
{
  for(int i = 0; i < n; i++)
    m_Problem.SetVariableValue(i, x[i]);

  m_Problem.SetWarmStartLambda(m, lambda);
  m_Problem.SetWarmStartZL(n, z_L);
  m_Problem.SetWarmStartZU(n, z_U);

  m_Problem.Finalize();
}
