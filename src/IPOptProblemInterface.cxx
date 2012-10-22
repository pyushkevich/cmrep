#include "IPOptProblemInterface.h"

#include "itkTimeProbe.h"

using namespace gnlp;
using namespace Ipopt;

/******************************************************************
  IPOPT STUFF
  *****************************************************************/

IPOptProblemInterface::IPOptProblemInterface(
    gnlp::ConstrainedNonLinearProblem *p)
{
  m_Problem = p;
}

bool IPOptProblemInterface::get_nlp_info(
    Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag,
    IndexStyleEnum &index_style)
{
  n = m_Problem->GetNumberOfVariables();
  m = m_Problem->GetNumberOfConstraints();
  nnz_jac_g = m_Problem->GetConstraintsJacobian().GetNumberOfSparseValues();
  nnz_h_lag = m_Problem->GetHessianOfLagrangean().GetNumberOfSparseValues();
  index_style = TNLP::C_STYLE;
  return true;
}

bool IPOptProblemInterface::get_bounds_info(
    Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u)
{
  for(int i = 0; i < n; i++)
    m_Problem->GetVariableBounds(i, x_l[i], x_u[i]);

  for(int j = 0; j < m; j++)
    m_Problem->GetConstraintBounds(j, g_l[j], g_u[j]);

  return true;
}

bool IPOptProblemInterface::get_starting_point(
    Index n, bool init_x, Number *x, bool init_z,
    Number *z_L, Number *z_U, Index m, bool init_lambda, Number *lambda)
{
  assert(init_z == false && init_lambda == false);
  for(int i = 0; i < n; i++)
    x[i] = m_Problem->GetVariableValue(i);

  return true;
}

bool IPOptProblemInterface::eval_f(
    Index n, const Number *x, bool new_x, Number &obj_value)
{
  // Set the values of all the variables
  if(new_x)
    m_Problem->SetVariableValues(x);

  // Compute the objective
  obj_value = m_Problem->GetObjective()->Evaluate();

  return true;
}

bool IPOptProblemInterface::eval_grad_f(
    Index n, const Number *x, bool new_x, Number *grad_f)
{
  // Set the values of all the variables
  if(new_x)
    m_Problem->SetVariableValues(x);

  // Get the partial derivatives of the objective function
  for(unsigned int i = 0; i < n; i++)
    {
    Expression *pd = m_Problem->GetObjectivePD(i);
    grad_f[i] = pd ? pd->Evaluate() : 0.0;
    }

  return true;
}

bool IPOptProblemInterface::eval_g(
    Index n, const Number *x, bool new_x, Index m, Number *g)
{
  // Set the values of all the variables
  if(new_x)
    m_Problem->SetVariableValues(x);

  // Get the partial derivatives of the objective function
  for(unsigned int i = 0; i < m; i++)
    g[i] = m_Problem->GetConstraint(i)->Evaluate();


  return true;
}

bool IPOptProblemInterface::eval_jac_g(
    Index n, const Number *x, bool new_x, Index m,
    Index nele_jac, Index *iRow, Index *jCol, Number *values)
{
  // Get the Jacobian sparse matrix
  typedef ConstrainedNonLinearProblem::SparseExpressionMatrix SparseMat;
  SparseMat DG = m_Problem->GetConstraintsJacobian();

  if(iRow && jCol)
    {
    // A request for the sparsity structure. We need to iterate the non-NULL
    // entries of the sparse array
    int iEntry = 0;
    for(int row = 0; row < DG.GetNumberOfRows(); row++)
      {
      for(SparseMat::RowIterator it = DG.Row(row); !it.IsAtEnd(); ++it)
        {
        iRow[iEntry] = row;
        jCol[iEntry] = it.Column();
        iEntry++;
        }
      }
    }
  else
    {
    // Set the values of all the variables
    if(new_x)
      m_Problem->SetVariableValues(x);

    // A request for the Jacobian values
    int iEntry = 0;
    for(int row = 0; row < DG.GetNumberOfRows(); row++)
      {
      for(SparseMat::RowIterator it = DG.Row(row); !it.IsAtEnd(); ++it)
        {
        values[iEntry++] = it.Value()->Evaluate();
        }
      }
    }

  return true;
}


bool IPOptProblemInterface::eval_h(
    Index n, const Number *x, bool new_x, Number obj_factor, Index m,
    const Number *lambda, bool new_lambda, Index nele_hess,
    Index *iRow, Index *jCol, Number *values)
{
  // Get the Jacobian sparse matrix
  typedef ConstrainedNonLinearProblem::SparseExpressionMatrix SparseMat;
  SparseMat H = m_Problem->GetHessianOfLagrangean();

  if(iRow && jCol)
    {
    // A request for the sparsity structure. We need to iterate the non-NULL
    // entries of the sparse array
    int iEntry = 0;
    for(int row = 0; row < H.GetNumberOfRows(); row++)
      {
      for(SparseMat::RowIterator it = H.Row(row); !it.IsAtEnd(); ++it)
        {
        iRow[iEntry] = row;
        jCol[iEntry] = it.Column();
        iEntry++;
        }
      }
    }
  else
    {
    // Set the values of all the variables
    if(new_x)
      m_Problem->SetVariableValues(x);

    // Set the values of the lamdas
    if(new_lambda)
      m_Problem->SetLambdaValues(lambda);

    // Set the objective factor
    m_Problem->SetSigma(obj_factor);

    // A request for the Jacobian values
    int iEntry = 0;
    for(int row = 0; row < H.GetNumberOfRows(); row++)
      {
      for(SparseMat::RowIterator it = H.Row(row); !it.IsAtEnd(); ++it)
        {
        values[iEntry++] = it.Value()->Evaluate();
        }
      }
    }

  return true;
}

void IPOptProblemInterface::finalize_solution(
    SolverReturn status, Index n, const Number *x,
    const Number *z_L, const Number *z_U, Index m,
    const Number *g, const Number *lambda, Number obj_value,
    const IpoptData *ip_data, IpoptCalculatedQuantities *ip_cq)
{
  m_Problem->SetVariableValues(x);
  m_Problem->SetLambdaValues(lambda);
}