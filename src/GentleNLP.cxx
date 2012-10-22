#include "GentleNLP.h"
#include <cassert>

namespace gnlp
{

void Problem::AddChildExpression(Expression *ex)
{
  m_ChildExpressions.insert(ex);
}

Expression *Problem::GetPartialDerivative(Expression *ex, Variable *v)
{
  PartialDerivative pd = std::make_pair(ex, v);
  PartialDerivativeMap::const_iterator it = m_Partials.find(pd);
  if(it != m_Partials.end())
    {
    return it->second;
    }
  else
    {
    return m_Partials[pd] =
        ex->DependsOn(v) ?
          ex->MakePartialDerivative(v) : NULL;
    }
}

Problem::~Problem()
{
  for(ExpressionSet::iterator it = m_ChildExpressions.begin();
      it != m_ChildExpressions.end(); it++)
    {
    delete *it;
    }
}

Expression::Expression(Problem *problem)
{
  // Register so the problem can manage memory of this expression
  m_Problem = problem;
  m_Problem->AddChildExpression(this);
}

Expression * Variable::MakePartialDerivative(Variable *variable)
{
  // If the derivative is zero, this should not be created!
  return (variable == this) ? new Constant(m_Problem, 1) : NULL;
}

Expression *Constant::MakePartialDerivative(Variable *)
{
  return NULL;
}

std::string Constant::GetName()
{
  std::ostringstream oss;
  oss << m_Value;
  return oss.str();
}

Expression *NegateOperatorTraits::Differentiate(
    Problem *p, Expression *self, Expression *a, Expression *dA)
{
  return new Negation(p, dA);
}

std::string NegateOperatorTraits::GetName(Expression *a)
{
  std::ostringstream oss;
  oss << "-(" << a->GetName() << ")";
  return oss.str();
}

Expression *SquareOperatorTraits::Differentiate(
    Problem *p, Expression *self, Expression *a, Expression *dA)
{
  return new BinaryProduct(p, new BinaryProduct(p, a, dA), new Constant(p, 2));
}

std::string SquareOperatorTraits::GetName(Expression *a)
{
  std::ostringstream oss;
  oss << "(" << a->GetName() << ")^2";
  return oss.str();
}

Expression *SquareRootOperatorTraits::Differentiate(
    Problem *p, Expression *self, Expression *a, Expression *dA)
{
  if(dA)
    return new BinaryFraction(p,
                              new BinaryProduct(p, new Constant(p, 0.5), dA),
                              self);
  else return NULL;
}

std::string SquareRootOperatorTraits::GetName(Expression *a)
{
  std::ostringstream oss;
  oss << "sqrt(" << a->GetName() << ")";
  return oss.str();
}

Expression *PlusOperatorTraits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *dA, Expression *dB)
{
  if(dA && dB)
    return new BinarySum(p, dA, dB);
  if(dA) return dA;
  if(dB) return dB;
  return NULL;
}

std::string PlusOperatorTraits::GetName(Expression *a, Expression *b)
{
  std::ostringstream oss;
  oss << "(" << a->GetName() << ") + (" << b->GetName() << ")";
  return oss.str();
}

Expression *MinusOperatorTraits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *dA, Expression *dB)
{
  if(dA && dB)
    return new BinaryDifference(p, dA, dB);
  if(dA) return dA;
  if(dB) return new Negation(p, dB);
  return NULL;
}

std::string MinusOperatorTraits::GetName(Expression *a, Expression *b)
{
  std::ostringstream oss;
  oss << "(" << a->GetName() << ") - (" << b->GetName() << ")";
  return oss.str();
}

Expression *ProductOperatorTraits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *dA, Expression *dB)
{
  if(dA && dB)
    return new BinarySum(p,
                         new BinaryProduct(p, a, dB),
                         new BinaryProduct(p, b, dA));
  if(dA) return new BinaryProduct(p, b, dA);
  if(dB) return new BinaryProduct(p, a, dB);
  return NULL;
}

std::string ProductOperatorTraits::GetName(Expression *a, Expression *b)
{
  std::ostringstream oss;
  oss << "(" << a->GetName() << ") * (" << b->GetName() << ")";
  return oss.str();
}

Expression *RatioOperatorTraits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *dA, Expression *dB)
{
  if(dB)
    {
    // Q = b'(a/b)
    Expression *Q = new BinaryProduct(p, dB, self);

    // (a' - Q) / b
    if(dA)
      return new BinaryFraction(p, new BinaryDifference(p, dA, Q), b);
    else
      return new BinaryFraction(p, new Negation(p, Q), b);
    }
  else if(dA)
    return new BinaryFraction(p, dA, b);
  return NULL;
}

std::string RatioOperatorTraits::GetName(Expression *a, Expression *b)
{
  std::ostringstream oss;
  oss << "(" << a->GetName() << ") / (" << b->GetName() << ")";
  return oss.str();
}

Expression *MakeSum(Problem *p, std::vector<Expression *> &expr)
{
  switch(expr.size())
    {
    case 0: return NULL;
    case 1: return expr[0];
    case 2: return new BinarySum(p, expr[0], expr[1]);
    case 3: return new TernarySum(p, expr[0], expr[1], expr[2]);
    }

  BigSum *bs = new BigSum(p);
  for(int i = 0; i < expr.size(); i++)
    bs->AddSummand(expr[i]);

  return bs;
}


Expression *GradientMagnitude3Traits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *c,
    Expression *dA, Expression *dB, Expression *dC)
{
  // Ignore trivial case
  if(!dA && !dB && !dC)
    return NULL;

  // Compute the numerator pieces
  std::vector<Expression *> N;
  if(dA) N.push_back(new BinaryProduct(p, a, dA));
  if(dB) N.push_back(new BinaryProduct(p, b, dB));
  if(dC) N.push_back(new BinaryProduct(p, c, dC));

  // Create the numerator
  Expression *num = MakeSum(p, N);

  // Create the expression
  return new BinaryFraction(p, num, self);
}

std::string GradientMagnitude3Traits::GetName(
    Expression *a, Expression *b, Expression *c)
{
  std::ostringstream oss;
  oss << "|" << a->GetName() << "," << b->GetName() << "," << c->GetName() << "|";
  return oss.str();
}


Expression *GradientMagnitudeSqr3Traits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *c,
    Expression *dA, Expression *dB, Expression *dC)
{
  // Ignore trivial case
  if(!dA && !dB && !dC)
    return NULL;

  // Compute the numerator pieces
  std::vector<Expression *> N;
  if(dA) N.push_back(new BinaryProduct(p, a, dA));
  if(dB) N.push_back(new BinaryProduct(p, b, dB));
  if(dC) N.push_back(new BinaryProduct(p, c, dC));

  // Create the numerator
  Expression *num = MakeSum(p, N);

  // Create the expression
  return new BinaryProduct(p, new Constant(p, 2), num);
}

std::string GradientMagnitudeSqr3Traits::GetName(
    Expression *a, Expression *b, Expression *c)
{
  std::ostringstream oss;
  oss << "|" << a->GetName() << "," << b->GetName() << "," << c->GetName() << "|^2";
  return oss.str();
}


Expression *SumOperator3Traits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *c,
    Expression *dA, Expression *dB, Expression *dC)
{
  // Ignore trivial case
  if(!dA && !dB && !dC)
    return NULL;

  // Combine the non-null pieces
  std::vector<Expression *> N;
  if(dA) N.push_back(dA);
  if(dB) N.push_back(dB);
  if(dC) N.push_back(dC);

  // Create the numerator
  return MakeSum(p, N);
}

std::string SumOperator3Traits::GetName(
    Expression *a, Expression *b, Expression *c)
{
  std::ostringstream oss;
  oss << "(" << a->GetName() << ") + (" << b->GetName() << ") + (" << c->GetName() << ")";
  return oss.str();
}


Expression *ProductOperator3Traits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *c,
    Expression *dA, Expression *dB, Expression *dC)
{
  // Ignore trivial case
  if(!dA && !dB && !dC)
    return NULL;

  // Combine the non-null pieces
  std::vector<Expression *> N;
  if(dA) N.push_back(new TernaryProduct(p, dA, b, c));
  if(dB) N.push_back(new TernaryProduct(p, a, dB, c));
  if(dC) N.push_back(new TernaryProduct(p, a, b, dC));

  // Create the numerator
  return MakeSum(p, N);
}

std::string ProductOperator3Traits::GetName(
    Expression *a, Expression *b, Expression *c)
{
  std::ostringstream oss;
  oss << "(" << a->GetName() << "+" << b->GetName() << "+" << c->GetName() << ")";
  return oss.str();
}

Expression *BigSum::MakePartialDerivative(Variable *variable)
{
  std::vector<Expression *> pd;
  for(Iterator it = m_A.begin(); it != m_A.end(); it++)
    {
    if((*it)->DependsOn(variable))
      pd.push_back(m_Problem->GetPartialDerivative(*it, variable));
    }

  return MakeSum(m_Problem, pd);
}

std::string BigSum::GetName()
{
  std::ostringstream oss;
  for(Iterator it = m_A.begin(); it != m_A.end(); it++)
    {
    if(it != m_A.begin())
      oss << " + ";
    oss << "(" << (*it)->GetName() << ")";
    }
  return oss.str();
}

/******************************************************************
  Generic Constrained PRoblem STUFF
  *****************************************************************/
const double ConstrainedNonLinearProblem::LBINF = -1e100;
const double ConstrainedNonLinearProblem::UBINF = +1e100;

ConstrainedNonLinearProblem::ConstrainedNonLinearProblem()
{
  m_SigmaF = new Variable(this, "SigmaF");
  m_SigmaF->SetValue(1.0);
}

void ConstrainedNonLinearProblem
::AddVariable(Variable *var, double cmin, double cmax)
{
  // Store the index in the variable for later use
  var->SetIndex(m_X.size());

  // Store the variable in the array
  m_X.push_back(var);
  m_LowerBoundX.push_back(cmin);
  m_UpperBoundX.push_back(cmax);
}

Variable *ConstrainedNonLinearProblem
::AddVariable(std::string name, double val, double cmin, double cmax)
{
  Variable *v = new Variable(this, name);
  this->AddVariable(v, cmin, cmax);
  v->SetValue(val);
  return v;
}

void ConstrainedNonLinearProblem
::AddConstraint(Expression *ex, double cmin, double cmax)
{
  // Store the constraint
  m_G.push_back(ex);
  m_LowerBoundG.push_back(cmin);
  m_UpperBoundG.push_back(cmax);

  // Create a lambda for this constraint
  std::string name = "Lambda "; name += m_Lambda.size();
  Variable *lambda = new Variable(this, name);
  lambda->SetValue(1.0);
  m_Lambda.push_back(lambda);
}

void ConstrainedNonLinearProblem
::SetObjective(Expression *ex)
{
  m_F = ex;
}

void ConstrainedNonLinearProblem
::SetupProblem()
{
  typedef std::set<Variable *> VarSet;

  // Only a single initialization is allowed!
  assert(m_GradF.size() == 0);

  // This is where all the work takes place. This function is responsible to
  // set up the Jacobian and Hessian matrices

  // Start with the gradient
  std::cout << "Computing partial derivatives of the objective" << std::endl;
  for(int i = 0; i < m_X.size(); i++)
    {
    m_GradF.push_back(this->GetPartialDerivative(m_F, m_X[i]));
    std::cout << "." << std::flush;
    if((i+1) % 50 == 0)
      std::cout << " " << (i+1) << std::endl;
    }
  std::cout << std::endl;

  // Now the Jacobian of G
  SparseExpressionMatrix::STLSourceType stl_DG;
  for(int j = 0; j < m_G.size(); j++)
    {
    // Find the dependent variables of the constraint
    VarSet vars;
    m_G[j]->GetDependentVariables(vars);

    // Create a row for the sparse matrix
    SparseExpressionMatrix::STLRowType stl_row;

    // For each dependent variable, create a sparse entry.
    for(VarSet::iterator it = vars.begin(); it != vars.end(); it++)
      {
      Variable *v = *it;
      Expression *dg_dv = this->GetPartialDerivative(m_G[j], v);
      stl_row.push_back(std::make_pair(v->GetIndex(), dg_dv));
      }

    // Add the row
    stl_DG.push_back(stl_row);
    }

  // Form a proper sparse matrix from the STL construct
  m_DG.SetFromSTL(stl_DG, m_X.size());

  // Now, for the trickiest part, we need the Hessian! To build it up, we
  // need a different data structure: a map from row/col pair to Expression
  typedef std::pair<int, int> IndexPair;
  typedef std::map<IndexPair, Expression *> SparseMap;
  SparseMap Hmap;

  // We begin by adding entries from the objective function
  std::cout << "Computing Hessian entries" << std::endl;
  for(int i = 0; i < m_X.size(); i++)
    {
    Expression *df = m_GradF[i];
    if(df)
      {
      VarSet depvar;
      df->GetDependentVariables(depvar);
      for(VarSet::iterator it = depvar.begin(); it != depvar.end(); it++)
        {
        Variable *v = *it;
        if(v->GetIndex() <= i)   // Hessian must be lower triangular!
          {
          Expression *ddf = this->GetPartialDerivative(df, v);
          IndexPair idx = std::make_pair(i, v->GetIndex());
          BigSum *sum = new BigSum(this);
          sum->AddSummand(new BinaryProduct(this, m_SigmaF, ddf));
          Hmap[idx] = sum;
          }
        }
      }
    std::cout << "." << std::flush;
    if((i+1) % 50 == 0)
      std::cout << " " << (i+1) << std::endl;
    }
  std::cout << std::endl;

  // Now add entries for each of the constraints
  for(int j = 0; j < m_G.size(); j++)
    {
    // Iterate over all partial derivatives of G
    for(SparseExpressionMatrix::RowIterator rit = m_DG.Row(j); !rit.IsAtEnd(); ++rit)
      {
      Expression *dg = rit.Value();
      int row = rit.Column(); // Column in DG = Row in the hessian matrix
      VarSet depvar;
      dg->GetDependentVariables(depvar);

      for(VarSet::iterator it = depvar.begin(); it != depvar.end(); it++)
        {
        Variable *v = *it;
        if(v->GetIndex() <= row)   // Hessian must be lower triangular!
          {
          Expression *ddg = this->GetPartialDerivative(dg, v);
          IndexPair idx = std::make_pair(row, v->GetIndex());
          BigSum *sum;
          SparseMap::iterator mit = Hmap.find(idx);
          if(mit != Hmap.end())
            {
            sum = (BigSum *) mit->second;
            }
          else
            {
            sum = new BigSum(this);
            Hmap[idx] = sum;
            }

          sum->AddSummand(new BinaryProduct(this, m_Lambda[j], ddg));
          }
        }
      }
    }

  // Now we need to populate a sparse matrix
  SparseExpressionMatrix::STLSourceType stl_H;
  stl_H.resize(m_X.size());
  for(SparseMap::iterator it = Hmap.begin(); it != Hmap.end(); it++)
    {
    int row = it->first.first;
    int col = it->first.second;
    Expression *ex = it->second;
    stl_H[row].push_back(std::make_pair(col, ex));
    }

  m_Hessian.SetFromSTL(stl_H, m_X.size());
}

void ConstrainedNonLinearProblem::SetLambdaValues(const double *l)
{
  for(int i = 0; i < m_Lambda.size(); i++)
    m_Lambda[i]->SetValue(l[i]);
}

void ConstrainedNonLinearProblem::GetVariableBounds(int i, double &lb, double &ub)
{
  lb = m_LowerBoundX[i];
  ub = m_UpperBoundX[i];
}

void ConstrainedNonLinearProblem::SetVariableValues(const double *x)
{
  for(int i = 0; i < m_X.size(); i++)
    m_X[i]->SetValue(x[i]);
}

void ConstrainedNonLinearProblem::GetConstraintBounds(int i, double &lb, double &ub)
{
  lb = m_LowerBoundG[i];
  ub = m_UpperBoundG[i];
}

unsigned int ConstrainedNonLinearProblem::GetNumberOfVariables()
{
  return m_X.size();
}

unsigned int ConstrainedNonLinearProblem::GetNumberOfConstraints()
{
  return m_Lambda.size();
}

void ConstrainedNonLinearProblem::SetSigma(double value)
{
  m_SigmaF->SetValue(value);
}


}





