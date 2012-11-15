#include "GentleNLP.h"
#include <cassert>

namespace gnlp
{

void Problem::AddChildExpression(Expression *ex)
{
  m_ChildExpressions.insert(ex);
}

const Problem::Dependency &
Problem::GetDependentVariables(Expression *ex)
{
  DependencyMap::const_iterator it = m_DependencyMap.find(ex);
  if(it != m_DependencyMap.end())
    {
    return it->second;
    }
  else
    {
    // Ask the expression to collect the dependent variables
    ex->MakeDependencyList(m_DependencyMap[ex]);
    return m_DependencyMap[ex];
    }
}

void
Problem::AppendDependentVariables(Expression *ex, Dependency &target)
{
  const Dependency &dep = this->GetDependentVariables(ex);
  target.insert(dep.begin(), dep.end());
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
    const Dependency &dep = this->GetDependentVariables(ex);
    return m_Partials[pd] =
        dep.find(v) != dep.end() ?
          ex->MakePartialDerivative(v) : NULL;
    }
}

void Problem::MakeChildrenDirty()
{
  for(ExpressionSet::iterator it = m_ChildExpressions.begin();
      it != m_ChildExpressions.end(); it++)
    {
    (*it)->MakeDirty();
    }
}

void Problem::ClearDerivativeCaches()
{
  m_DependencyMap.clear();
  m_Partials.clear();
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
  return new ScalarProduct(p, new BinaryProduct(p, a, dA), 2.0);
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
    return new BinaryFraction(p, new ScalarProduct(p, dA, 0.5), self);
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

Expression *ProductOperatorTraits::DiffProductHelper(
    Problem *p, Expression *dA, Expression *b)
{
  // Compute the term b * dA, simplifying if dA is a constant
  if(!dA)
    return NULL;

  Constant *ctry = dynamic_cast<Constant *>(dA);

  if(ctry && ctry->Evaluate() == 1.0)
    return b;
  else if(ctry)
    return new ScalarProduct(p, b, ctry->Evaluate());
  else
    return new BinaryProduct(p, b, dA);
}

Expression *ProductOperatorTraits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *dA, Expression *dB)
{
  // We should be careful about cases where the derivatives of some
  // expressions are constants, since we can fold things into a
  // simpler expression
  Expression *AdB = DiffProductHelper(p, dB, a), *BdA = DiffProductHelper(p, dA, b);

  if(AdB && BdA)
    return new BinarySum(p,AdB,BdA);
  else if(BdA) return BdA;
  else if(AdB) return AdB;
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


std::string
ScalarProduct::GetName()
{
  std::ostringstream oss;
  oss << "(" << m_Const << " * " << m_A->GetName() << ")";
  return oss.str();
}

Expression *
ScalarProduct::MakePartialDerivative(Variable *variable)
{
  Expression *dA = m_Problem->GetPartialDerivative(m_A, variable);
  return dA ? new ScalarProduct(m_Problem, dA, m_Const) : NULL;
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
  return new ScalarProduct(p, num, 2.0);
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
    Expression *pd_i = m_Problem->GetPartialDerivative(*it, variable);
    if(pd_i)
      pd.push_back(pd_i);
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
  Vector helper functions
  *****************************************************************/
VarVec CrossProduct(Problem *p, const VarVec &a, const VarVec &b)
{
  VarVec out(3, NULL);

  out[0] = new BinaryDifference(p,
                                new BinaryProduct(p, a[1], b[2]),
                                new BinaryProduct(p, a[2], b[1]));

  out[1] = new BinaryDifference(p,
                                new BinaryProduct(p, a[2], b[0]),
                                new BinaryProduct(p, a[0], b[2]));

  out[2] = new BinaryDifference(p,
                                new BinaryProduct(p, a[0], b[1]),
                                new BinaryProduct(p, a[1], b[0]));

  return out;
}

Expression *DotProduct(Problem *p, const VarVec &a, const VarVec &b)
{
  VarVec sum;
  for(int i = 0; i < a.size(); i++)
    sum.push_back(new BinaryProduct(p, a[i], b[i]));
  return MakeSum(p, sum);
}

// Return an expression for (2*A)^2, where A is area of the triangle formed by a,b,c
Expression *TriangleTwiceAreaSqr(Problem *p, const VarVec &a, const VarVec &b, const VarVec &c)
{
  VarVec U = VectorApplyPairwise<BinaryDifference>(p, b, a);
  VarVec V = VectorApplyPairwise<BinaryDifference>(p, c, a);
  VarVec W = CrossProduct(p, U, V);
  return new TernaryGradientMagnitudeSqr(p, W[0], W[1], W[2]);
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

Variable* ConstrainedNonLinearProblem
::AddExpressionAsConstrainedVariable(Expression *exp)
{
  Variable *v = this->AddVariable("dummy", exp->Evaluate());
  this->AddConstraint(new BinaryDifference(this, v, exp), 0, 0);
  return v;
}

void ConstrainedNonLinearProblem
::SetObjective(Expression *ex)
{
  m_F = ex;
}

void ConstrainedNonLinearProblem
::SetupProblem(bool hessian)
{
  typedef std::set<Variable *> VarSet;

  // Only a single initialization is allowed!
  assert(m_GradF.size() == 0);

  // This is where all the work takes place. This function is responsible to
  // set up the Jacobian and Hessian matrices

  // Start with the gradient
  for(int i = 0; i < m_X.size(); i++)
    {
    m_GradF.push_back(this->GetPartialDerivative(m_F, m_X[i]));
    }

  // Now the Jacobian of G
  std::cout << "Computing Jacobian entries" << std::endl;
  SparseExpressionMatrix::STLSourceType stl_DG;
  for(int j = 0; j < m_G.size(); j++)
    {
    // Find the dependent variables of the constraint
    const Dependency &vars = this->GetDependentVariables(m_G[j]);

    // Create a row for the sparse matrix
    SparseExpressionMatrix::STLRowType stl_row;

    // For each dependent variable, create a sparse entry.
    for(Dependency::const_iterator it = vars.begin(); it != vars.end(); it++)
      {
      Variable *v = *it;
      Expression *dg_dv = this->GetPartialDerivative(m_G[j], v);
      stl_row.push_back(std::make_pair(v->GetIndex(), dg_dv));
      }

    // Add the row
    stl_DG.push_back(stl_row);

    if((j+1) % 50 == 0)
      std::cout << "." << std::flush;
    if((j+1) % 4000 == 0)
      std::cout << " " << (j+1) << std::endl;
    }

  std::cout << std::endl;

  // Form a proper sparse matrix from the STL construct
  m_DG.SetFromSTL(stl_DG, m_X.size());

  // Now, for the trickiest part, we need the Hessian! To build it up, we
  // need a different data structure: a map from row/col pair to Expression
  typedef std::pair<int, int> IndexPair;
  typedef std::map<IndexPair, Expression *> SparseMap;
  SparseMap Hmap;

  if(!hessian)
    return;

  // We begin by adding entries from the objective function
  std::cout << "Computing Hessian entries" << std::endl;
  for(int i = 0; i < m_X.size(); i++)
    {
    Expression *df = m_GradF[i];
    if(df)
      {
      const Dependency &depvar = this->GetDependentVariables(df);
      for(Dependency::const_iterator it = depvar.begin(); it != depvar.end(); it++)
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
    if((i+1) % 50 == 0)
      std::cout << "." << std::flush;
    if((i+1) % 4000 == 0)
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

      const Dependency &depvar = this->GetDependentVariables(dg);
      for(Dependency::const_iterator it = depvar.begin(); it != depvar.end(); it++)
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

void ConstrainedNonLinearProblem::SetGradientSmoothingKernel(SparseRealMatrix K)
{
  m_GradientKernel = K;
}


void ConstrainedNonLinearProblem::SetLambdaValues(const double *l)
{
  for(int i = 0; i < m_Lambda.size(); i++)
    m_Lambda[i]->SetValue(l[i]);

  // Make all the expressions dirty
  this->MakeChildrenDirty();
}

void ConstrainedNonLinearProblem::GetVariableBounds(int i, double &lb, double &ub)
{
  lb = m_LowerBoundX[i];
  ub = m_UpperBoundX[i];
}

void ConstrainedNonLinearProblem::SetVariableValues(const double *x)
{
  // Set the variable values
  for(int i = 0; i < m_X.size(); i++)
    m_X[i]->SetValue(x[i]);

  // Make all the expressions dirty
  this->MakeChildrenDirty();
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

  // Make all the expressions dirty
  this->MakeChildrenDirty();
}


}





