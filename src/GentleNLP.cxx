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
  // Combined find/insert command
  std::pair<DependencyMap::iterator, bool> ret =
      m_DependencyMap.insert(std::make_pair(ex, Dependency()));

  if(ret.second)
    {
    // Insertion occurred, need to set the RHS to real value
    ex->MakeDependencyList(ret.first->second);
    }

  return ret.first->second;
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

  // Combined find/insert command
  std::pair<PartialDerivativeMap::iterator, bool> ret =
      m_Partials.insert(std::make_pair(pd, static_cast<Expression *>(NULL)));

  if(ret.second)
    {
    // Insertion occurred, need to set the RHS to real value
    const Dependency &dep = this->GetDependentVariables(ex);
    if(dep.find(v) != dep.end())
      ret.first->second = ex->MakePartialDerivative(v);
    }

  // Return the partial derivative expression
  return ret.first->second;

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

void Problem::PrintMemoryStatistics()
{
  // Print the size of all expressions

  // Total size of all dependencies
  unsigned long dep_size;
  for(DependencyMap::iterator jt = m_DependencyMap.begin();
      jt != m_DependencyMap.end(); jt++)
    {
    dep_size += jt->second.size();
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





// Make a product of N expressions
Expression *MakeProduct(Problem *p, int n, Expression **src)
{
  assert(n <= 3);
  double partConst = 1.0;
  std::vector<Expression *> partExpr;

  for(int i = 0; i < n; i++)
    {
    if(src[i]->IsConstant())
      partConst *= src[i]->Evaluate();
    else
      partExpr.push_back(src[i]);
    }

  if(partExpr.size() == 0)
    {
    return new Constant(p, partConst);
    }

  else if(partExpr.size() == 1)
    {
    if(partConst == 1.0)
      return partExpr[0];
    else
      return new ScalarProduct(p, partExpr[0], partConst);
    }

  else if(partExpr.size() == 2)
    {
    if(partConst == 1.0)
      return new BinaryProduct(p, partExpr[0], partExpr[1]);
    else
      return new ScalarProduct(p, new BinaryProduct(p, partExpr[0], partExpr[1]), partConst);
    }

  else
    {
    return new TernaryProduct(p, partExpr[0], partExpr[1], partExpr[2]);
    }
}

Expression* MakeProduct(Problem *p, Expression *a, Expression *b)
{
  Expression *ex[] = {a, b};
  return MakeProduct(p, 2, ex);
}

Expression* MakeProduct(Problem *p, Expression *a, Expression *b, Expression *c)
{
  Expression *ex[] = {a, b, c};
  return MakeProduct(p, 3, ex);
}

Expression* MakeProduct(Problem *p, Expression *a, double scalar)
{
  if(a->IsConstant())
    return new Constant(p, a->Evaluate() * scalar);
  else
    return new ScalarProduct(p, a, scalar);
}

Expression* MakeNegation(Problem *p, Expression *a)
{
  if(a->IsConstant())
    return new Constant(p, -a->Evaluate());
  else
    return new Negation(p, a);
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
  return MakeNegation(p, dA);
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
  return MakeProduct(p, a, dA, new Constant(p, 2.0));
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
    {
    return new BinaryFraction(p, MakeProduct(p, dA, 0.5), self);
    }
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
  WeightedSumGenerator wsg(p);
  if(dA) wsg.AddTerm(dA);
  if(dB) wsg.AddTerm(dB);
  return wsg.GenerateSum();
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
  WeightedSumGenerator wsg(p);
  if(dA) wsg.AddTerm(dA);
  if(dB) wsg.AddTerm(dB, -1);
  return wsg.GenerateSum();
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
  // We should be careful about cases where the derivatives of some
  // expressions are constants, since we can fold things into a
  // simpler expression
  WeightedSumGenerator wsg(p);
  if(dA) wsg.AddTerm(MakeProduct(p, dA, b));
  if(dB) wsg.AddTerm(MakeProduct(p, dB, a));
  return wsg.GenerateSum();
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
    Expression *Q = MakeProduct(p, dB, self);

    // (a' - Q) / b
    if(dA)
      return new BinaryFraction(p, new BinaryDifference(p, dA, Q), b);
    else
      return new BinaryFraction(p, MakeNegation(p, Q), b);
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
  // Some simplification to avoid expressions with lots of constants
  Expression *dA = m_Problem->GetPartialDerivative(m_A, variable);
  if(dA)
    {
    if(dA->IsConstant())
      return new Constant(m_Problem, m_Const * dA->Evaluate());
    else
      return new ScalarProduct(m_Problem, dA, m_Const);
    }
  else
    return NULL;
}

Expression *MakeSum(Problem *p, std::vector<Expression *> &expr)
{
  // We want to find all the unique expressions and all the constants
  typedef std::map<Expression *, double> Map;
  Map exmap;

  // Group expressions by recurrence
  for(int i = 0; i < expr.size(); i++)
    {
    Expression *e = expr[i];
    if(e->IsConstant())
      {
      Map::iterator it = exmap.find(NULL);
      if(it != exmap.end())
        it->second += e->Evaluate();
      else
        exmap[NULL] = e->Evaluate();
      }
    else
      {
      Map::iterator it = exmap.find(e);
      if(it != exmap.end())
        it->second++;
      else
        exmap[e] = 1.0;
      }
    }

  // Add up the expressions
  std::vector<Expression *> compact;
  for(Map::iterator it = exmap.begin(); it!=exmap.end(); ++it)
    {
    if(it->second != 0.0)
      {
      if(it->first == NULL)
        {
        compact.push_back(new Constant(p, it->second));
        }
      else if(it->second == 1.0)
        {
        compact.push_back(it->first);
        }
      else
        {
        compact.push_back(new ScalarProduct(p, it->first, it->second));
        }
      }
    }

  // Now turn into a sum
  switch(compact.size())
    {
    case 0: return NULL;
    case 1: return compact[0];
    case 2: return new BinarySum(p, compact[0], compact[1]);
    case 3: return new TernarySum(p, compact[0], compact[1], compact[2]);
    }

  BigSum *bs = new BigSum(p);
  for(int i = 0; i < compact.size(); i++)
    bs->AddSummand(compact[i]);

  return bs;
}

/*
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
*/

Expression *GradientMagnitude3Traits::Differentiate(
    Problem *p, Expression *self,
    Expression *a, Expression *b, Expression *c,
    Expression *dA, Expression *dB, Expression *dC)
{
  // Ignore trivial case
  if(!dA && !dB && !dC)
    return NULL;

  // Compute the numerator pieces
  WeightedSumGenerator wsg(p);
  if(dA) wsg.AddTerm(MakeProduct(p, a, dA));
  if(dB) wsg.AddTerm(MakeProduct(p, b, dB));
  if(dC) wsg.AddTerm(MakeProduct(p, c, dC));

  // Create the numerator
  Expression *num = wsg.GenerateSum();

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
  WeightedSumGenerator wsg(p);
  if(dA) wsg.AddTerm(MakeProduct(p, a, dA));
  if(dB) wsg.AddTerm(MakeProduct(p, b, dB));
  if(dC) wsg.AddTerm(MakeProduct(p, c, dC));

  // Create the numerator
  Expression *num = wsg.GenerateSum();

  // Create the expression
  return MakeProduct(p, num, 2.0);
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
  WeightedSumGenerator wsg(p);
  if(dA) wsg.AddTerm(dA);
  if(dB) wsg.AddTerm(dB);
  if(dC) wsg.AddTerm(dC);

  // Create the numerator
  return wsg.GenerateSum();
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
  WeightedSumGenerator wsg(p);
  if(dA) wsg.AddTerm(MakeProduct(p, dA, b, c));
  if(dB) wsg.AddTerm(MakeProduct(p, a, dB, c));
  if(dC) wsg.AddTerm(MakeProduct(p, a, b, dC));

  // Create the numerator
  return wsg.GenerateSum();
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
  WeightedSumGenerator wsg(m_Problem);
  for(Iterator it = m_A.begin(); it != m_A.end(); it++)
    {
    Expression *pd_i = m_Problem->GetPartialDerivative(*it, variable);
    if(pd_i)
      wsg.AddTerm(pd_i);
    }

  return wsg.GenerateSum();
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

Expression *MagnitudeSqr(Problem *p, const VarVec &a)
{
  WeightedSumGenerator wsg(p);
  for(int i = 0; i < a.size(); i++)
    wsg.AddTerm(new Square(p, a[i]));
  return wsg.GenerateSum();
}

Expression *DistanceSqr(Problem *p, const VarVec &a, const VarVec &b)
{
  WeightedSumGenerator wsg(p);
  for(int i = 0; i < a.size(); i++)
    wsg.AddTerm(new Square(p, new BinaryDifference(p, a[i], b[i])));
  return wsg.GenerateSum();
}

Expression *DotProduct(Problem *p, const VarVec &a, const VarVec &b)
{
  WeightedSumGenerator wsg(p);
  for(int i = 0; i < a.size(); i++)
    wsg.AddTerm(MakeProduct(p, a[i], b[i]));
  return wsg.GenerateSum();
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
::AddConstraint(Expression *ex, const std::string &category, double cmin, double cmax)
{
  // Store the constraint
  m_G.push_back(ex);
  m_LowerBoundG.push_back(cmin);
  m_UpperBoundG.push_back(cmax);
  m_ConstraintCategory.push_back(category);

  // Create a lambda for this constraint
  std::string name = "Lambda "; name += m_Lambda.size();
  Variable *lambda = new Variable(this, name);
  lambda->SetValue(1.0);
  m_Lambda.push_back(lambda);
}

Variable* ConstrainedNonLinearProblem
::AddExpressionAsConstrainedVariable(Expression *exp, const std::string &category)
{
  Variable *v = this->AddVariable(category, exp->Evaluate());
  this->AddConstraint(new BinaryDifference(this, v, exp), category, 0, 0);
  return v;
}

void ConstrainedNonLinearProblem
::SetObjective(Expression *ex)
{
  m_F = ex;
}

struct counter {
  int value;
  operator int () const { return value; }
  counter &operator ++() { value ++; return *this; }
  counter() : value(0) {}
};

double TestCentralDifference(Expression *ex, Variable *v, double delta=1e-5)
{
  double val = v->Evaluate();

  // Do the perturbation
  ex->MakeTreeDirty();
  v->SetValue(val + delta);
    double f2 = ex->Evaluate();

  ex->MakeTreeDirty();
  v->SetValue(val - delta);
  double f1 = ex->Evaluate();

  v->SetValue(val);

  return (f2 - f1) / (2 * delta);
}

void TestDeriv(Problem *p, Expression *ex, Variable *v, Expression *deriv, std::string what)
{
  if(deriv) deriv->MakeTreeDirty();
  double dAnalytic = deriv ? deriv->Evaluate() : 0.0;
  double dNumeric = TestCentralDifference(ex, v);
  double diffAbs = fabs(dAnalytic - dNumeric);
  double diffRel = 2 * (dAnalytic - dNumeric) / fabs(dAnalytic + dNumeric);
  if(diffAbs > 1e-6)
    {
    std::cout << "--- BAD DERIVATIVE in " << what << " --- " << std::endl;
    std::cout << "EXPRESSION: " << ex->GetName() << std::endl;
    std::cout << "VARIABLE:" << v->GetName() << std::endl;
    std::cout << "DERIVATIVE" << v->GetName() << std::endl;
    std::cout << "NUMERIC: "<< dNumeric << std::endl;
    std::cout << "ANALYTC: "<< dAnalytic << std::endl;
    std::cout << "DIFFERN: "<< diffAbs << std::endl;
    }
}

void ConstrainedNonLinearProblem
::SetupProblem(bool hessian, bool deriv_test)
{
  typedef std::set<Variable *> VarSet;

  // Only a single initialization is allowed!
  m_GradF.clear();

  // This is where all the work takes place. This function is responsible to
  // set up the Jacobian and Hessian matrices

  // Start with the gradient
  for(int i = 0; i < m_X.size(); i++)
    {
    m_GradF.push_back(this->GetPartialDerivative(m_F, m_X[i]));
    if(deriv_test)
      TestDeriv(this, m_F, m_X[i], m_GradF.back(), "Objective");
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
      if(dg_dv)
        stl_row.push_back(std::make_pair(v->GetIndex(), dg_dv));
      if(deriv_test)
        TestDeriv(this, m_G[j], v, dg_dv, "Constraints");
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

  // Keep track of hessian entries
  double nSecondDeriv = 0;
  std::map<int, counter> histDepth, histSize;

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
          IndexPair idx = std::make_pair(i, v->GetIndex());
          BigSum *sum = new BigSum(this);
          Hmap[idx] = sum;

          Expression *ddf = this->GetPartialDerivative(df, v);

          if(deriv_test)
            TestDeriv(this, df, v, ddf, "Hessian-F");

          if(ddf)
            {
            sum->AddSummand(MakeProduct(this, ddf, m_SigmaF));

            // Do some stats on the hessian entries
            nSecondDeriv++;
            int depth, size;
            ddf->GetTreeStats(depth, size);
            ++histDepth[depth];
            ++histSize[size];
            if(depth > 1)
              {
              std::cout << ddf->GetName() << std::endl;
              }
            }
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
          IndexPair idx = std::make_pair(row, v->GetIndex());
          std::pair<SparseMap::iterator, bool> ret =
              Hmap.insert(std::make_pair(idx, (BigSum *) NULL));

          if(ret.second)
            {
            // Actual insertion took place
            ret.first->second = new BigSum(this);
            }

          BigSum *sum = (BigSum *) ret.first->second;

          Expression *ddg = this->GetPartialDerivative(dg, v);
          if(deriv_test)
            TestDeriv(this, dg, v, ddg, "Hessian-G");
          if(ddg)
            {
            // Do some stats on the hessian entries
            nSecondDeriv++;
            int depth, size;
            ddg->GetTreeStats(depth, size);
            ++histDepth[depth];
            ++histSize[size];

            sum->AddSummand(MakeProduct(this, m_Lambda[j], ddg));

            if(depth > 1)
              {
              std::cout << "DERIV ISSUE" << std::endl;
              std::cout << "Con = " << m_G[j]->GetName() << std::endl;
              std::cout << "Var1 =  " << m_X[row]->GetName() << std::endl;
              std::cout << "Var2 =  " << v->GetName() << std::endl;
              std::cout << ddg->GetName() << std::endl;
              }
            }
          }
        }
      }
    if((j+1) % 50 == 0)
      std::cout << "." << std::flush;
    if((j+1) % 4000 == 0)
      std::cout << " " << (j+1) << std::endl;
    }

  std::cout << " " << std::endl;
  std::cout << "Hessian Histogram" << std::endl;
  std::cout << "  total entries = " << nSecondDeriv << std::endl;
  std::cout << "  depth histo: " << std::endl;
  for(std::map<int, counter>::const_iterator ith = histDepth.begin(); ith!=histDepth.end(); ith++)
    std::cout << "\t" << ith->first << "\t" << (int) ith->second << std::endl;
  std::cout << "  size histo: " << std::endl;
  for(std::map<int, counter>::const_iterator ith = histSize.begin(); ith!=histSize.end(); ith++)
    std::cout << "\t" << ith->first << "\t" << (int) ith->second << std::endl;

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

WeightedSumGenerator::WeightedSumGenerator(Problem *p)
{
  m_Problem = p;
  m_Constant = 0.0;
}

void WeightedSumGenerator::AddConstant(double value)
{
  m_Constant += value;
}

void WeightedSumGenerator::AddTerm(Expression *expr, double weight)
{

  // Constant expressions get added up
  if(expr->IsConstant())
    {
    m_Constant += expr->Evaluate() * weight;
    return;
    }

  // Negations
  Negation *neg = dynamic_cast<Negation *>(expr);
  if(neg)
    {
    this->AddTerm(neg->GetOperand(0), -weight);
    return;
    }

  // Scalar products
  ScalarProduct *sp = dynamic_cast<ScalarProduct *>(expr);
  if(sp)
    {
    this->AddTerm(sp->GetOperand(0), weight * sp->GetScalar());
    return;
    }

  // Binary sums
  BinarySum *bs = dynamic_cast<BinarySum *>(expr);
  if(bs)
    {
    this->AddTerm(bs->GetOperand(0), weight);
    this->AddTerm(bs->GetOperand(1), weight);
    return;
    }

  // Binary differences
  BinaryDifference *bd = dynamic_cast<BinaryDifference *>(expr);
  if(bd)
    {
    this->AddTerm(bd->GetOperand(0), weight);
    this->AddTerm(bd->GetOperand(1), -weight);
    return;
    }

  // Ternary sums
  TernarySum *ts = dynamic_cast<TernarySum *>(expr);
  if(ts)
    {
    this->AddTerm(ts->GetOperand(0), weight);
    this->AddTerm(ts->GetOperand(1), weight);
    this->AddTerm(ts->GetOperand(2), weight);
    return;
    }

  // Big sums
  BigSum *bigs = dynamic_cast<BigSum *>(expr);
  if(bigs)
    {
    for(int i = 0; i < bigs->GetNumberOfOperands(); i++)
      this->AddTerm(bigs->GetOperand(i), weight);
    return;
    }

  // Append the weight for this expression
  std::pair<WeightMap::iterator, bool> q =
      m_WeightMap.insert(std::pair<Expression *,double>(expr, weight));
  if(!q.second)
    q.first->second += weight;
}

Expression *WeightedSumGenerator::GenerateSum()
{
  Expression *ex = this->DoGenerateSum();
  return ex;
}

Expression *WeightedSumGenerator::DoGenerateSum()
{
  // Check the constant
  if(m_Constant != 0.0)
    {
    Expression *c = new Constant(m_Problem, m_Constant);
    m_WeightMap.insert(std::pair<Expression *, double>(c, 1.0));
    }

  // Drop all the zero cases
  for(WeightMap::iterator it = m_WeightMap.begin(); it != m_WeightMap.end();)
    {
    if(it->second == 0)
      m_WeightMap.erase(it++);
    else
      it++;
    }

  // Handle special cases of 0 terms
  if(m_WeightMap.size() == 0)
    return NULL;

  // Handle special cases of 1 or 2 terms
  else if(m_WeightMap.size() == 1)
    {
    WeightMap::iterator it1 = m_WeightMap.begin();
    if(it1->second == 1.0)
      return it1->first;
    else if(it1->second == -1.0)
      return new Negation(m_Problem, it1->first);
    else
      return new ScalarProduct(m_Problem, it1->first, it1->second);
    }

  else if(m_WeightMap.size() == 2)
    {
    WeightMap::iterator it1 = m_WeightMap.begin();
    WeightMap::iterator it2 = m_WeightMap.begin(); ++it2;

    // Map to scalar product if necessary
    Expression *we1 = it1->second == 1.0 || it1->second == -1.0
        ? it1->first : new ScalarProduct(m_Problem, it1->first, it1->second);

    // Map to scalar product if necessary
    Expression *we2 = it2->second == 1.0 || it2->second == -1.0
        ? it2->first : new ScalarProduct(m_Problem, it2->first, it2->second);

    // Handle -1 cases
    if(it1->second == -1.0 && it2->second == -1.0)
      return new Negation(m_Problem, new BinarySum(m_Problem, we1, we2));
    else if(it1->second == -1.0)
      return new BinaryDifference(m_Problem, we2, we1);
    else if(it2->second == -1.0)
      return new BinaryDifference(m_Problem, we1, we2);
    else
      return new BinarySum(m_Problem, we1, we2);
    }

  else if(m_WeightMap.size() == 3)
    {
    // Here there is only a triple sum to worry about
    WeightMap::iterator it1 = m_WeightMap.begin();
    WeightMap::iterator it2 = m_WeightMap.begin(); ++it2;
    WeightMap::iterator it3 = m_WeightMap.begin(); ++it3; ++it3;

    // Map to scalar product if necessary
    Expression *we1 = it1->second == 1.0
        ? it1->first : new ScalarProduct(m_Problem, it1->first, it1->second);

    // Map to scalar product if necessary
    Expression *we2 = it2->second == 1.0
        ? it2->first : new ScalarProduct(m_Problem, it2->first, it2->second);

    // Map to scalar product if necessary
    Expression *we3 = it3->second == 1.0
        ? it3->first : new ScalarProduct(m_Problem, it3->first, it3->second);

    return new TernarySum(m_Problem, we1, we2, we3);
    }

  else
    {
    BigSum *bs = new BigSum(m_Problem);
    for(WeightMap::iterator it = m_WeightMap.begin(); it != m_WeightMap.end(); ++it)
      {
      Expression *we = it->second == 1.0
          ? it->first : new ScalarProduct(m_Problem, it->first, it->second);
      bs->AddSummand(we);
      }
    return bs;
    }
}




}





