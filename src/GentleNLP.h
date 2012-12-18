#ifndef GENTLENLP_H
#define GENTLENLP_H

#include <algorithm>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <sstream>

#include "SparseMatrix.h"

namespace gnlp
{

class Expression;
class Variable;


/**
  A representation of a non-linear problem.
 */
class Problem
{
public:
  typedef std::set<Variable *> Dependency;
  typedef std::set<Expression *> ExpressionSet;

  Problem() {}
  virtual ~Problem();

  /** Get the set of dependent variables for a given registered expression.
      Like the partial derivatives, these sets are cached, so repeated calls
      to this method will be quick. */
  const Dependency &GetDependentVariables(Expression *ex);

  /** Similar to above, but the dependent variables will be added to
      the contents of the list deps passed in by reference */
  void AppendDependentVariables(Expression *ex, Dependency &target);

  /** Get a partial derivative for an expression with repsect to a
    variable. The problem keeps track of these expressions so that
    repeat calls to GetPartialDerivative only create child expressions
    when the derivative has not previously been computed */
  Expression *GetPartialDerivative(Expression *, Variable *);

  /** Register a child expression with the problem. This will make sure the
    child expression is deleted. */
  void AddChildExpression(Expression *);

  /** Dirty all the child expressions */
  void MakeChildrenDirty();

  /** Get the set of all stored expressions */
  const ExpressionSet& GetChildExpressions() { return m_ChildExpressions; }

  /** Clear the caches of partial variables and dependencies. Call this when
    the problem has been fully configured and you want to just evaluate it
    on different datasets */
  void ClearDerivativeCaches();

  void PrintMemoryStatistics();

protected:

  // List of known expressions
  ExpressionSet m_ChildExpressions;

  // List of partial derivatives of known expressions
  typedef std::pair<Expression *, Variable *> PartialDerivative;
  typedef std::map<PartialDerivative, Expression *> PartialDerivativeMap;
  PartialDerivativeMap m_Partials;

  // List of dependencies for each expression
  typedef std::map<Expression *, Dependency> DependencyMap;
  DependencyMap m_DependencyMap;
};

/**
  Abstract expression class
 */
class Expression
{
public:

  /** Create an expression. The memory for the expression is managed by
    the problem passed in to the constructor. Do not delete the expression */
  Expression(Problem *problem);
  virtual ~Expression() {}

  /** Evaluate the expression */
  virtual double Evaluate() = 0;

  /** Populate a list of dependent variables */
  virtual void MakeDependencyList(Problem::Dependency &vars) = 0;

  /** Create a partial derivative of the expression with respect to variable.
      This should only be called once per expression/variable, since each call
      will make a new() expression. */
  virtual Expression *MakePartialDerivative(Variable *variable) = 0;

  /** Debugging purposes : print the expression using LateX */
  virtual std::string GetName() = 0;

  /** Dirty the expression. Does nothing by default */
  virtual void MakeDirty() {}

  /** Dirty the expression and the children */
  virtual void MakeTreeDirty()
  {
    this->MakeDirty();
    for(int i = 0; i < this->GetNumberOfOperands(); i++)
      this->GetOperand(i)->MakeTreeDirty();
  }

  /** Whether the expression is a constant */
  virtual bool IsConstant() = 0;

  /** Compute the statistics of the expression (three depth, size) */
  virtual void GetTreeStats(int &depth, int &size) = 0;

  /** Get the operands of the expression */
  virtual int GetNumberOfOperands() const = 0;
  virtual Expression* GetOperand(int i) const = 0;

protected:

  // The problem that owns this expression
  Problem *m_Problem;
};

/**
  A variable that is optimized over. The variable has a value. It can also
  have an optional name and index. These are used by the problem class and
  are not required for the Variable class itself to operate properly
  */
class Variable : public Expression
{
public:
  Variable(Problem *parent, std::string name)
    : Expression(parent), m_Name(name), m_Index(-1),
      m_Value(0), m_IndexInProblem(-1) {}

  void SetValue(double value) { m_Value = value; }
  std::string GetName() { return m_Name; }

  void SetIndex(int value) { m_Index = value; }
  int GetIndex() { return m_Index; }

  virtual double Evaluate() { return m_Value; }

  virtual Expression *MakePartialDerivative(Variable *variable);

  virtual void MakeDependencyList(Problem::Dependency &vars)
  {
    vars.insert(this);
  }

  virtual bool IsConstant() { return false; }

  virtual void GetTreeStats(int &depth, int &size)
  {
    depth = 1;
    size = 1;
  }

  virtual int GetNumberOfOperands() const { return 0; }
  virtual Expression* GetOperand(int i) const { return NULL; }

protected:
  std::string m_Name;
  int m_Index;
  int m_IndexInProblem;
  double m_Value;
};


/**
  A constant expression (in the context of optimization)
  */
class Constant : public Expression
{
public:
  Constant(Problem *parent, double value = 0.0)
    : Expression(parent), m_Value(value) {}

  virtual double Evaluate() { return m_Value; }
  std::string GetName();
  virtual Expression *MakePartialDerivative(Variable *variable);

  virtual void MakeDependencyList(Problem::Dependency &vars) {}

  virtual bool IsConstant() { return true; }

  virtual void GetTreeStats(int &depth, int &size)
  {
    depth = 1;
    size = 1;
  }

  virtual int GetNumberOfOperands() const { return 0; }
  virtual Expression* GetOperand(int i)  const { return NULL; }

protected:
  double m_Value;
};


class CachingExpression : public Expression
{
public:
  CachingExpression(Problem *parent)
    : Expression(parent), m_Value(0.0), m_Dirty(true) {}

  /** Evaluate returns cached value if possible */
  virtual double Evaluate()
  {
    if(m_Dirty)
      {
      m_Value = this->ComputeValue();
      m_Dirty = false;
      }
    return m_Value;
  }

  /** Do the actual evaluation */
  virtual double ComputeValue() = 0;

  /** Dirty the expression */
  virtual void MakeDirty() { m_Dirty = true; }

  virtual bool IsConstant() { return false; }

protected:

  // The cached value of the expression
  double m_Value;

  // Whether the cached value is valid
  bool m_Dirty;
};

class NegateOperatorTraits
{
public:
  static double Operate(double a) { return -a; }
  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *dA);
  static std::string GetName(Expression *a);
};

class SquareOperatorTraits
{
public:
  static double Operate(double a) { return a * a; }
  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *dA);
  static std::string GetName(Expression *a);
};

class SquareRootOperatorTraits
{
public:
  static double Operate(double a) { return sqrt(a); }
  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *dA);
  static std::string GetName(Expression *a);
};


template<class TOperatorTraits>
class UnaryExpression : public CachingExpression
{
public:
  UnaryExpression(Problem *parent, Expression *a)
    : CachingExpression(parent), m_A(a) {}

  virtual double ComputeValue()
  {
    return TOperatorTraits::Operate(m_A->Evaluate());
  }

  virtual void MakeDependencyList(Problem::Dependency &vars)
  {
    m_Problem->AppendDependentVariables(m_A, vars);
  }

  virtual Expression *MakePartialDerivative(Variable *variable)
  {
    Expression *dA = m_Problem->GetPartialDerivative(m_A, variable);
    return dA ? TOperatorTraits::Differentiate(m_Problem, this, m_A, dA) : NULL;
  }

  std::string GetName() { return TOperatorTraits::GetName(m_A); }

  virtual void GetTreeStats(int &depth, int &size)
  {
    int da, sa;
    m_A->GetTreeStats(da, sa);
    depth = 1 + da;
    size = 1 + sa;
  }

  virtual int GetNumberOfOperands() const { return 1; }
  virtual Expression* GetOperand(int i) const { return i == 0 ? m_A : NULL; }



protected:
  Expression *m_A;
};

typedef UnaryExpression<NegateOperatorTraits> Negation;
typedef UnaryExpression<SquareOperatorTraits> Square;
typedef UnaryExpression<SquareRootOperatorTraits> SquareRoot;




class ScalarProduct : public CachingExpression
{
public:
  ScalarProduct(Problem *parent, Expression *a, double value)
    : CachingExpression(parent), m_A(a), m_Const(value) {}

  virtual double ComputeValue()
  {
    return m_Const * m_A->Evaluate();
  }

  virtual void MakeDependencyList(Problem::Dependency &vars)
  {
    m_Problem->AppendDependentVariables(m_A, vars);
  }

  virtual Expression *MakePartialDerivative(Variable *variable);

  std::string GetName();

  virtual void GetTreeStats(int &depth, int &size)
  {
    int da, sa;
    m_A->GetTreeStats(da, sa);
    depth = 1 + da;
    size = 1 + sa;
  }

  virtual int GetNumberOfOperands() const { return 1; }
  virtual Expression* GetOperand(int i) const { return i == 0 ? m_A : NULL; }

  double GetScalar() { return m_Const; }

protected:
  Expression *m_A;
  double m_Const;
};




class PlusOperatorTraits
{
public:
  static double Operate(double a, double b) { return a + b; }
  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *b,
                                   Expression *dA, Expression *dB);
  static std::string GetName(Expression *a, Expression *b);
};

class MinusOperatorTraits
{
public:
  static double Operate(double a, double b) { return a - b; }
  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *b,
                                   Expression *dA, Expression *dB);
  static std::string GetName(Expression *a, Expression *b);
};

class ProductOperatorTraits
{
public:
  static double Operate(double a, double b) { return a * b; }
  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *b,
                                   Expression *dA, Expression *dB);
  static std::string GetName(Expression *a, Expression *b);
};

class RatioOperatorTraits
{
public:
  static double Operate(double a, double b) { return a / b; }
  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *b,
                                   Expression *dA, Expression *dB);
  static std::string GetName(Expression *a, Expression *b);
};


/**
  A sum, difference, product of two expressions
  */
template <class TOperatorTraits>
class BinaryExpression : public CachingExpression
{
public:
  BinaryExpression(Problem *parent, Expression *a, Expression *b)
    : CachingExpression(parent), m_A(a), m_B(b) {}

  virtual double ComputeValue()
  {
    return TOperatorTraits::Operate(m_A->Evaluate(), m_B->Evaluate());
  }

  virtual void MakeDependencyList(Problem::Dependency &vars)
  {
    m_Problem->AppendDependentVariables(m_A, vars);
    m_Problem->AppendDependentVariables(m_B, vars);
  }

  virtual Expression *MakePartialDerivative(Variable *variable)
  {
    return TOperatorTraits::Differentiate(
          m_Problem, this,
          m_A, m_B,
          m_Problem->GetPartialDerivative(m_A, variable),
          m_Problem->GetPartialDerivative(m_B, variable));
  }

  std::string GetName() { return TOperatorTraits::GetName(m_A, m_B); }


  virtual void GetTreeStats(int &depth, int &size)
  {
    int da, sa, db, sb;
    m_A->GetTreeStats(da, sa);
    m_B->GetTreeStats(db, sb);
    depth = 1 + std::max(da, db);
    size = 1 + sa + sb;
  }

  virtual int GetNumberOfOperands() const { return 2; }
  virtual Expression* GetOperand(int i) const
  {
    switch(i)
      {
      case 0 : return m_A;
      case 1 : return m_B;
      default: return NULL;
      }
  }

protected:
  Expression *m_A, *m_B;

};

typedef BinaryExpression<PlusOperatorTraits> BinarySum;
typedef BinaryExpression<MinusOperatorTraits> BinaryDifference;
typedef BinaryExpression<ProductOperatorTraits> BinaryProduct;
typedef BinaryExpression<RatioOperatorTraits> BinaryFraction;


class GradientMagnitude3Traits
{
public:
  static double Operate(double a, double b, double c)
  {
    return sqrt(a*a + b*b + c*c);
  }

  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *b, Expression *c,
                                   Expression *dA, Expression *dB, Expression *dC);

  static std::string GetName(Expression *a, Expression *b, Expression *c);
};

class GradientMagnitudeSqr3Traits
{
public:
  static double Operate(double a, double b, double c)
  {
    return a*a + b*b + c*c;
  }

  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *b, Expression *c,
                                   Expression *dA, Expression *dB, Expression *dC);

  static std::string GetName(Expression *a, Expression *b, Expression *c);
};

class SumOperator3Traits
{
public:
  static double Operate(double a, double b, double c)
  {
    return a + b + c;
  }

  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *b, Expression *c,
                                   Expression *dA, Expression *dB, Expression *dC);

  static std::string GetName(Expression *a, Expression *b, Expression *c);
};

class ProductOperator3Traits
{
public:
  static double Operate(double a, double b, double c)
  {
    return a * b * c;
  }

  static Expression *Differentiate(Problem *p, Expression *self,
                                   Expression *a, Expression *b, Expression *c,
                                   Expression *dA, Expression *dB, Expression *dC);

  static std::string GetName(Expression *a, Expression *b, Expression *c);
};



/**
  A function of three expressions
  */
template <class TOperatorTraits>
class TernaryExpression : public CachingExpression
{
public:
  TernaryExpression(Problem *parent, Expression *a, Expression *b, Expression *c)
    : CachingExpression(parent), m_A(a), m_B(b), m_C(c) {}

  virtual double ComputeValue()
  {
    return TOperatorTraits::Operate(m_A->Evaluate(), m_B->Evaluate(), m_C->Evaluate());
  }

  virtual void MakeDependencyList(Problem::Dependency &vars)
  {
    m_Problem->AppendDependentVariables(m_A, vars);
    m_Problem->AppendDependentVariables(m_B, vars);
    m_Problem->AppendDependentVariables(m_C, vars);
  }

  virtual Expression *MakePartialDerivative(Variable *variable)
  {
    return TOperatorTraits::Differentiate(
          m_Problem, this,
          m_A, m_B, m_C,
          m_Problem->GetPartialDerivative(m_A, variable),
          m_Problem->GetPartialDerivative(m_B, variable),
          m_Problem->GetPartialDerivative(m_C, variable));
  }

  std::string GetName() { return TOperatorTraits::GetName(m_A, m_B, m_C); }

  virtual void GetTreeStats(int &depth, int &size)
  {
    int da, sa, db, sb, dc, sc;
    m_A->GetTreeStats(da, sa);
    m_B->GetTreeStats(db, sb);
    m_C->GetTreeStats(dc, sc);
    depth = 1 + std::max(dc, std::max(da, db));
    size = 1 + sa + sb + sc;
  }

  virtual int GetNumberOfOperands() const { return 3; }
  virtual Expression* GetOperand(int i) const
  {
    switch(i)
      {
      case 0 : return m_A;
      case 1 : return m_B;
      case 2 : return m_C;
      default: return NULL;
      }
  }


protected:
  Expression *m_A, *m_B, *m_C;
};


typedef TernaryExpression<GradientMagnitude3Traits> TernaryGradientMagnitude;
typedef TernaryExpression<GradientMagnitudeSqr3Traits> TernaryGradientMagnitudeSqr;
typedef TernaryExpression<SumOperator3Traits> TernarySum;
typedef TernaryExpression<ProductOperator3Traits> TernaryProduct;



Expression *MakeSum(Problem *p, std::vector<Expression *> &expr);
Expression *MakeProduct(Problem *p, Expression *e1, Expression *e2);
Expression *MakeProduct(Problem *p, Expression *e1, Expression *e2, Expression *e3);

/**
  This helper function creates a sum of several expressions based on their
  number, i.e., BinarySum, TernarySum or BigSum
  */
// Expression *MakeSum(Problem *p, VarVec &expr);

/**
 * A class that generates an expression for a weighted sum of terms,
 * simplifying as much as possible by pulling together constants and
 * similar expressions
 */
class WeightedSumGenerator
{
public:
  WeightedSumGenerator(Problem *p);

  void AddTerm(Expression *expr, double weight = 1.0);
  void AddConstant(double value);

  Expression *GenerateSum();

protected:
  Problem *m_Problem;

  Expression *DoGenerateSum();

  typedef std::map<Expression *, double> WeightMap;
  WeightMap m_WeightMap;

  double m_Constant;
};



/**
  Sample some external fixed function at a location. The function is specified
  by the traits object, which must support the Evaluate(x,y,z,order_x,order_y,order_z)
  operator.
  */
template <class TFunctionTraits>
class SampleFunction3 : public CachingExpression
{
public:
  typedef SampleFunction3<TFunctionTraits> Self;

  SampleFunction3(Problem *parent,
                  TFunctionTraits *traits,
                  Expression *x, Expression *y, Expression *z)
    : CachingExpression(parent), m_Traits(traits), m_X(x), m_Y(y), m_Z(z)
  {
    // Set to zero-th order derivatives by default
    m_OrderX = 0; m_OrderY = 0; m_OrderZ = 0;
    m_DfDx = NULL, m_DfDy = NULL, m_DfDz = NULL;
  }

  void SetOrder(int ox, int oy, int oz)
  {
    m_OrderX = ox; m_OrderY = oy; m_OrderZ = oz;
  }

  virtual double ComputeValue()
  {
    return m_Traits->Evaluate(m_X->Evaluate(), m_Y->Evaluate(), m_Z->Evaluate(),
                              m_OrderX, m_OrderY, m_OrderZ);
  }

  virtual void MakeDependencyList(Problem::Dependency &vars)
  {
    m_Problem->AppendDependentVariables(m_X, vars);
    m_Problem->AppendDependentVariables(m_Y, vars);
    m_Problem->AppendDependentVariables(m_Z, vars);
  }

  virtual Expression *MakePartialDerivative(Variable *variable)
  {
    // Get the partials of X, Y, Z
    Expression *dX = m_Problem->GetPartialDerivative(m_X, variable);
    Expression *dY = m_Problem->GetPartialDerivative(m_Y, variable);
    Expression *dZ = m_Problem->GetPartialDerivative(m_Z, variable);

    // Define the gradient if undefined - notice it is only created once
    if(!m_DfDx)
      {
      m_DfDx = new SampleFunction3(m_Problem, m_Traits, m_X, m_Y, m_Z);
      m_DfDx->SetOrder(m_OrderX+1, m_OrderY, m_OrderZ);
      m_DfDy = new SampleFunction3(m_Problem, m_Traits, m_X, m_Y, m_Z);
      m_DfDy->SetOrder(m_OrderX, m_OrderY+1, m_OrderZ);
      m_DfDz = new SampleFunction3(m_Problem, m_Traits, m_X, m_Y, m_Z);
      m_DfDz->SetOrder(m_OrderX, m_OrderY, m_OrderZ+1);
      }

    // Create the necessary expressions - chain rule
    WeightedSumGenerator wsg(m_Problem);
    if(dX) wsg.AddTerm(MakeProduct(m_Problem, m_DfDx, dX));
    if(dY) wsg.AddTerm(MakeProduct(m_Problem, m_DfDy, dY));
    if(dZ) wsg.AddTerm(MakeProduct(m_Problem, m_DfDz, dZ));
    return wsg.GenerateSum();
  }

  std::string GetName() { return "f(x,y,z)"; }

  virtual void GetTreeStats(int &depth, int &size)
  {
    int da, sa, db, sb, dc, sc;
    m_X->GetTreeStats(da, sa);
    m_Y->GetTreeStats(db, sb);
    m_Z->GetTreeStats(dc, sc);
    depth = 1 + std::max(dc, std::max(da, db));
    size = 1 + sa + sb + sc;
  }

  virtual int GetNumberOfOperands() const { return 3; }
  virtual Expression* GetOperand(int i) const
  {
    switch(i)
      {
      case 0 : return m_X;
      case 1 : return m_Y;
      case 2 : return m_Z;
      default: return NULL;
      }
  }

protected:
  Expression *m_X, *m_Y, *m_Z;
  Self *m_DfDx, *m_DfDy, *m_DfDz;
  int m_OrderX, m_OrderY, m_OrderZ;
  TFunctionTraits *m_Traits;
};


/**
  A sum of several elements
  */
class BigSum : public CachingExpression
{
public:
  BigSum(Problem *parent) : CachingExpression(parent) {}

  void AddSummand(Expression *ex)
  {
    m_A.push_back(ex);
  }

  virtual double ComputeValue()
  {
    double val = 0;
    for(Iterator it = m_A.begin(); it!=m_A.end(); ++it)
      val += (*it)->Evaluate();
    return val;
  }

  virtual void MakeDependencyList(Problem::Dependency &vars)
  {
    for(Iterator it = m_A.begin(); it!=m_A.end(); ++it)
      m_Problem->AppendDependentVariables(*it, vars);
  }

  virtual Expression *MakePartialDerivative(Variable *variable);

  std::string GetName();

  virtual void GetTreeStats(int &depth, int &size)
  {
    int di, si;
    depth = 0; size = 1;
    for(Iterator it = m_A.begin(); it!=m_A.end(); ++it)
      {
      Expression *ei = *it;
      ei->GetTreeStats(di, si);
      depth = std::max(depth, di);
      size += si;
      }
  }

  virtual int GetNumberOfOperands() const { return m_A.size(); }
  virtual Expression* GetOperand(int i) const
  {
    return (i < m_A.size()) ? m_A[i] : NULL;
  }

protected:
  typedef std::vector<Expression *>::iterator Iterator;
  std::vector<Expression *> m_A;
};


/** Typedefs for arrays of expressions */
typedef std::vector<Expression *> VarVec;
typedef std::vector<VarVec> VarVecArray;


/**
  A helper function to apply expressions to each element in a vector of
  expressions, useful for vector addition, subtraction, dot-product.
  */
template<class TExpression>
VarVec VectorApplyPairwise(Problem *p, const VarVec &a, const VarVec &b)
{
  VarVec out(a.size(), NULL);
  for(int i = 0; i < a.size(); i++)
    out[i] = new TExpression(p, a[i], b[i]);
  return out;
}

inline vnl_vector<double> VectorEvaluate(VarVec &vec)
{
  vnl_vector<double> v(vec.size());
  for(int i = 0; i < vec.size(); i++)
    v[i] = vec[i]->Evaluate();
  return v;
}

VarVec CrossProduct(Problem *p, const VarVec &a, const VarVec &b);
Expression *DotProduct(Problem *p, const VarVec &a, const VarVec &b);
Expression *TriangleTwiceAreaSqr(Problem *p, const VarVec &a, const VarVec &b, const VarVec &c);
Expression *DistanceSqr(Problem *p, const VarVec &a, const VarVec &b);
Expression *MagnitudeSqr(Problem *p, const VarVec &a);

/**
  A problem for IPOpt
  */
class ConstrainedNonLinearProblem : public Problem
{
public:
  static const double LBINF, UBINF;

  // The Jacobian of G is stored as a sparse array of Expression pointers
  typedef ImmutableSparseArray<Expression *> SparseExpressionMatrix;
  typedef ImmutableSparseMatrix<double> SparseRealMatrix;

  ConstrainedNonLinearProblem();
  virtual ~ConstrainedNonLinearProblem() {}

  // Add a variable to the NLP. Even though the variables are already known
  // to the problem when they are constructed, this method needs to be called
  // to actually include the variable in the optimization
  void AddVariable(Variable *var, double cmin = LBINF, double cmax = UBINF);

  // Create and add a variable in a single call, same as
  // Varirable v = new Variable(p, name); v.SetValue(x); p.AddVariable(v, ...)
  Variable *AddVariable(std::string name, double val = 0,
                        double cmin = LBINF, double cmax = UBINF);

  // Set the function to be minimized (f)
  void SetObjective(Expression *ex);

  // Add a constraint (g)
  void AddConstraint(Expression *ex, const std::string &category, double cmin, double cmax);

  // Add a constrained variable, i.e. a variable is forced to equal expression
  Variable* AddExpressionAsConstrainedVariable(Expression *exp, const std::string &category);

  // Setup the problem for optimization. This must be called after all the variables
  // constraints, and the objective have been specified
  void SetupProblem(bool useHessian, bool deriv_test = false);

  // As an option, we can provide a matrix for regularizing the gradient. The
  // gradient will be multiplied by this matrix before being passed to the
  // optimizer. This may have some benefit for convergence
  void SetGradientSmoothingKernel(SparseRealMatrix K);

  // Get the kernel for gradient smoothing
  SparseRealMatrix &GetGradientSmoothingKernel() { return m_GradientKernel; }

  // Get the number of variables
  unsigned int GetNumberOfVariables();

  // Get the number of constraints
  unsigned int GetNumberOfConstraints();

  // Get the variable
  Variable *GetVariable(int i) { return m_X[i]; }

  // Get the value of a variable
  double GetVariableValue(int i) { return m_X[i]->Evaluate(); }

  // Get the bounds of a variable
  void GetVariableBounds(int i, double &lb, double &ub);

  // Set the values of the variables
  void SetVariableValues(const double *x);

  // Set the values of the Lambdas
  void SetLambdaValues(const double *l);

  // Set the value of Sigma (see IPOpt)
  void SetSigma(double value);

  // Get the objective function
  Expression *GetObjective() { return m_F; }

  // Get the partial derivative of the objective with respect to the i-th
  // variable. This may return NULL if the objective does not depend on some
  // of the variables.
  Expression *GetObjectivePD(unsigned int i) { return m_GradF[i]; }

  // Get the i-th constraint
  Expression *GetConstraint(unsigned int j) { return m_G[j]; }

  // Get the bounds of a constraint
  void GetConstraintBounds(int i, double &lb, double &ub);

  // Get the category of a constraint
  const std::string &GetConstraintCategory(int i) { return m_ConstraintCategory[i]; }

  // Get the Jacobian sparse matrix
  SparseExpressionMatrix &GetConstraintsJacobian() { return m_DG; }

  // Get the Hessian of the Lagrangean
  SparseExpressionMatrix &GetHessianOfLagrangean() { return m_Hessian; }

protected:

  Expression *m_F;
  std::vector<Variable *> m_X, m_Lambda;
  std::vector<Expression *> m_G;
  std::vector<Expression *> m_GradF;

  std::vector<double> m_LowerBoundX, m_UpperBoundX;
  std::vector<double> m_LowerBoundG, m_UpperBoundG;
  std::vector<std::string> m_ConstraintCategory;

  // The matrix for the Jacobian
  SparseExpressionMatrix m_DG;

  // The matrix for the Hessian of the Lagrangean
  SparseExpressionMatrix m_Hessian;

  // For IPOpt, during Hessian computation, there is a scaling factor
  Variable *m_SigmaF;

  // Kernel for the gradient
  SparseRealMatrix m_GradientKernel;
};

}

#include "SparseMatrix.txx"
template class ImmutableSparseArray<gnlp::Expression *>;


#endif // GENTLENLP_H
