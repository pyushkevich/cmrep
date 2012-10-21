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

  Problem() {}
  virtual ~Problem();

  /** Get a partial derivative for an expression with repsect to a
    variable. The problem keeps track of these expressions so that
    repeat calls to GetPartialDerivative only create child expressions
    when the derivative has not previously been computed */
  Expression *GetPartialDerivative(Expression *, Variable *);

  void AddChildExpression(Expression *);

protected:

  // List of known expressions
  typedef std::set<Expression *> ExpressionSet;
  ExpressionSet m_ChildExpressions;

  // List of partial derivatives of known expressions
  typedef std::pair<Expression *, Variable *> PartialDerivative;
  typedef std::map<PartialDerivative, Expression *> PartialDerivativeMap;
  PartialDerivativeMap m_Partials;
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

  /** Evaluate the expression */
  virtual double Evaluate() = 0;

  /** Determine if the expression depends on a variable */
  virtual bool DependsOn(Variable *variable) = 0;

  /** Populate a list of dependent variables */
  virtual void GetDependentVariables(std::set<Variable *> &vars) = 0;

  /** Create a partial derivative of the expression with respect to variable.
      This should only be called once per expression/variable, since each call
      will make a new() expression. */
  virtual Expression *MakePartialDerivative(Variable *variable) = 0;

  /** Debugging purposes : print the expression using LateX */
  virtual std::string GetName() = 0;

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
    : Expression(parent), m_Name(name), m_Index(-1), m_Value(0) {}

  void SetValue(double value) { m_Value = value; }
  std::string GetName() { return m_Name; }

  void SetIndex(int value) { m_Index = value; }
  int GetIndex() { return m_Index; }

  virtual double Evaluate() { return m_Value; }
  virtual bool DependsOn(Variable *variable) { return this == variable; }
  virtual Expression *MakePartialDerivative(Variable *variable);

  virtual void GetDependentVariables(std::set<Variable *> &vars)
  {
    vars.insert(this);
  }

protected:
  std::string m_Name;
  int m_Index;
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
  virtual bool DependsOn(Variable *variable) { return false; }
  virtual Expression *MakePartialDerivative(Variable *variable);

  virtual void GetDependentVariables(std::set<Variable *> &vars) {}

protected:
  double m_Value;
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


template<class TOperatorTraits>
class UnaryExpression : public Expression
{
public:
  UnaryExpression(Problem *parent, Expression *a)
    : Expression(parent), m_A(a) {}

  virtual double Evaluate()
  {
    return TOperatorTraits::Operate(m_A->Evaluate());
  }

  virtual bool DependsOn(Variable *variable)
  {
    return m_A->DependsOn(variable);
  }

  virtual void GetDependentVariables(std::set<Variable *> &vars)
  {
    m_A->GetDependentVariables(vars);
  }

  virtual Expression *MakePartialDerivative(Variable *variable)
  {
    Expression *dA = m_Problem->GetPartialDerivative(m_A, variable);
    return dA ? TOperatorTraits::Differentiate(m_Problem, this, m_A, dA) : NULL;
  }

  std::string GetName() { return TOperatorTraits::GetName(m_A); }


protected:
  Expression *m_A;
};

typedef UnaryExpression<NegateOperatorTraits> Negation;
typedef UnaryExpression<SquareOperatorTraits> Square;



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
class BinaryExpression : public Expression
{
public:
  BinaryExpression(Problem *parent, Expression *a, Expression *b)
    : Expression(parent), m_A(a), m_B(b) {}

  virtual double Evaluate()
  {
    return TOperatorTraits::Operate(m_A->Evaluate(), m_B->Evaluate());
  }

  virtual bool DependsOn(Variable *variable)
  {
    return m_A->DependsOn(variable) || m_B->DependsOn(variable);
  }

  virtual void GetDependentVariables(std::set<Variable *> &vars)
  {
    m_A->GetDependentVariables(vars);
    m_B->GetDependentVariables(vars);
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

protected:
  Expression *m_A, *m_B;

};

typedef BinaryExpression<PlusOperatorTraits> BinarySum;
typedef BinaryExpression<MinusOperatorTraits> BinaryDifference;
typedef BinaryExpression<ProductOperatorTraits> BinaryProduct;
typedef BinaryExpression<RatioOperatorTraits> BinaryFraction;


/**
  A sum of several elements
  */
class BigSum : public Expression
{
public:
  BigSum(Problem *parent) : Expression(parent) {}

  void AddSummand(Expression *ex)
  {
    m_A.push_back(ex);
  }

  virtual double Evaluate()
  {
    double val = 0;
    for(Iterator it = m_A.begin(); it!=m_A.end(); ++it)
      val += (*it)->Evaluate();
    return val;
  }

  virtual bool DependsOn(Variable *variable)
  {
    for(Iterator it = m_A.begin(); it!=m_A.end(); ++it)
      if((*it)->DependsOn(variable))
        return true;
    return false;
  }

  virtual void GetDependentVariables(std::set<Variable *> &vars)
  {
    for(Iterator it = m_A.begin(); it!=m_A.end(); ++it)
      (*it)->GetDependentVariables(vars);
  }

  virtual Expression *MakePartialDerivative(Variable *variable);

  std::string GetName();

protected:
  typedef std::list<Expression *>::iterator Iterator;
  std::list<Expression *> m_A;

};


/**
  A problem for IPOpt
  */
class ConstrainedNonLinearProblem : public Problem
{
public:
  static const double LBINF, UBINF;

  // The Jacobian of G is stored as a sparse array of Expression pointers
  typedef ImmutableSparseArray<Expression *> SparseExpressionMatrix;

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
  void AddConstraint(Expression *ex, double cmin, double cmax);

  // Setup the problem for optimization. This must be called after all the variables
  // constraints, and the objective have been specified
  void SetupProblem();

  // Get the number of variables
  unsigned int GetNumberOfVariables();

  // Get the number of constraints
  unsigned int GetNumberOfConstraints();

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


  // The matrix for the Jacobian
  SparseExpressionMatrix m_DG;

  // The matrix for the Hessian of the Lagrangean
  SparseExpressionMatrix m_Hessian;

  // For IPOpt, during Hessian computation, there is a scaling factor
  Variable *m_SigmaF;
};

}

#include "SparseMatrix.txx"
template class ImmutableSparseArray<gnlp::Expression *>;


#endif // GENTLENLP_H
