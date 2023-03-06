#ifndef QCQPPROBLEM_H
#define QCQPPROBLEM_H

#include <string>
#include <vector>
#include <map>
#include <vnl_vector.h>
#include <vnl_vector_fixed.h>

namespace qcqp {

static constexpr double LBINF = -1e100;
static constexpr double UBINF = +1e100;

template <unsigned int VDim>
struct Index
{
  int first;
  Index<VDim-1u> rest;

  template<class... T>
  Index(int in_first, T... in_rest) : first(in_first), rest(in_rest...) {}

  Index(int in_first = 0) : first(in_first) {}

  int product() const { return first * rest.product(); }

  Index<VDim> make_offset_table() const
    {
    return Index<VDim>({rest.product(), rest.make_offset_table()});
    }

  int flatten(const Index<VDim> &offset_table) const
    {
    return offset_table.first * first + rest.flatten(offset_table.rest);
    }
};

template<>
struct Index<1u>
{
  int first;
  Index(int in_first = 0) : first(in_first) {}
  int product() const { return first; }
  Index<1u> make_offset_table() const { return Index<1u>(1); }
  int flatten(const Index<1> &offset_table) const { return first; }
};

Index<5> test(1,2,3,4,5);
Index<5> test_ot = test.make_offset_table();
int z = test.product();
int y = test.flatten(test_ot);

// Partial specification, same as test_partial(1,2,3,0,0)
Index<5> test_partial(1,2,3);


class VariableBase
{
public:
  VariableBase(const std::string &name, double cmin, double cmax, int start_index)
    : name(name), cmin(cmin), cmax(cmax), start_index(start_index) {}

  friend class ElementRef;
  friend class QuadraticExpression;
  friend class LinearExpression;

protected:

  // Name of the variable
  std::string name;

  // Upper and lower bounds
  double cmin, cmax;

  // Index of the variable in complete variable array
  int start_index;

  // Data array
  vnl_vector<double> data;
};

/**
 * A reference to a position in the variable - used to access data and also
 * to set elements in constraint matrices/vectors
 */
class ElementRef
{
public:
  ElementRef(VariableBase &v, int offset) : v(v), offset(offset) {}

  operator double() const { return v.data[offset]; }
  ElementRef &operator = (double x) { v.data[offset] = x; return *this; }

  template <unsigned int N>
  vnl_vector_fixed<double, N> as_vector() const
    {
    return vnl_vector_fixed<double, N>(v.data.data_block() + offset);
    }

  template <unsigned int N>
  ElementRef &from_vector(vnl_vector_fixed<double, N> &x)
    {
    for(unsigned int i = 0; i < N; i++)
      v.data[offset + i] = x[i];
    return *this;
    }

  friend class QuadraticExpression;
  friend class LinearExpression;

protected:
  VariableBase &v;
  int offset;
};

/**
 * @brief Storage for an optimization variable, represented as a tensor of doubles
 */
template <unsigned int VDim>
struct Variable : public VariableBase
{
public:

  // Constructor
  Variable(const std::string &name, const Index<VDim> &size, double cmin, double cmax, int start_index)
    : VariableBase({name, cmin, cmax, start_index}), size(size), offset_table(size.make_offset_table())
    {
    data.set_size(size.product());
    }

  // Data access
  ElementRef operator()(const Index<VDim> &idx)
    {
    int offset = idx.flatten(offset_table);
    return ElementRef(*this, offset);
    }

  // Data access
  const ElementRef operator()(const Index<VDim> &idx) const
    {
    int offset = idx.flatten(offset_table);
    return ElementRef(*this, offset);
    }

protected:
  // Size of the variable
  Index<VDim> size;

  // Offset table
  Index<VDim> offset_table;
};

/**
 * Entry in the sparse matrices when constructing them, just a double
 * that takes default value 0
 */
struct Weight
{
  Weight() : value(0.0) {}
  Weight(double v) : value(v) {}
  double value;
  operator double() const { return value; }
  Weight &operator= (double v) { value = v; return *this; }
};

/**
 * A linear expression of the form b^t x + c. Can be used to build up
 * quadratic expressions, e.g., using dot product operation
 */
class LinearExpression
{
public:
  friend class Problem;
  friend class QuadraticExpression;

  double &b(ElementRef e1)
    {
    int pos1 = e1.v.start_index + e1.offset;
    return m_b[pos1].value;
    }

  double &c() { return m_c; }

protected:

  using VectorMap = std::map< int, Weight >;
  VectorMap m_b;
  double m_c = 0.0;
};

/**
 * A quadratic expression of the form x^t A x + b^t x + c. It represents either a
 * constraint or an objective function term. It can be built up by adding weights
 * corresponding to variables in a problem
 */
class QuadraticExpression
{
public:
  friend class Problem;

  double &A(ElementRef e1, ElementRef e2)
    {
    int pos1 = e1.v.start_index + e1.offset;
    int pos2 = e2.v.start_index + e2.offset;
    return m_A[{pos1, pos2}].value;
    }

  double &b(ElementRef e1)
    {
    int pos1 = e1.v.start_index + e1.offset;
    return m_b[pos1].value;
    }

  double &c() { return m_c; }

  void add_dotted(const LinearExpression &l1, const LinearExpression &l2)
    {
    // We need to compute the outer product of the two linear expressions
    for(auto i1 : l1.m_b)
      {
      for(auto i2 : l2.m_b)
        {
        m_A[{i1.first, i2.first}].value += i1.second * i2.second;
        }
      m_b[i1.first].value += i1.second * l2.m_c;
      }
    for(auto i2 : l2.m_b)
      {
      m_b[i2.first].value += i2.second * l1.m_c;
      }
    m_c += l1.m_c * l2.m_c;
    }

protected:

  using MatrixMap = std::map< std::tuple<int, int>, Weight>;
  using VectorMap = std::map< int, Weight >;
  MatrixMap m_A;
  VectorMap m_b;
  double m_c = 0.0;
};

/**
 * Storage for constraints. Each constraint is associated with a sparse matrix C_k, a
 * sparse vector d_k, a constant term f_k, as well as lower and upper bounds.
 */
class Constraint : public QuadraticExpression
{
public:
  Constraint(const std::string &name, double lower, double upper)
    : name(name), lower_bound(lower), upper_bound(upper) {}

protected:

  // Name of constraint
  std::string name;

  // Bounds
  double lower_bound, upper_bound;
};


/**
 * Storage for objective terms. Each term is a quadratic expression with an
 * additional weight, useful for configuration
 */
class Loss : public QuadraticExpression
{
public:
  Loss(const std::string &name, double weight)
    : name(name), weight(weight) {}

protected:

  // Name of constraint
  std::string name;

  // Bounds
  double weight;
};



/**
 * @brief Helper for building a quadratically constrained quadratic programming problem.
 *
 */
class Problem
{
public:


  /**
   * Add a variable to the optimization. The variable is really just an index into the
   * array of optimization variables that can then be referenced by name conveniently.
   *
   * The method supports adding an array of indexed variables by setting the size value
   * to the array size.
   */
  template <unsigned int VDim>
  Variable<VDim> &AddVariable(const char *id, Index<VDim> size, double cmin = LBINF, double cmax = UBINF)
    {
    // Create the variable
    Variable<VDim> *v = new Variable<VDim>({ std::string(id), size, cmin, cmax, (int) m_Size });
    m_Size += size.product();

    // Store the pointer
    m_Variables.push_back(v);

    // Return a const reference to the variable
    return *v;
    }

  /**
   * Add a constraint or set of similar constraints. Each constraint is associated with
   * a quadratic expression alpha <= X^t C_k X + d^t X <= beta. Once the constrainet
   */
  Constraint &AddConstraint(const char *id, double cmin, double cmax)
    {
    Constraint *c = new Constraint(id, cmin, cmax);
    m_Constraints.push_back(c);
    return *c;
    }

  /**
   * Add a loss term. Each loss is also quadratic and has a weight
   */
  Loss &AddLoss(const char *id, double weight)
    {
    Loss *l = new Loss(id, weight);
    m_Losses.push_back(l);
    return *l;
    }

protected:

  // Known variables
  std::vector<VariableBase *> m_Variables;

  // Known constraints
  std::vector<Constraint *> m_Constraints;

  // Known losses
  std::vector<Loss *> m_Losses;

  // Total number of variables so far
  unsigned int m_Size = 0;

};

}


#endif // QCQPPROBLEM_H
