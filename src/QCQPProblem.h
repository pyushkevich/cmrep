#ifndef QCQPPROBLEM_H
#define QCQPPROBLEM_H

#include <string>
#include <vector>
#include <map>
#include <vnl_vector.h>
#include <vnl_vector_fixed.h>
#include <tuple>
#include <unordered_set>
#include "SparseMatrix.h"

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
  int flatten(const Index<1> &) const { return first; }
};

class VariableRefBase
{
public:
  VariableRefBase(const std::string &name, double cmin, double cmax,
                  int start_index, int flat_size)
    : name(name), cmin(cmin), cmax(cmax),
      start_index(start_index), flat_size(flat_size) {}

  VariableRefBase()
    : cmin(LBINF), cmax(UBINF), start_index(-1), flat_size(0) {}

  const char *GetName() const { return name.c_str(); }
  int GetFlatSize() const { return flat_size; }
  int GetStartIndex() const { return start_index; }

  friend class ElementConstRef;
  friend class ElementRef;
  friend class QuadraticExpression;
  friend class LinearExpression;
  friend class Problem;

protected:

  // Name of the variable
  std::string name;

  // Upper and lower bounds
  double cmin, cmax;

  // Index of the variable in complete variable array
  int start_index, flat_size;

  // Pointer to the data for this variable, managed by the problem
  double *data = nullptr;
};

/**
 * A reference to a position in the variable - used to access data and also
 * to set elements in constraint matrices/vectors
 */
class ElementConstRef
{
public:
  ElementConstRef(VariableRefBase &v, int offset) : v(v), offset(offset) {}

  operator double() const { return v.data[offset]; }

  template <unsigned int N>
  vnl_vector_fixed<double, N> as_vector() const
    {
    return vnl_vector_fixed<double, N>(v.data + offset);
    }

  friend class QuadraticExpression;
  friend class LinearExpression;

protected:
  VariableRefBase &v;
  int offset;
};

class ElementRef : public ElementConstRef
{
public:
  ElementRef(VariableRefBase &v, int offset) : ElementConstRef(v, offset) {}

  ElementRef &operator = (double x) { v.data[offset] = x; return *this; }

  template <unsigned int N>
  ElementRef &from_vector(const vnl_vector_fixed<double, N> &x)
    {
    for(unsigned int i = 0; i < N; i++)
      v.data[offset + i] = x[i];
    return *this;
    }

  friend class QuadraticExpression;
  friend class LinearExpression;
};

/**
 * @brief Storage for an optimization variable, represented as a tensor of doubles
 */
template <unsigned int VDim>
struct VariableRef : public VariableRefBase
{
public:

  // Constructor
  VariableRef(const std::string &name, const Index<VDim> &size, double cmin, double cmax, int start_index)
    : VariableRefBase({name, cmin, cmax, start_index, size.product()}), size(size), offset_table(size.make_offset_table())
    {
    }

  VariableRef() {}

  template<class... T>
  ElementRef operator()(int in_first, T... in_rest)
    {
    int offset = Index<VDim>(in_first, in_rest...).flatten(offset_table);
    return ElementRef(*this, offset);
    }

  template<class... T>
  ElementConstRef operator()(int in_first, T... in_rest) const
    {
    int offset = Index<VDim>(in_first, in_rest...).flatten(offset_table);
    return ElementConstRef(*this, offset);
    }

protected:
  // Size of the variable
  Index<VDim> size;

  // Offset table
  Index<VDim> offset_table;
};

/**
 * A wrapper around an atomic type that always initializes to zero
 */
template <class T>
struct zero_init
{
  zero_init() : value(0) {}
  zero_init(T v) : value(v) {}
  operator T() const { return value; }
  zero_init<T> &operator= (T v) { value = v; return *this; }

  T value;
};

/**
 * Entry in the sparse matrices when constructing them, just a double
 * that takes default value 0
 */
using Weight = zero_init<double>;

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
 *
 * The matrix A is triangular; regardless of the order in which variables
 * are specified, they will be ordered to insert into the upper triangle (row <= column)
 */
class QuadraticExpression
{
public:
  friend class Problem;

  double &A(ElementRef e1, ElementRef e2)
    {
    int pos1 = e1.v.start_index + e1.offset;
    int pos2 = e2.v.start_index + e2.offset;
    return A(pos1, pos2);
    }

  double &A(int pos1, int pos2)
    {
    m_Dependencies.insert(pos1);
    m_Dependencies.insert(pos2);
    if(pos1 <= pos2)
      return m_Aij[{pos1, pos2}];
    else
      return m_Aij[{pos2, pos1}];
    }

  double &b(ElementRef e1)
    {
    int pos1 = e1.v.start_index + e1.offset;
    return b(pos1);
    }

  double &b(int pos1)
    {
    m_Dependencies.insert(pos1);
    return m_bi[pos1];
    }

  double &c() { return m_c; }

  void add_dotted(const LinearExpression &l1, const LinearExpression &l2)
    {
    // We need to compute the outer product of the two linear expressions
    for(auto i1 : l1.m_b)
      {
      for(auto i2 : l2.m_b)
        {
        int pos1 = std::min(i1.first, i2.first);
        int pos2 = std::max(i1.first, i2.first);
        A(pos1, pos2) += i1.second * i2.second;
        }
      b(i1.first) += i1.second * l2.m_c;
      }
    for(auto i2 : l2.m_b)
      {
      b(i2.first) += i2.second * l1.m_c;
      }
    m_c += l1.m_c * l2.m_c;
    }

  void add(const QuadraticExpression &qe)
    {
    for(auto a : qe.m_Aij)
      m_Aij[a.first] += a.second;
    for(auto b : qe.m_bi)
      m_bi[b.first] += b.second;
    m_c += qe.m_c;
    }

  template <typename T>
  double Evaluate(const T* x)
  {
    // Evaluating the function means computing the quadratic form xAx+bx+c
    // for each of the loss terms. In the future we might want to create more
    // efficient data structures for these constraints
    double y = m_c;
    for(auto it_a : m_Aij)
      {
      int i = std::get<0>(it_a.first);
      int j = std::get<1>(it_a.first);
      y += x[i] * x[j] * it_a.second;
      }
    for(auto it_b : m_bi)
      {
      int i = it_b.first;
      y += x[i] * it_b.second;
      }

    return y;
  }

  template <typename T>
  void ComputeGradient(const T* x, T* grad_x, double weight)
  {
    for(auto it_a : m_Aij)
      {
      int i = std::get<0>(it_a.first);
      int j = std::get<1>(it_a.first);
      grad_x[i] += x[j] * it_a.second * weight;
      grad_x[j] += x[i] * it_a.second * weight;
      }
    for(auto it_b : m_bi)
      {
      int i = it_b.first;
      grad_x[i] += it_b.second * weight;
      }
  }

protected:

  // Sparse storage for the components a_ij
  std::map< std::tuple<int,int>, double> m_Aij;
  std::map< int, double> m_bi;
  double m_c = 0.0;

  // Set of all variables that we depend on
  std::unordered_set<int> m_Dependencies;

  /*

  // For converting to sparse matrices, it is more efficient to
  // index by row and column separately, rather than by a row/col
  // index
  using MatrixMap = std::map< int, LinearExpression >;
  MatrixMap m_Ab;

  // The only thing stored separately is the offset c
  double m_c = 0.0;

  */
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

  friend class Problem;

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

  friend class Problem;

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
  // The Jacobian and the Hessian of Lagrangean are represented as sparse
  // 2D arrays, where each entry is used to compute the expression w^t x + z
  // or w^t lambda + z.
  struct TensorRow
  {
    std::vector< std::tuple<size_t, double> > w;
    double z;
  };

  using SparseTensor = ImmutableSparseArray<TensorRow>;

  /**
   * Add a variable tensor to the optimization. The variable is really just an index into the
   * array of optimization variables that can then be referenced by name conveniently.
   *
   * The method supports adding an array of indexed variables by setting the size value
   * to the array size.
   */
  template <unsigned int VDim>
  VariableRef<VDim> AddVariableTensor(const char *id, Index<VDim> size, double cmin = LBINF, double cmax = UBINF)
    {
    // Create the variable reference
    VariableRef<VDim> *v = new VariableRef<VDim>({ std::string(id), size, cmin, cmax, (int) m_Size });

    // Create the data for this variable and place in the variable ref
    m_VariableData[v->name].set_size(v->flat_size);
    v->data = m_VariableData[v->name].data_block();

    // Update the total problem size
    m_Size += v->flat_size;

    // Store the pointer
    m_Variables.push_back(v);

    // Return a const reference to the variable
    return *v;
    }

  /** Get the number of added variable tensors */
  int GetNumberOfVariableTensors() const { return m_Variables.size(); }

  /** Get the poinster to the i-th variable tensor */
  const VariableRefBase *GetVariableTensor(int i) const { return m_Variables[i]; }

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

  /** Get number of registered variables */
  unsigned int GetNumberOfVariables() const { return m_Size; }

  /** Get number of registered constraints */
  unsigned int GetNumberOfConstraints() const { return m_Constraints.size(); }

  /** Get number of registered constraints */
  unsigned int GetNumberOfLosses() const { return m_Losses.size(); }

  /** Get the sparse array used to compute the Jacobian */
  SparseTensor &GetConstraintsJacobian() { return m_Jacobian; }

  /** Get the sparse array used to compute the Jacobian */
  SparseTensor &GetHessianOfLagrangean() { return m_Hessian; }

  /** Get the bounds on variables */
  void GetVariableBounds(int i_var, double &lb, double &ub)
  {
    lb = m_VariableLB[i_var];
    ub = m_VariableUB[i_var];
  }

  /** Get the bounds on constraints */
  void GetConstraintBounds(int i_con, double &lb, double &ub)
  {
    lb = m_Constraints[i_con]->lower_bound;
    ub = m_Constraints[i_con]->upper_bound;
  }

  /** Get the initial value of a variable */
  double GetVariableValue(int i)
  {
    return m_VariableValue[i];
  }

  /** Update the value of a variable */
  void SetVariableValue(int i, double x)
  {
    m_VariableValue[i] = x;
  }

  /** Get all variables as a flat array */
  vnl_vector<double> GetVariableValues()
    {
    vnl_vector<double> x(m_Size);
    for(auto v : m_Variables)
      {
      for(int i = 0; i < v->flat_size; i++)
        x[v->start_index+i] = v->data[i];
      }
    return x;
    }

  /** Compute the objective function for given x */
  template <typename T>
  double EvaluateLoss(const T* x)
  {
    // Evaluating the function means computing the quadratic form xAx+bx+c
    // for each of the loss terms. In the future we might want to create more
    // efficient data structures for these constraints
    double total_loss = 0;
    for(auto *l : m_Losses)
      {
      double value = l->Evaluate(x);
      total_loss += value * l->weight;
      }

    return total_loss;
  }

  using LossReport = std::map< std::string, std::tuple<double, double> >;

  /** Print the loss report */
  template <typename T>
  LossReport GetLossReport(const T* x)
  {
    LossReport report;
    for(auto *l : m_Losses)
      report[l->name] = { l->weight, l->Evaluate(x) };
    return report;
  }


  /** Compute the gradient of the objective function for given x */
  template <typename T>
  void EvaluateLossGradient(const T* x, T *grad_x)
  {
    // Initialize the gradient to zero
    for(unsigned int i = 0; i < m_Size; i++)
      grad_x[i] = 0.0;

    // Evaluating the function means computing the quadratic form xAx+bx+c
    // for each of the loss terms. In the future we might want to create more
    // efficient data structures for these constraints
    for(auto *l : m_Losses)
      l->ComputeGradient(x, grad_x, l->weight);
  }

  /** Compute the individual constraints */
  template <typename T>
  double EvaluateConstraint(const T* x, int k)
  {
    // Evaluating the function means computing the quadratic form xAx+bx+c
    // for each of the loss terms. In the future we might want to create more
    // efficient data structures for these constraints
    return m_Constraints[k]->Evaluate(x);
  }

  /** Get the name of the constraint */
  std::string GetConstraintName(int k) const
  {
    return m_Constraints[k]->name;
  }


  /**
   * Setup the problem
   */
  void SetupProblem()
  {
    // Compute the constraints Jacobian. Each constraint is in the form
    // lb <= g(x) <= ub, where g(x) = sum_ij a_ij x_i x_i + sum_i b_i x_i + c
    //
    // so its derivative with respect to x_k is
    //   sum_i a_ik x_i + sum_i a_ki x_i + b_k
    //
    // The Jacobian matrix is a matrix where every row corresponds
    // to a constraint, and every column to a variable, i.e., columns are
    // the z-values.

    // Allocate the row array
    size_t ncon = m_Constraints.size();
    size_t *jac_row_arr = new size_t[ncon + 1];
    jac_row_arr[0] = 0;

    // The Jacobian for constraint g(x) wrt x_k is non-zero if k appears in the
    // non-zero entries for g as either row or column

    // For each constraint, we need to compute the gradient entries. Right now we
    // do this using a map, which might be slow
    struct TensorRowSource { double z; std::map<int, double> w; };
    std::vector< std::map<int, TensorRowSource> > con_jac_map(ncon);
    for(size_t i = 0; i < ncon; i++)
      {
      auto &c = *m_Constraints[i];
      auto &cm = con_jac_map[i];

      for(auto it_a : c.m_Aij)
        {
        int i = std::get<0>(it_a.first);
        int j = std::get<1>(it_a.first);
        cm[i].w[j] += it_a.second;
        cm[j].w[i] += it_a.second;
        }
      for(auto it_b : c.m_bi)
        {
        int i = it_b.first;
        cm[i].z += it_b.second;
        }
      }

    // For each row, count the number of nonzero entries in the Jacobian
    for(size_t i = 0; i < ncon; i++)
      jac_row_arr[i+1] = jac_row_arr[i] + con_jac_map[i].size();

    // Allocate storage for the entries
    size_t *jac_col_arr = new size_t[jac_row_arr[ncon]];
    TensorRow *jac_val = new TensorRow[jac_row_arr[ncon]];

    // Second pass - allocate individual entries
    for(size_t i = 0; i < ncon; i++)
      {
      // Reference to the constraint
      auto &c = *m_Constraints[i];

      // Record the indices of the non_zeros in the Jacobian matrix
      int j = jac_row_arr[i];

      // Iterate over the rows of the constraint matrix m_Ab
      for(auto &it: con_jac_map[i])
        {
        jac_col_arr[j] = it.first;
        jac_val[j].w.reserve(it.second.w.size());
        for(auto it2: it.second.w)
          jac_val[j].w.push_back({it2.first, it2.second});
        jac_val[j].z = it.second.z;
        j++;
        }
      }

    // Create the Jacobian matrix
    m_Jacobian.SetArrays(ncon, m_Size, jac_row_arr, jac_col_arr, jac_val);

    // Now compute the Hessian of the Lagrangian. This has the dimensions
    // n x n where n is the number of variables, and when computing we need
    // to iterate over the Lagrange multipliers. To form this matrix, we need
    // to 'splat' each of the constraints' A matrices

    // What we essentially have to do is to transpose the way in which the
    // constraint information is stored so that it is row/col/constraint order
    // as opposed to the constraint/row/col order. One way to do this is using
    // maps of maps - might be slow but let's try
    struct SrcEntry
    {
      double h_loss = 0.;             // Hessian of the losses
      std::map<int, double> h_con;    // Constraint Hessians
    };

    std::vector< std::map<int, SrcEntry > > hol_src;

    // No sparsity in the rows
    hol_src.resize(m_Size);

    // Iterate over the constraints
    for(unsigned int k = 0; k < ncon; k++)
      {
      // Within the constraint, iterate over NNZ rows
      for(auto it : m_Constraints[k]->m_Aij)
        {
        int i = std::get<0>(it.first);
        int j = std::get<1>(it.first);
        if(i <= j)
          hol_src[i][j].h_con[k] += it.second;
        if(j <= i)
          hol_src[j][i].h_con[k] += it.second;
        }
      }

    // We also need to add the losses to the hessian, we can add them
    // using index -1
    for(auto *l : m_Losses)
      {
      // Within the loss, iterate over NNZ rows
      for(auto it: l->m_Aij)
        {
        int i = std::get<0>(it.first);
        int j = std::get<1>(it.first);
        if(i <= j)
          hol_src[i][j].h_loss += it.second * l->weight;
        if(j <= i)
          hol_src[j][i].h_loss += it.second * l->weight;
        }
      }

    // Set the row index
    size_t *hol_row_arr = new size_t[m_Size + 1];
    hol_row_arr[0] = 0;
    for(size_t i = 0; i < m_Size; i++)
      hol_row_arr[i+1] = hol_row_arr[i] + hol_src[i].size();

    // We can now allocate the HoL column index and data
    size_t *hol_col_arr = new size_t[hol_row_arr[m_Size]];
    TensorRow *hol_val = new TensorRow[hol_row_arr[m_Size]];

    // Now set the values of the column arrays
    for(size_t i = 0; i < m_Size; i++)
      {
      size_t j = hol_row_arr[i];
      for(auto &it_row : hol_src[i])
        {
        hol_col_arr[j] = it_row.first;
        hol_val[j].z = it_row.second.h_loss;
        hol_val[j].w.reserve(it_row.second.h_con.size());
        for(auto &it_w : it_row.second.h_con)
          hol_val[j].w.push_back({it_w.first, it_w.second});
        j++;
        }
      }

    // Create the sparse matrix
    m_Hessian.SetArrays(m_Size, m_Size, hol_row_arr, hol_col_arr, hol_val);

    // Create the bounds arrays
    m_VariableLB.set_size(m_Size);
    m_VariableUB.set_size(m_Size);
    m_VariableValue.set_size(m_Size);
    for(auto *v : m_Variables)
      {
      for(int i = 0; i < v->flat_size; i++)
        {
        m_VariableLB[i + v->start_index] = v->cmin;
        m_VariableUB[i + v->start_index] = v->cmax;
        m_VariableValue[i + v->start_index] = m_VariableData[v->name][i];
        }
      }

  }

  template <class TNumber>
  void SetWarmStartLambda(int n, const TNumber *values)
  {
    m_WarmStartLambda.set_size(n);
    for(int i = 0; i < n; i++)
      m_WarmStartLambda[i] = values[i];
  }

  template <class TNumber>
  void SetWarmStartZL(int n, const TNumber *values)
  {
    m_WarmStartZL.set_size(n);
    for(int i = 0; i < n; i++)
      m_WarmStartZL[i] = values[i];
  }

  template <class TNumber>
  void SetWarmStartZU(int n, const TNumber *values)
  {
    m_WarmStartZU.set_size(n);
    for(int i = 0; i < n; i++)
      m_WarmStartZU[i] = values[i];
  }

  const vnl_vector<double> &GetWarmStartLambda() const { return m_WarmStartLambda; }
  const vnl_vector<double> &GetWarmStartZL() const { return m_WarmStartZL; }
  const vnl_vector<double> &GetWarmStartZU() const { return m_WarmStartZU; }

  /** Maps variable values from flat array back to the per-variable arrays */
  void Finalize()
  {
    for(auto *v : m_Variables)
      {
      for(int i = 0; i < v->flat_size; i++)
        {
        m_VariableData[v->name][i] = m_VariableValue[i + v->start_index];
        }
      }
  }

protected:

  // Known variables
  std::vector<VariableRefBase *> m_Variables;

  // Known constraints
  std::vector<Constraint *> m_Constraints;

  // Known losses
  std::vector<Loss *> m_Losses;

  // Variables, referenced by name
  std::map<std::string, vnl_vector<double> > m_VariableData;

  // Total number of variables so far
  unsigned int m_Size = 0;

  // Jacobian and Hessian computed for this problem
  SparseTensor m_Jacobian, m_Hessian;

  // Variable bounds - available after SetupProblem
  vnl_vector<double> m_VariableLB, m_VariableUB, m_VariableValue;

  // Storage for warm start data - Lagrange multipliers
  vnl_vector<double> m_WarmStartLambda, m_WarmStartZL, m_WarmStartZU;


};

}


#endif // QCQPPROBLEM_H
