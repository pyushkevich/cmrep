#include "GentleNLP.h"
#include "PointSetHamiltonianSystem.h"

/**
 * This expression represents the function
 *
 * f(p0,y) = |q1(p0) - y|^2 + H(p0) * lambda
 *
 * In other words, the SSD from a geodesic shooting result to a set of
 * target points. The expression is evaluated using geodesic shooting
 */
template <class TFloat, unsigned int VDim>
class GeodesicShootingToTargetExpression : public Expression
{
public:
  typedef PointSetHamiltonianSystem<TFLoat,VDim> HamiltonianSystem;
  typedef std::vector< std::vector< Constant * > > ConstVecArray;
  typedef std::vector< std::vector< Variable * > > VarVecArray;

  /**
   * @param parent
   * @param q0: should be constants (initial template does not move)!
   * @param p0: a set of variables
   * @param y: another set of variables
   * @param sigma
   * @param n_steps
   * @param lambda
   */
  GeodesicShootingToTargetExpression(
      Problem *parent,
      const ConstVecArray &q0, const VarVecArray &p0, const VarVecArray &y,
      double sigma, int n_steps, double lambda)
    : Expression(parent)
  {
    typename HamiltonianSystem::Matrix q0_mtx(q0.size(), VDim);
    for(int j = 0; j < q0.size; j++)
      for(int a = 0; a < VDim; a++)
        q0_mtx(j,a) = q0[j][a]->Evaluate();

    m_Shooting = new HamiltonianSystem(q0_mtx, sigma, n_steps);
    m_Dirty = m_GradientDirty = true;
    m_P0 = p0;
    m_Y = y;
    m_Lambda = lambda;
  }

  ~GeodesicShootingExpression() { delete m_System; }

  virtual double Evaluate()
  {
    if(m_Dirty)
      {
      typename HamiltonianSystem::Matrix p0_mtx(m_P0.size(), VDim);
      typename HamiltonianSystem::Matrix q1, p1;
      for(int j = 0; j < m_P0.size; j++)
        for(int a = 0; a < VDim; a++)
          p0_mtx(j,a) = m_P0[j][a]->Evaluate();

      double H = m_Shooting->FlowHamiltonian(p0_mtx,q1,p1);

      double d_targ = 0.0;
      for(int j = 0; j < m_P0.size; j++)
        {
        for(int a = 0; a < VDim; a++)
          {
          double dja = q1(j,a) - m_Y[j][a]->Evaluate();
          d_targ += dja * dja;
          }
        }

      m_Dirty = false;
      m_CachedEnergy = m_Lambda * H + d_targ;
      }

    return m_CachedEnergy;
  }

  virtual void MakeDirty()
  {
    m_Dirty = true;
    m_GradientDirty = true;
  }

  virtual void MakeDependencyList(Problem::Dependency &vars)
  {
    for(int j = 0; j < m_P0.size; j++)
      for(int a = 0; a < VDim; a++)
        {
        m_Problem->AppendDependentVariables(m_P0[j][a], vars);
        m_Problem->AppendDependentVariables(m_Y[j][a], vars);
        }
  }

  virtual Expression *MakePartialDerivative(Variable *variable)
  {
    // We need to find which index corresponds to this variable. Unlike
    // other expressions in Gentle NLP, we only allow input expressions to
    // be variables. This is because having momentum depend on other vars
    // would require a full jacobian computation instead of simple backprop
    for(int j = 0; j < m_P0.size; j++)
      {
      for(int a = 0; a < VDim; a++)
        {
        if(m_P0[j][a] == variable)
          {
          return new MomentumDerivativeExpression(this, j, a);
          }
        else if(m_Y[j][a] == variable)
          {
          return new TargetDerivativeExpression(this, j, a);
          }
        }
      }

    // Some other variable was passed in - this is not right!
    return NULL;
  }

  virtual std::string GetName()
  {
    return std::string("GS_to_trg");
  }




protected:
  HamiltonianSystem *m_Shooting;
  VarVecArray m_P0, m_Y;
  bool m_Dirty, m_GradientDirty;
  double m_Lambda;

  double m_CachedEnergy;
};

/**
 * This expression is a wrapper that returns the value of the geodesic
 * shooting expression q_1[p_0]
 */
class GeodesicShootingExpression : public Expression
{
public:
  GeodesicShootingExpression(
      Problem *parent, GeodesicShootingWrapper *gs, int index)
    : Expression(parent), m_Shooting(gs), m_Index(index)
  {
  }

  virtual double Evaluate()
  {
    // Get Q1 from the shooting
    m_Shooting->DoShooting();
    return m_Shooting()->GetQ1(m_Index);
  }

protected:
  GeodesicShootingWrapper *m_Shooting;
  int m_Index;
};

class GeodesicShootingWrapper
{
public:
  GeodesicShootingWrapper()
  {
    m_HamiltonianSystem->
  }


protected:
  // The hamiltonian system used to perform the shooting
  PointSetHamiltonianSystem *m_HamiltonianSystem;

  // Whether the shooting itself is dirty



}
