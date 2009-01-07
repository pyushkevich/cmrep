#ifndef __AffineMapping_h_
#define __AffineMapping_h_

#include "GenericMedialModel.h"
#include <algorithm>

class PrincipalComponents;
class SubdivisionMedialModel;

/*
 * Coefficient Transformations
 *
 * Coefficient transformations are mappings of the form C' = T(C, P).  Here C'
 * and C are coefficient vectors in R^m and P is a parameter vector in R^n. An
 * example of the coefficient mapping is the affine transform. C are
 * untransformed coefficients of a medial model, P are the parameters of the
 * affine transform and C' are the coefficients corresponding to the transformed
 * model. 
 *
 * When working with coefficient transformations, we are interested in applying
 * them (i.e., computing C') and also in applying their Jacobians to arbitrary
 * variations (i.e., v = J_T(C) * v_C or v = J_T(P) * v_P). For example, to
 * compute the derivative of f(C') with respect to variation v_P, we need to
 * take the variational derivative of f with respect to variation J_T(P) * v_P.
 *
 * Finally, these coefficient mappings are designed to be used in
 * optimization. The starting point for this optimization should be P = 0.
 * During optimization, P_opt will be found such that T(C, P_opt) is
 * maximized.
 */

/**
 * This class represents a generic transformation of coefficients of 
 * a medial model. This transformation has the form
 *   C' = T(C, P)
 * where C' and C are the model's coefficients (M-vector) and P are 
 * the parameters of the transformation.
 */
class CoefficientMapping
{
public:
  typedef vnl_vector<double> Vec;

  /** Apply the coefficient mapping C' = T(C, P) */
  virtual Vec Apply(const Vec &C, const Vec &p) = 0;

  /** Compute the variation J_T(P) * v_P given C, P and v_P */
  virtual Vec ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP) = 0;

  /** Compute the variation J_T(C) * v_C given C, P and v_C */
  virtual Vec ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC) = 0;

  /** Compute the variation corresponding to i-th parameter */
  virtual Vec GetVariationForParameter(const Vec &C, const Vec &P, size_t i)
    {
    Vec ei(n, 0.0); ei[i] = 1.0; 
    return ApplyJacobianInParameters(C, P, ei);
    }

  /** Get the number of coefficients in the coefficient space |C| */
  virtual size_t GetNumberOfCoefficients() const { return m; }

  /** Get the number of transformation parameters |Z| */
  virtual size_t GetNumberOfParameters() const { return n; }

  /** Whether this mapping is linear in P and C */
  virtual bool IsLinear() const { return true; }

protected:
  // The number of coefficients
  size_t m;

  // The number of parameters
  size_t n;
};

/**
 * This class represents a concatenation of two coefficient mappings 
 * h(C, Z_h) = f( g(C, Z_g), Z_f), where Z_h is formed by appending
 * Z_g to Z_f, ie, Z_h = [Z_f Z_g]
 */
class CompositionCoefficientMapping : public CoefficientMapping
{
public:
  typedef vnl_vector<double> Vec;

  CompositionCoefficientMapping(CoefficientMapping *f, CoefficientMapping *g)
    {
    assert(f->GetNumberOfCoefficients() == g->GetNumberOfCoefficients());
    this->f = f;
    this->g = g;
    this->nf = f->GetNumberOfParameters();
    this->ng = g->GetNumberOfParameters();
    this->n = ng + nf;
    this->m = f->GetNumberOfCoefficients();
    }

  // Destructor - virtual for child methods
  virtual ~CompositionCoefficientMapping() {};

  Vec Apply(const Vec &C, const Vec &Z) 
    { return f->Apply(g->Apply(C, Z.extract(ng,nf)), Z.extract(nf)); }

  /**
   * Expansion by chain rule:
   * U(X; Y,Z) = T(S(X;Z);Y)
   * J_U(y,z) * [v_y|v_z] = J_T(X) * J_S(Z) * v_z + J_T(Y_ * v_y;
   */
  Vec ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP) 
    {
    // Split the parameters and the variation into f and g components
    Vec Pf = P.extract(nf), Pg = P.extract(ng, nf);
    Vec varPf = varP.extract(nf), varPg = varP.extract(ng, nf);

    // Apply the inner transformation 
    Vec gCP = g->Apply(C, Pg);

    // Compute the combined variation
    return 
      f->ApplyJacobianInParameters(gCP, Pf, varPf) + 
      f->ApplyJacobianInCoefficients(gCP, Pf, g->ApplyJacobianInParameters(C, Pg, varPg));
    }

  /**
   * Expansion by chain rule:
   * U(X; Y,Z) = T(S(X;Z);Y)
   * J_U(x) * [v_x] = J_T(X) * J_S(X) * v_x;
   */
  Vec ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
    {
    // Split the parameters and the variation into f and g components
    Vec Pf = P.extract(nf), Pg = P.extract(ng, nf);

    // Apply the inner transformation 
    Vec gCP = g->Apply(C, Pg);

    // Combine the two transformations
    return f->ApplyJacobianInCoefficients(gCP, Pf, g->ApplyJacobianInCoefficients(C, Pg, varC));
    }

  /** The mapping is linear if its components are */
  bool IsLinear() const { return f->IsLinear() && g->IsLinear(); }

protected:
  size_t nf, ng;
  CoefficientMapping *f, *g;
};


/**
 * Coefficient mapping corresponding to the affine transform. Since medial
 * models are generic, we do not know a priori how the affine transform will
 * affect the coefficients of a medial model. Thus this class takes a pointer
 * to the model itself and calls virtual methods within the model
 */
class AffineTransformCoefficientMapping : public CoefficientMapping
{
public:
  typedef vnl_vector<double> Vec;

  AffineTransformCoefficientMapping(GenericMedialModel *model)
    {
    xTransformDescriptor = model->GetAffineTransformDescriptor();
    xCenter = xTransformDescriptor->GetCenterOfRotation(model->GetCoefficientArray());
    m = model->GetNumberOfCoefficients();
    n = 12;
    }

  /** Apply the coefficient mapping C' = T(C, P) */
  Vec Apply(const Vec &C, const Vec &p)
    {
    // Split the parameters into matrix A and vector b
    AffineMatrix A(p.data_block());
    AffineVector b = p.extract(3, 9);

    // Add an identity matrix to A (so that p is zero-based)
    A[0][0] += 1.0; A[1][1] += 1.0; A[2][2] += 1.0;

    // Apply the affine transform to the coefficients
    AffineTransformDescriptor::Mat AA(3,3);
    AA = A;
    return xTransformDescriptor->ApplyAffineTransform(C, AA, b, xCenter);
    }

  /** Compute the variation J_T(P) * v_P given C, P and v_P */
  Vec ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP)
    {
    // Split the variation into matrix A and vector b
    AffineMatrix varA(varP.data_block());
    AffineVector varb = varP.extract(3, 9);

    // Compute the corresponding variation using the model
    AffineTransformDescriptor::Mat AA = varA;
    return xTransformDescriptor->ApplyJacobianInParameters(C, AA, varb, xCenter);
    }

  /** Compute the variation J_T(C) * v_C given C, P and v_C */
  Vec ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
    {
    // Split the variation into matrix A and vector b
    AffineMatrix A(P.data_block());
    AffineVector b = P.extract(3, 9);

    // Add an identity matrix to A (so that p is zero-based)
    A[0][0] += 1.0; A[1][1] += 1.0; A[2][2] += 1.0;

    // Compute the corresponding variation using the model
    AffineTransformDescriptor::Mat AA = A;
    return xTransformDescriptor->ApplyJacobianInCoefficients(varC, AA, b, xCenter);
    }

  /** Affine transform is linear */
  bool IsLinear() const { return true; }

private:
  // Typedefs for the affine transform
  typedef vnl_matrix_fixed<double, 3, 3> AffineMatrix;
  typedef vnl_vector_fixed<double, 3> AffineVector;
  
  // The center of rotation in the affine transform
  AffineVector xCenter;
  
  // The model-specific descriptor of how the affine transform acts on the
  // model's coefficients
  const AffineTransformDescriptor *xTransformDescriptor;
};  

/** 
 * This class is used to select a subset of the coefficients for optimization.
 * The mapping is defined as C'[i] = C[i] + M[i] * P[j[i]], where M[i] is a
 * mask vector of zeros and ones and j[i] = Sum_{k=1}^j M[k].
 */
class SubsetCoefficientMapping : public CoefficientMapping
{
public:
  typedef vnl_vector<double> Vec;

  /** 
   * This constructor takes a mask of the selected coefficients. The mapping is
   * defined as follows:
   * j = 0; C'[i] = (mask[i] == 0) ? C[i] : C[i] + P[j++]
   */
  SubsetCoefficientMapping(const vnl_vector<size_t> &mask)
    {
    this->mask = mask;
    this->m = mask.size();
    this->n = m - std::count(mask.data_block(), mask.data_block() + m, 0);
    }

  Vec Apply(const Vec &C, const Vec &P) 
    {
    Vec X = C;
    for(size_t i = 0, j = 0; i < m; i++)
      if(mask[i] != 0)
        X[i] += P[j++];
    return X;
    }

  Vec ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP) 
    {
    Vec X(m);
    for(size_t i = 0, j = 0; i < m; i++)
      X[i] = (mask[i] == 0) ? 0.0 : varP[j++];
    return X;
    }

  Vec ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
    { return varC; }

  /** This mapping is linear */
  bool IsLinear() const { return true; }


private:
  vnl_vector<size_t> mask;  
};

/**
 * This is the identity mapping. It is only used when no other mapping is
 * defined. By this mapping, C'(C, P) = C + P.
 */
class IdentityCoefficientMapping : public CoefficientMapping
{
public:
  typedef vnl_vector<double> Vec;
  
  IdentityCoefficientMapping(size_t nCoefficients)
    { this->n = this->m = nCoefficients; }
  
  IdentityCoefficientMapping(GenericMedialModel *model)
    { this->n = this->m = model->GetNumberOfCoefficients(); }
  
  Vec Apply(const Vec &C, const Vec &P)
    { return C + P; }

  Vec ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP) 
    { return varP; }

  Vec ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC)
    { return varC; }

  /** This mapping is linear */
  bool IsLinear() const { return true; }
};

/** 
 * PCA Coefficient mapping. This mapping takes the form 
 * 
 * C' = C + Sqrt[L[i]] E[i] P[i]
 *
 * where L[i] are the PCE eigenvalues, E[i] are the eigenvectors and P[i] are
 * the parameters. This mapping is initialized with a PCA object. It is
 * expected that C is the mean of the PCA, but it does not have to be.
 */
class PCACoefficientMapping : public CoefficientMapping
{
public:
  typedef vnl_vector<double> Vec;
  
  PCACoefficientMapping(PrincipalComponents *pca, size_t nModes);
  Vec Apply(const Vec &C, const Vec &P);
  Vec ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP);
  Vec ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC);

  /** This mapping is linear */
  bool IsLinear() const { return true; }

private:
  PrincipalComponents *pca;
};

/**
 * A combined PCA/Affine mapping. In other words
 *
 * C' = AffineTransform[C + PCA[p1]], p2]
 */
class PCAPlusAffineCoefficientMapping : public CompositionCoefficientMapping
{
public:
  PCAPlusAffineCoefficientMapping(
    GenericMedialModel *model, PrincipalComponents *pca, size_t nModes);
  
  virtual ~PCAPlusAffineCoefficientMapping();
};

/** 
 * Coefficient mapping representing a reflection of a medial model about
 * a specified plane. This is only defined for non-branching medial models.
 * The control point pairs {(u,v),(u,-v)} are forced to be symmetric about
 * some plane P
 */
class ReflectionCoefficientMapping : public CoefficientMapping
{
public:
  typedef vnl_vector<double> Vec;

  /** This mapping is currently only defined for subdivision models. The plane
   * of reflection is defined as n.x - b = 0 */
  ReflectionCoefficientMapping(SubdivisionMedialModel *model, SMLVec3d &n, double b);

  /** Apply the coefficient mapping C' = T(C, P) */
  Vec Apply(const Vec &C, const Vec &p);

  /** Compute the variation J_T(P) * v_P given C, P and v_P */
  Vec ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP);

  /** Compute the variation J_T(C) * v_C given C, P and v_C */
  Vec ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC);

  /** Whether this mapping is linear in P and C */
  bool IsLinear() const { return true; }

protected:
  // The plane of reflection
  SMLVec3d p;
  double b;
  size_t k, nc;

  // The mapping from points to opposite points
  std::vector<size_t> opp;
};


class TriangleMesh;

/**
 * This class uses a set of N eigenfunctions of the LBO operator to 
 * create a basis for mesh deformation. The eigenfunctions are global,
 * orthonormal and increase in frequency.
 */
class MeshBasisCoefficientMapping : public CoefficientMapping
{
public:
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  /** Create the basis coefficient mapping */
  MeshBasisCoefficientMapping(
    const TriangleMesh *mesh, size_t basisSize, size_t nComponents);

  /** Apply the coefficient mapping C' = T(C, P) */
  Vec Apply(const Vec &C, const Vec &p);

  /** Compute the variation J_T(P) * v_P given C, P and v_P */
  Vec ApplyJacobianInParameters(const Vec &C, const Vec &P, const Vec &varP);

  /** Compute the variation J_T(C) * v_C given C, P and v_C */
  Vec ApplyJacobianInCoefficients(const Vec &C, const Vec &P, const Vec &varC);

  /** Access the basis */
  double GetBasisComponent(size_t iBasis, size_t iVertex)
    { return V[iBasis][iVertex]; }


protected:
  // The number of basis vectors
  size_t nb;

  // The number of components per vertex
  size_t nc;

  // The number of vertices
  size_t nv;

  // Eigenvectors of the Laplace operator
  
  Mat V;
};





#endif // __AffineMapping_h_

