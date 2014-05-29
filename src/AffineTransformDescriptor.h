#ifndef __AffineTransformDescriptor_h_
#define __AffineTransformDescriptor_h_

#include <vnl/vnl_det.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_trace.h>

class AffineTransformDescriptor
{
public:
  virtual ~AffineTransformDescriptor() {}

  // Vector and matrix definition
  typedef vnl_vector<double> Vector;
  typedef vnl_matrix_fixed<double, 3, 3> Mat33;
  typedef vnl_vector_fixed<double, 3> Vec3;

  /** Get center of rotation for affine transforms */
  virtual SMLVec3d GetCenterOfRotation(const Vector &C) const = 0;

  /** 
   * Apply the affine transform to the coefficients. That is, solve the
   * following problem for C':
   * 
   *   Model(C') = AffineTransform[ Model(C); A, b, ctr ];
   *
   * The solution depends on the definition of the model. Hence this method is
   * virtual and implemented differently for different surfaces.
   */
  virtual Vector ApplyAffineTransform(
    const Vector &C, const Mat33 &A, const Vec3 &b, const Vec3 &ctr) const = 0;

  /** 
   * Use the Jacobian of the transformation C'(A,b) to map a direction in the
   * affine transform space into the corresponding direction in the
   * coefficient space
   */
  virtual Vector ApplyJacobianInParameters(
    const Vector &C, const Mat33 &A, const Mat33 &dA, const Vec3 &db, const Vec3 &ctr) const = 0;

  /**
   * Use the Jacobian of the transformation C'(C) to map a direction in
   * the source coefficient space to target coefficient space
   */
  virtual Vector ApplyJacobianInCoefficients(
    const Vector &dC, const Mat33 &A, const Vec3 &b, const Vec3 &ctr) const = 0;
};

class FourierAffineTransformDescriptor : public AffineTransformDescriptor
{
public:
  virtual ~FourierAffineTransformDescriptor() {}
  typedef vnl_vector<double> Vector;
  typedef vnl_matrix_fixed<double, 3, 3> Mat33;
  typedef vnl_vector_fixed<double, 3> Vec3;

  /** 
   * Get the center of rotation given a coefficient vector. This is the first
   * coefficient C[0,0]
   */
  SMLVec3d GetCenterOfRotation(const Vector &C) const
    { return SMLVec3d(C.extract(3)); }

  /** 
   * Apply the transform to the coefficients. The rule is:
   *   C'[ij] = A * (C[ij]-c) + b + c        if i = j = 0
   *   C'[ij] = A * C[ij]                    o/w
   */
  Vector ApplyAffineTransform(const Vector &C, const Mat33 &A, const Vec3 &b, const Vec3 &ctr) const
    {
    Vector CPrime = C;

    // The first coefficient undergoes the complete affine transform
    CPrime.update(A * (C.extract(3, 0) - ctr) + ctr + b, 0);

    // For the rest of the coefficients, we only multiply by A
    for(size_t i = 4; i < C.size(); i+=4)
      CPrime.update(A * C.extract(3, i), i);

    return CPrime;
    }

  /**
   * Compute direction in coefficient space corresponding to a direction in
   * the affine transform space. 
   */
  Vector ApplyJacobianInParameters(
    const Vector &C, const Mat33 &A, const Mat33 &dA, const Vec3 &db, const Vec3 &ctr) const
    {
    Vector dCPrime(C.size(), 0.0);

    // Apply the variation in the first coefficient
    dCPrime.update(dA * (C.extract(3, 0) - ctr) + db, 0);

    // Apply the variation in the rest of the coefficients
    for(size_t i = 4; i < C.size(); i+=4)
      dCPrime.update(dA * C.extract(3, i), i);

    return dCPrime;
    }

  /** 
   * Compute direction in output coefficient space corresponding to a
   * direction in input coefficient space (v_C' = J(C'(C)) * v_C)
   */
  Vector ApplyJacobianInCoefficients(
    const Vector &dC, const Mat33 &A, const Vec3 &b, const Vec3 &ctr) const
    {
    // Start with C' = C
    Vector dCPrime = dC;

    // Apply the variation in the affine-transformed coefficients
    for(size_t i = 0; i < dC.size(); i+=4)
      dCPrime.update(A * dC.extract(3, i), i);

    return dCPrime;
    }
};

/**
 * This is an affine transform descriptor for models that are defined by a
 * list of points (control points in a subdivision surface, for example). Each
 * point is represented by four doubles.
 *
 * The center of rotation is the mean of the points. This may be a bad idea in
 * some applications, i.e., where the points are unevenly spaced over the
 * surface. However, the center of rotation is not that important, and simply
 * needs to be somewhere close to the object
 *
 * The constructor takes as a parameter the number of components per point
 *
 * In certain cases, the coefficients may contain values that should scale with
 * the size of the affine transform. For example, in the Biharmonic model, the 
 * radius along the boundary is one of the coefficients, and should scale with 
 * the size. The user can specify a list of coefficients that scale. 
 *
 * The scale is defined as the Jacobian of the affine transform to the power 1/3. 
 * This is somewhat arbitrary for the full affine transform, but at least it makes
 * sense for the similarity transform. 
 */
class PointArrayAffineTransformDescriptor : public AffineTransformDescriptor
{
public:
  virtual ~PointArrayAffineTransformDescriptor() {}
  typedef vnl_vector<double> Vector;
  typedef vnl_matrix_fixed<double, 3, 3> Mat33;
  typedef vnl_vector_fixed<double, 3> Vec3;
  typedef std::list<size_t> Index;

  typedef Index::const_iterator IndexCIt;

  /** Constructor, pass #components */
  PointArrayAffineTransformDescriptor(size_t nComponents)
    : AffineTransformDescriptor()
    {
    this->nc = nComponents;
    }

  /** 
   * Constructor, pass number of components, index of components that scale.
   * Each value should be between 0 and ncomponents - 1
   */
  PointArrayAffineTransformDescriptor(size_t nComponents, Index scalingComponents)
    : AffineTransformDescriptor()
    {
    this->nc = nComponents;
    this->si = scalingComponents;
    }

  /** 
   * Get the center of rotation given a coefficient vector. This is the first
   * mean of the points' XYZ components
   */
  SMLVec3d GetCenterOfRotation(const Vector &C) const
    { 
    SMLVec3d xCenter(0.0);
    for(size_t i = 0; i < C.size(); i+=nc)
      xCenter += C.extract(3, i);
    return xCenter * (nc * 1.0 / C.size());
    }

  /** 
   * Apply the transform to the coefficients. The rule is:
   *   C'[ij] = A * (C[ij] - ctr) + (ctr + b)
   */
  Vector ApplyAffineTransform(const Vector &C, const Mat33 &A, const Vec3 &b, const Vec3 &ctr) const
    {
    Vector CPrime = C;
    double scale = 1;

    // Compute ctr + b for speed
    SMLVec3d v = ctr + b;

    // Compute scale factor if needed
    if(si.size())
      {
      // std::cout << A << std::endl;
      double detA = vnl_det(A);
      scale = exp(log(fabs(detA)) / 3);
      // printf("DetA = %f, Scale = %f\n", detA, scale);
      }

    // Update the coefficients
    for(size_t i = 0; i < C.size(); i+=nc)
      {
      // For the coordinates, multiply by A
      CPrime.update(A * (C.extract(3, i) - ctr) + v, i);

      // For the scaling components, multiply by index
      for(IndexCIt it = si.begin(); it != si.end(); it++)
        {
        size_t c = *it;
        CPrime[i+c] = C[i+c] * scale;
        }
      }

    return CPrime;
    }

  /**
   * Compute direction in coefficient space corresponding to a direction in
   * the affine transform space. By differentiating the expression above we
   * have:
   *
   *   dC'[ij]/db[k]  = 1         if (j = k)
   *                  = 0         o/w
   *
   *   dC'[ij]/dA[kl] = C[il]     if j=k
   *                  = 0         o/w
   */
  Vector ApplyJacobianInParameters(
    const Vector &C, const Mat33 &A, const Mat33 &dA, const Vec3 &db, const Vec3 &ctr) const
    {
    Vector dCPrime(C.size(), 0.0);
    double scale = 1.0, dscale = 0;

    // If scaling, compute the derivative of the scale
    if(si.size())
      {
      double detA = vnl_det(A);
      scale = exp(log(fabs(detA)) / 3);
      if(scale > 0)
        dscale = scale * vnl_trace(vnl_inverse(A) * dA) / 3.0;
      }

    // Apply the variation in the rest of the coefficients
    for(size_t i = 0; i < C.size(); i+=nc)
      {
      // Update the coordinates
      dCPrime.update(dA * (C.extract(3, i) - ctr) + db, i);

      // Update the scaling components
      for(IndexCIt it = si.begin(); it != si.end(); it++)
        {
        size_t c = *it;
        dCPrime[i+c] = C[i+c] * dscale;
        }
      }

    return dCPrime;
    }

  /** 
   * Compute direction in output coefficient space corresponding to a
   * direction in input coefficient space (v_C' = J(C'(C)) * v_C)
   */
  Vector ApplyJacobianInCoefficients(
    const Vector &dC, const Mat33 &A, const Vec3 &b, const Vec3 &ctr) const
    {
    // Start out with C' = C
    Vector dCPrime = dC;
    double scale = 1.0;

    if(si.size())
      {
      double detA = vnl_det(A);
      scale = exp(log(fabs(detA)) / 3);
      }

    // Update the coefficients
    for(size_t i = 0; i < dC.size(); i+=nc)
      {
      // For the coordinates, multiply by A
      dCPrime.update(A * dC.extract(3, i), i);

      // For the scaling components, multiply by index
      for(IndexCIt it = si.begin(); it != si.end(); it++)
        {
        size_t c = *it;
        dCPrime[i+c] = dC[i+c] * scale;
        }
      }

    return dCPrime;
    }

private:
  size_t nc;
  Index si;

};

#endif // __AffineTransformDescriptor_h_
