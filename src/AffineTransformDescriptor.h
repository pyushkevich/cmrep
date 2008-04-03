#ifndef __AffineTransformDescriptor_h_
#define __AffineTransformDescriptor_h_

class AffineTransformDescriptor
{
public:
  // Vector and matrix definition
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  /** Get center of rotation for affine transforms */
  virtual SMLVec3d GetCenterOfRotation(const Vec &C) const = 0;

  /** 
   * Apply the affine transform to the coefficients. That is, solve the
   * following problem for C':
   * 
   *   Model(C') = AffineTransform[ Model(C); A, b, ctr ];
   *
   * The solution depends on the definition of the model. Hence this method is
   * virtual and implemented differently for different surfaces.
   */
  virtual Vec ApplyAffineTransform(
    const Vec &C, const Mat &A, const Vec &b, const Vec &ctr) const = 0;

  /** 
   * Use the Jacobian of the transformation C'(A,b) to map a direction in the
   * affine transform space into the corresponding direction in the
   * coefficient space
   */
  virtual Vec ApplyJacobianInParameters(
    const Vec &C, const Mat &dA, const Vec &db, const Vec &ctr) const = 0;

  /**
   * Use the Jacobian of the transformation C'(C) to map a direction in
   * the source coefficient space to target coefficient space
   */
  virtual Vec ApplyJacobianInCoefficients(
    const Vec &dC, const Mat &A, const Vec &b, const Vec &ctr) const = 0;
};

class FourierAffineTransformDescriptor : public AffineTransformDescriptor
{
public:
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  /** 
   * Get the center of rotation given a coefficient vector. This is the first
   * coefficient C[0,0]
   */
  SMLVec3d GetCenterOfRotation(const Vec &C) const
    { return SMLVec3d(C.extract(3)); }

  /** 
   * Apply the transform to the coefficients. The rule is:
   *   C'[ij] = A * (C[ij]-c) + b + c        if i = j = 0
   *   C'[ij] = A * C[ij]                    o/w
   */
  Vec ApplyAffineTransform(const Vec &C, const Mat &A, const Vec &b, const Vec &ctr) const
    {
    Vec CPrime = C;

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
  Vec ApplyJacobianInParameters(
    const Vec &C, const Mat &dA, const Vec &db, const Vec &ctr) const
    {
    Vec dCPrime(C.size(), 0.0);

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
  Vec ApplyJacobianInCoefficients(
    const Vec &dC, const Mat &A, const Vec &b, const Vec &ctr) const
    {
    // Start with C' = C
    Vec dCPrime = dC;

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
 */
class PointArrayAffineTransformDescriptor : public AffineTransformDescriptor
{
public:
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  /** 
   * Get the center of rotation given a coefficient vector. This is the first
   * mean of the points' XYZ components
   */
  SMLVec3d GetCenterOfRotation(const Vec &C) const
    { 
    SMLVec3d xCenter(0.0);
    for(size_t i = 0; i < C.size(); i+=4)
      xCenter += C.extract(3, i);
    return xCenter * (4.0 / C.size());
    }

  /** 
   * Apply the transform to the coefficients. The rule is:
   *   C'[ij] = A * (C[ij] - ctr) + (ctr + b)
   */
  Vec ApplyAffineTransform(const Vec &C, const Mat &A, const Vec &b, const Vec &ctr) const
    {
    Vec CPrime = C;

    // Compute ctr + b for speed
    SMLVec3d v = ctr + b;

    // For the rest of the coefficients, multiply by A
    for(size_t i = 0; i < C.size(); i+=4)
      CPrime.update(A * (C.extract(3, i) - ctr) + v, i);

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
  Vec ApplyJacobianInParameters(
    const Vec &C, const Mat &dA, const Vec &db, const Vec &ctr) const
    {
    Vec dCPrime(C.size(), 0.0);

    // Apply the variation in the rest of the coefficients
    for(size_t i = 0; i < C.size(); i+=4)
      dCPrime.update(dA * (C.extract(3, i) - ctr) + db, i);

    return dCPrime;
    }

  /** 
   * Compute direction in output coefficient space corresponding to a
   * direction in input coefficient space (v_C' = J(C'(C)) * v_C)
   */
  Vec ApplyJacobianInCoefficients(
    const Vec &dC, const Mat &A, const Vec &b, const Vec &ctr) const
    {
    // Start out with C' = C
    Vec dCPrime = dC;

    // Apply the variation in the rest of the coefficients
    for(size_t i = 0; i < dC.size(); i+=4)
      dCPrime.update(A * dC.extract(3, i), i);

    return dCPrime;
    }
};

#endif // __AffineTransformDescriptor_h_
