#ifndef _MY_SML_MATH_
#define _MY_SML_MATH_

#include <vnl/vnl_cross.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_matops.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
// #include "align.h"

typedef vnl_vector_fixed<float,2> SMLVec2f;
typedef vnl_vector_fixed<float,3> SMLVec3f;
typedef vnl_vector_fixed<double,2> SMLVec2d;
typedef vnl_vector_fixed<double,3> SMLVec3d;
typedef vnl_vector_fixed<double,4> SMLVec4d;
typedef vnl_matrix_fixed<float,3,3> SMLMatrix3f;
typedef vnl_matrix_fixed<float,4,4> SMLMatrix4f;
typedef vnl_matrix_fixed<double,3,3> SMLMatrix3d;
typedef vnl_matrix_fixed<double,4,4> SMLMatrix4d;

class SMLVec4f : public vnl_vector_fixed<float,4> 
{
public:
  typedef vnl_vector_fixed<float,4> Superclass;
  
  SMLVec4f() : Superclass() {};
  SMLVec4f(float a, float b, float c, float d) :
    Superclass(a,b,c,d) {};
  SMLVec4f(float *x) : Superclass(x) {};

  SMLVec4f &operator =(const Superclass &src) 
    { *((Superclass *)this) = src; return *this; }
  
};

// Use a 4x4 matrix to transform a point in 3D euclidean space
template<class T>
inline vnl_vector_fixed<T,3> 
TransformVector(
  const vnl_matrix_fixed<T,4,4> &M,const vnl_vector_fixed<T,3> &V)
{
  return vnl_vector_fixed<T,3> 
    (M[0][0] * V[0] + M[0][1] * V[1] + M[0][2] * V[2],
     M[1][0] * V[0] + M[1][1] * V[1] + M[1][2] * V[2],
     M[2][0] * V[0] + M[2][1] * V[1] + M[2][2] * V[2]);
}

template<class T>
inline vnl_vector_fixed<T,3> 
TransformPoint(
  const vnl_matrix_fixed<T,4,4> &M,const vnl_vector_fixed<T,3> &V)
{
  return vnl_vector_fixed<T,3> 
    (M[0][0] * V[0] + M[0][1] * V[1] + M[0][2] * V[2] + M[0][3],
     M[1][0] * V[0] + M[1][1] * V[1] + M[1][2] * V[2] + M[1][3],
     M[2][0] * V[0] + M[2][1] * V[1] + M[2][2] * V[2] + M[2][3]);
}

inline double dot(const SMLVec3d &v1, const SMLVec3d &v2)
{
  return
    v1.data_block()[0] * v2.data_block()[0] +
    v1.data_block()[1] * v2.data_block()[1] +
    v1.data_block()[2] * v2.data_block()[2];
}

inline void vmuladd(SMLVec3d &t, const SMLVec3d &v, double s)
{
  t.data_block()[0] += v.data_block()[0] * s;
  t.data_block()[1] += v.data_block()[1] * s;
  t.data_block()[2] += v.data_block()[2] * s;
}

inline void vmulset(SMLVec3d &t, const SMLVec3d &v, double s)
{
  t.data_block()[0] = v.data_block()[0] * s;
  t.data_block()[1] = v.data_block()[1] * s;
  t.data_block()[2] = v.data_block()[2] * s;
}

template<class T>
inline T
squared_distance(const vnl_vector_fixed<T, 3> &A, const vnl_vector_fixed<T, 3> &B)
{
  double d0 = A[0] - B[0], d1 = A[1] - B[1], d2 = A[2] - B[2];
  return d0 * d0 + d1 * d1 + d2 * d2;
}




#endif // _MY_SML_MATH_
