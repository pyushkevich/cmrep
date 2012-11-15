#ifndef CONSTRAINEDCMREPOBJECTIVES_H
#define CONSTRAINEDCMREPOBJECTIVES_H

#include "GentleNLP.h"
#include "itkImage.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include <itkBSplineDecompositionImageFilter.h>
#include <itkConstantPadImageFilter.h>

// Add the vriables to optimize over

struct ConstrainedBoundaryCMRepProblemData
{
  // Coorfinates of the boundary, medial axis, normals, U-vectors
  gnlp::VarVecArray X, M, N, U;

  // Radius values
  gnlp::VarVec R;

  // Boundary and medial triangle areas
  gnlp::VarVec taX, taM;
};




// My own hack to allow b-spline interpolation of an image, derivatives and
// Hessian. It's only for 3D images.
class ImageJetInterpolator
{
public:
  typedef itk::Image<float, 3> ImageType;
  typedef vnl_vector_fixed<double, 3> Vec3;
  typedef vnl_matrix_fixed<double, 3, 3> Mat3;

  ImageJetInterpolator();

  void SetInputImage(ImageType *image, float outside_value);

  // Interpolate using voxel coordinates
  void EvaluateAtContinuousIndex(Vec3 x, double &f, Vec3 &G, Mat3 &H);

  // Evaluate using physical coordinates
  void EvaluateAtPhysicalPoint(Vec3 x, double &f, Vec3 &G, Mat3 &H);

protected:

  // A mini-pipeline for the image
  typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadFilter;
  typedef itk::BSplineDecompositionImageFilter<ImageType, ImageType> CoeffFilter;

  PadFilter::Pointer m_PadFilter;
  CoeffFilter::Pointer m_CoeffFilter;

  // Transformation from physical coords to voxel coords
  Mat3 m_A;
  Vec3 m_b;

};


// This is a traits class for image interpolation. For efficiency, this
// function caches the gradient value for a fixed set of x, y and z. A
// separate traits object must be created for every sample point!
class FloatImageTraits
{
public:

  FloatImageTraits(ImageJetInterpolator *source)
  {
    m_Source = source;
    m_Dirty = true;
    m_LastPos.fill(1e-100);
  }

  double Evaluate(double x, double y, double z,
                  int ordx, int ordy, int ordz)
  {
    if(m_Dirty || m_LastPos[0] != x || m_LastPos[1] != y || m_LastPos[2] != z)
      {
      m_LastPos = SMLVec3d(x,y,z);

      // Perform the interpolation
      m_Source->EvaluateAtPhysicalPoint(m_LastPos, m_LastF[0][0][0], m_G, m_H);

      // Copy the results
      m_LastF[1][0][0] = m_G[0];
      m_LastF[0][1][0] = m_G[1];
      m_LastF[0][0][1] = m_G[2];

      m_LastF[2][0][0] = m_H[0][0];
      m_LastF[0][2][0] = m_H[1][1];
      m_LastF[0][0][2] = m_H[2][2];

      m_LastF[1][1][0] = m_H[0][1];
      m_LastF[1][0][1] = m_H[0][2];
      m_LastF[0][1][1] = m_H[1][2];

      m_Dirty = false;
      }

    return m_LastF[ordx][ordy][ordz];
  }

protected:

  bool m_Dirty;
  SMLVec3d m_LastPos;
  double m_LastF[3][3][3];

  ImageJetInterpolator *m_Source;
  ImageJetInterpolator::Vec3 m_G;
  ImageJetInterpolator::Mat3 m_H;
};

void
CreateOverlapObjective(
    gnlp::ConstrainedNonLinearProblem *p,
    GenericMedialModel *model,
    ImageJetInterpolator *imsrc,
    gnlp::VarVecArray &X,
    gnlp::VarVecArray &M,
    gnlp::VarVecArray &U,
    gnlp::Expression **outObjectIntegral,
    gnlp::Expression **outVolumeIntegral,
    gnlp::VarVecArray &outSampleX,
    gnlp::VarVec &outSampleF);


void CreateAreaElementScaledOverlapObjective(
    gnlp::ConstrainedNonLinearProblem *p,
    GenericMedialModel *model,
    ImageJetInterpolator *imsrc,
    gnlp::VarVecArray &X,
    gnlp::VarVecArray &M,
    gnlp::VarVec &R,
    gnlp::VarVec &AeltX,
    gnlp::VarVec &AeltM,
    gnlp::Expression **outObjectIntegral,
    gnlp::Expression **outVolumeIntegral,
    gnlp::VarVecArray &outSampleX,
    gnlp::VarVec &outSampleF);

#endif // CONSTRAINEDCMREPOBJECTIVES_H
