#include "IPOptProblemInterface.h"
#include "ScriptInterface.h"
#include "BasisFunctions2D.h"
#include "MedialAtom.h"
#include "CartesianMedialModel.h"
#include "OptimizationTerms.h"
#include "CoefficientMapping.h"
#include "MedialAtomGrid.h"
#include "PrincipalComponents.h"
#include "System.h"
#include "TestSolver.h"
#include "ITKImageWrapper.h"
#include "MedialModelIO.h"
#include "IpIpoptApplication.hpp"
#include "MedialAtomGrid.h"
#include "vtkPolyDataWriter.h"
#include "tetgen.h"

#include "vtkPolyData.h"
#include "vtkCellLocator.h"
#include <vector>
#include <map>
#include <utility>
#include "itk_to_nifti_xform.h"

#include "ConstrainedCMRepObjectives.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <itkBSplineInterpolateImageFunction.h>
#include <itkBSplineDecompositionImageFilter.h>
#include <itkConstantPadImageFilter.h>


using namespace gnlp;



// ------------------------------------------------------------------------
// Interpolation code
// ------------------------------------------------------------------------

ImageJetInterpolator::ImageJetInterpolator()
{
  // Image gets padded to handle external values
  m_PadFilter = PadFilter::New();
  ImageType::SizeType szPad = {{3, 3, 3}};
  m_PadFilter->SetPadBound(szPad);

  // B-spline coefficient computation
  m_CoeffFilter = CoeffFilter::New();
  m_CoeffFilter->SetInput(m_PadFilter->GetOutput());
  m_CoeffFilter->SetSplineOrder(3);
}

void ImageJetInterpolator::SetInputImage(ImageType *image, float outside_value)
{
  m_PadFilter->SetInput(image);
  m_PadFilter->SetConstant(outside_value);
  m_PadFilter->Update();
  m_CoeffFilter->Update();

  // Compute the world to voxel transform
  vnl_matrix_fixed<double, 4, 4> voxel_to_world =
      ConstructNiftiSform(image->GetDirection().GetVnlMatrix(),
                          image->GetOrigin().GetVnlVector(),
                          image->GetSpacing().GetVnlVector());

  vnl_matrix_fixed<double, 4, 4> world_to_voxel = vnl_inverse(voxel_to_world);
  m_A = world_to_voxel.extract(3, 3, 0, 0);
  m_At = m_A.transpose();
  m_b = world_to_voxel.get_column(3).extract(3, 0);
}

void ImageJetInterpolator::EvaluateAtPhysicalPoint(
    Vec3 x, double &f, Vec3 &G, Mat3 &H)
{
  Vec3 xv = m_A * x + m_b;
  Vec3 Gv;
  Mat3 Hv;

  this->EvaluateAtContinuousIndex(xv, f, Gv, Hv);

  G = m_At * Gv;
  H = m_At * Hv * m_A;
}

// Interpolate image, gradient and Hessian using uniform B-spline interpolation
// The parameter X should be in voxel coordinates, not physical coordinates.
void ImageJetInterpolator::EvaluateAtContinuousIndex(
    Vec3 x, double &f, Vec3 &G, Mat3 &H)
{
  // Coordinate indices for each position in this 4x4x4 box
  static int qbox[3][64];
  static bool init_static = false;

  // Get the image
  ImageType *image = m_CoeffFilter->GetOutput();

  // Allocate an image to hold the 4x4 region of support
  if(!init_static)
    {
    // Populate the indices
    for(int i = 0; i < 64; i++)
      {
      qbox[0][i] = i % 4;
      qbox[1][i] = (i / 4) % 4;
      qbox[2][i] = (i / 16);
      }

    init_static = true;
    }

  // Get the dimensions of the image
  ImageType::IndexType ix = image->GetBufferedRegion().GetIndex();
  ImageType::SizeType sz = image->GetBufferedRegion().GetSize();

  // The desired 4x4 region of support from the source image
  ImageType::IndexType ixSrc;
  ImageType::SizeType szBox = {{4, 4, 4}};

  // As we are doing the bounds check, we can also extract the t for each
  // coordinate direction
  double ti[3];

  // Do the bounds check
  for(int i = 0; i < 3; i++)
    {
    // Get the start point and offset
    int xi = floor(x[i]);
    ti[i] = x[i] - xi;

    // These are the desired bounds
    ixSrc[i] = xi - 1;

    // Check if the box is to the left of the image
    int dleft = ix[i] - ixSrc[i];
    int dright = (ixSrc[i] + szBox[i]) - (ix[i] + sz[i]);

    // Reject out of bounds cases - bail!
    if(dleft > 0 || dright > 0)
      {
      f = m_PadFilter->GetConstant(); G.fill(0); H.fill(0);
      return;
      }
    }

  // Image and region on which we will operate
  ImageType::RegionType regMain;
  regMain.SetIndex(ixSrc);
  regMain.SetSize(szBox);

  // Now we need to compute the weights and derivatives for all three dirs.
  // The weight array is indexed by (0) derivative order; (1) coordinate axis;
  // (3) control point.
  double W[3][3][4];
  for(int i = 0; i < 3; i++)
    {
    double t = ti[i], t2 = t * t, t3 = t * t2;

    W[0][i][0] = -1.0 * t3 + 3.0 * t2 - 3.0 * t + 1.0;
    W[0][i][1] =  3.0 * t3 - 6.0 * t2           + 4.0;
    W[0][i][2] = -3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0;
    W[0][i][3] =  1.0 * t3;

    W[1][i][0] = -3.0 * t2 +  6.0 * t - 3.0;
    W[1][i][1] =  9.0 * t2 - 12.0 * t;
    W[1][i][2] = -9.0 * t2 +  6.0 * t + 3.0;
    W[1][i][3] =  3.0 * t2;

    W[2][i][0] = -6.0  * t +  6.0;
    W[2][i][1] =  18.0 * t - 12.0;
    W[2][i][2] = -18.0 * t +  6.0;
    W[2][i][3] =  6.0  * t;
    }

  // Clear outputs
  f = 0.0;
  G.fill(0.0);
  H.fill(0.0);

  // Now, multiply through
  itk::ImageRegionIteratorWithIndex<ImageType> it(image, regMain);
  for(int i = 0; i < 64; ++i, ++it)
    {
    double v = it.Get();



    int qx = qbox[0][i], qy = qbox[1][i], qz = qbox[2][i];
    double wx0 = W[0][0][qx], wx1 = W[1][0][qx], wx2 = W[2][0][qx];
    double wy0 = W[0][1][qy], wy1 = W[1][1][qy], wy2 = W[2][1][qy];
    double wz0 = W[0][2][qz], wz1 = W[1][2][qz], wz2 = W[2][2][qz];

    f += wx0 * wy0 * wz0 * v;

    G[0] += wx1 * wy0 * wz0 * v;
    G[1] += wx0 * wy1 * wz0 * v;
    G[2] += wx0 * wy0 * wz1 * v;

    H[0][0] += wx2 * wy0 * wz0 * v;
    H[0][1] += wx1 * wy1 * wz0 * v;
    H[0][2] += wx1 * wy0 * wz1 * v;
    H[1][1] += wx0 * wy2 * wz0 * v;
    H[1][2] += wx0 * wy1 * wz1 * v;
    H[2][2] += wx0 * wy0 * wz2 * v;
    }

  // Fill the rest of the Hessian
  H[2][1] = H[1][2]; H[1][0] = H[0][1]; H[2][0] = H[0][2];

  // Scale everything by 1/6^3
  const double scale_factor = 1.0 / 216.0;
  f *= scale_factor;
  G *= scale_factor;
  H *= scale_factor;
}





// ------------------------------------------------------------------------
// Image overlap objective code
// ------------------------------------------------------------------------

/**
  This function copies the implementation of overlap integral from old
  cm-rep code. That has some sort of a quadratic form for estimating the
  integral of intensity over a wedge. It seems to work, although I don't
  rememeber exactly where I got it from and why it works.
  */
void
CreateOverlapObjective(
    ConstrainedNonLinearProblem *p,
    GenericMedialModel *model,
    ImageJetInterpolator *imsrc,
    VarVecArray &X,
    VarVecArray &M,
    VarVecArray &U,
    Expression **outObjectIntegral,
    Expression **outVolumeIntegral,
    gnlp::VarVecArray &outSampleX,
    gnlp::VarVec &outSampleF)
{
  // Accumulators for the objective and model volume
  BigSum *xVolumeIntegral = new BigSum(p);
  BigSum *xObjectIntegral = new BigSum(p);

  // Prepare the output sample arrays
  outSampleX.clear();
  outSampleF.clear();

  // Set up constants
  double cEIGHTEENTH = 1.0 / 18.0;

  // Get the number of points
  unsigned int nb = model->GetNumberOfBoundaryPoints();

  // For each spoke, we define a vector of three coefficients for computing the
  // interior volume element.
  VarVecArray xInteriorVolumeElement(nb, VarVec(3, NULL));

  // We initialize each entry to an empty sum
  for(int i = 0; i < nb; i++)
    for(int j = 0; j < 3; j++)
      xInteriorVolumeElement[i][j] = new BigSum(p);

  // For each triangle, we compute its contribution to the volume element
  for(MedialBoundaryTriangleIterator ibt = model->GetBoundaryTriangleIterator();
      !ibt.IsAtEnd() ; ++ibt)
    {
    // Access the four medial atoms
    size_t ia0 = ibt.GetAtomIndex(0);
    size_t ia1 = ibt.GetAtomIndex(1);
    size_t ia2 = ibt.GetAtomIndex(2);
    size_t ib0 = ibt.GetBoundaryIndex(0);
    size_t ib1 = ibt.GetBoundaryIndex(1);
    size_t ib2 = ibt.GetBoundaryIndex(2);

    // Access the medial points and the boundary points
    VarVec X0 = M[ia0]; VarVec X1 = M[ia1]; VarVec X2 = M[ia2];
    VarVec U0 = U[ib0]; VarVec U1 = U[ib1]; VarVec U2 = U[ib2];

    // Difference vectors
    VarVec X10 = VectorApplyPairwise<BinaryDifference>(p, X1, X0);
    VarVec X20 = VectorApplyPairwise<BinaryDifference>(p, X2, X0);
    VarVec U10 = VectorApplyPairwise<BinaryDifference>(p, U1, U0);
    VarVec U20 = VectorApplyPairwise<BinaryDifference>(p, U2, U0);

    // Average of the U's
    VarVec W(3, NULL);
    for(int j = 0; j < 3; j++)
      W[j] = new ScalarProduct(p, new TernarySum(p,U0[j],U1[j],U2[j]), cEIGHTEENTH);

    // Normals to various faces(?)
    VarVec Za = CrossProduct(p, X10, X20);
    VarVec Zb = VectorApplyPairwise<BinarySum>(p,
                                               CrossProduct(p,U10,X20),
                                               CrossProduct(p,X10,U20));
    VarVec Zc = CrossProduct(p, U10, U20);

    // These are the volume element factors to distribute
    Expression *v[3] = { DotProduct(p, Za, W), DotProduct(p, Zb, W), DotProduct(p, Zc, W) };

    // Add new variables for the volume coefficients at this triangle
    for(int j = 0; j < 3; j++)
      {
      // Create a variable that equals v[j]
      Variable *vc = p->AddExpressionAsConstrainedVariable(v[j], "VJ");

      // Add that variable to each spoke
      static_cast<BigSum *>(xInteriorVolumeElement[ib0][j])->AddSummand(vc);
      static_cast<BigSum *>(xInteriorVolumeElement[ib1][j])->AddSummand(vc);
      static_cast<BigSum *>(xInteriorVolumeElement[ib2][j])->AddSummand(vc);
      }
    }

  // Create a new set of variables to represent per-vertex weights
  /*
  VarVecArray V(nb, VarVec(3, NULL));
  for(MedialBoundaryPointIterator bip = model->GetBoundaryPointIterator();
      !bip.IsAtEnd(); ++bip)
    {
    int i = bip.GetIndex();

    for(int j = 0; j < 3; j++)
      {
      // Create the variable
      V[i][j] = p->AddExpressionAsConstrainedVariable(xInteriorVolumeElement[i][j]);
      }
    }
  */


  // Now, compute the samples
  int nCuts = 3;
  int nSamplesPerAtom = nCuts + 2;

  double delta = 1.0 / (1.0 + nCuts);
  double hd = 0.5 * delta;

  // Compute the xi values associated with the samples
  vector<double> tSamples(nSamplesPerAtom, 0.0);
  for(size_t i = 0; i < nSamplesPerAtom; i++)
    tSamples[i] = i / (nCuts + 1.0);

  // Compute the coefficients for volume element computation at each
  // depth level (xi)
  vnl_matrix<double> xSampleCoeff(nSamplesPerAtom, 3);

  xSampleCoeff[0][0] = hd;
  xSampleCoeff[0][1] = hd * hd;
  xSampleCoeff[0][2] = hd * hd * hd;

  xSampleCoeff[nCuts+1][0] = hd;
  xSampleCoeff[nCuts+1][1] = hd * (1-hd);
  xSampleCoeff[nCuts+1][2] = hd*(1-hd)*(1-hd);

  for(int i = 1; i <= nCuts; i++)
    {
    double t = tSamples[i];
    xSampleCoeff[i][0] = delta;
    xSampleCoeff[i][1] = delta * t;
    xSampleCoeff[i][2] = delta * (t * t + hd * hd);
    }

  // Iterate over the boundary atoms in the model
  for(MedialBoundaryPointIterator bip = model->GetBoundaryPointIterator();
      !bip.IsAtEnd(); ++bip)
    {
    int ibnd = bip.GetIndex(), iatom = bip.GetAtomIndex();

    // The coefficient vector
    VarVec &vvec = xInteriorVolumeElement[ibnd];
    // VarVec &vvec = V[ibnd];

    // Compute the volume element and image value for the intermediate points
    for(int j = 0; j < nSamplesPerAtom; j++)
      {
      // Compute the volume element for this sample
      Expression *volumeElt = new TernarySum(
            p,
            new ScalarProduct(p, vvec[0], xSampleCoeff[j][0]),
            new ScalarProduct(p, vvec[1], xSampleCoeff[j][1]),
            new ScalarProduct(p, vvec[2], xSampleCoeff[j][2]));

      // Compute the sample point where to sample the image
      VarVec S(3, NULL);
      for(int k = 0; k < 3; k++)
        S[k] = new BinarySum(p, M[iatom][k], new ScalarProduct(p, U[ibnd][k], tSamples[j]));

      // Create a caching traits object for this sample (TODO: deletion!)
      FloatImageTraits *tr = new FloatImageTraits(imsrc);

      // Create an expression for sampling the function
      SampleFunction3<FloatImageTraits> *fS =
          new SampleFunction3<FloatImageTraits>(p, tr, S[0], S[1], S[2]);

      /*
      printf("*NEW* Spoke %04d:%04d, Sample %02d: VE=%8.4f\t IM=%8.4f\t X=%f, %f, %f\n",
             ibnd, iatom, j,
             volumeElt->Evaluate(),
             fS->Evaluate(),
             S[0]->Evaluate(), S[1]->Evaluate(), S[2]->Evaluate()); */

      // Compute the contribution to the total
      xVolumeIntegral->AddSummand(volumeElt);
      xObjectIntegral->AddSummand(new BinaryProduct(p, volumeElt, fS));

      // Store the samples for output
      outSampleX.push_back(S);
      outSampleF.push_back(fS);
      }
    }

  *outVolumeIntegral = xVolumeIntegral;
  *outObjectIntegral = xObjectIntegral;
}



// ------------------------------------------------------------------------
// Derive an image match objective
// ------------------------------------------------------------------------
// For each boundary point, obtain a list of samples. Scale these by a volume
// element, computed as a weighted sum of area elements and radius
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
    gnlp::VarVec &outSampleF)
{
  // Accumulators for the objective and model volume
  BigSum *xVolumeIntegral = new BigSum(p);
  BigSum *xObjectIntegral = new BigSum(p);

  // Sample along each spoke
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int iBnd = it.GetIndex();
    int iAtom = it.GetAtomIndex();

    // The sum of samples along the spoke
    BigSum *spokeIntegralImage = new BigSum(p);

    // Obtain expressions at samples along the normal vector
    for(double delta = 0.0; delta <= 1.0; delta += 0.25)
      {
      // Interpolate the sample
      VarVec S(3, NULL);
      for(int j = 0; j < 3; j++)
        {
        S[j] = new BinarySum(
              p,
              new ScalarProduct(p, M[iAtom][j], 1.0-delta),
              new ScalarProduct(p, X[iBnd][j], delta));
        }

      // Create a traits object for this sample
      FloatImageTraits *tr = new FloatImageTraits(imsrc);

      // Create an expression for sampling the function
      SampleFunction3<FloatImageTraits> *fS =
          new SampleFunction3<FloatImageTraits>(p, tr, S[0], S[1], S[2]);

      // Store the samples
      outSampleX.push_back(S);
      outSampleF.push_back(fS);

      // Interpolate the area element. This is obviously an approximation, since
      // the correct form of the area element involves computing the triangle
      // areas for this sample. We don't want to do this for computational reasons
      Expression *aelt = new BinarySum(
            p,
            new ScalarProduct(p, AeltM[iAtom], 1-delta),
            new ScalarProduct(p, AeltX[iBnd], delta));

      // Compute the contribution to the trapesoid rule. It's a factor of R/2N
      // for endpoints and R/N for internal points. However, the factor of N
      // can be handled later.
      Expression *aeltWgt = (delta > 0.001 && delta < 0.999)
          ? aelt : new ScalarProduct(p, aelt, 0.5);

      // Accumulate the weighted sum for this spoke
      spokeIntegralImage->AddSummand(new BinaryProduct(p, fS, aeltWgt));
      }

    // Now we have a value along the spoke. We still have to scale it by
    // several factors: the length of the spoke (R), the sample interval
    // width, and the factor of 3 in the area element expressions.
    Expression *intImageFinal = new TernaryProduct(
          p, new Constant(p, 0.25 / 3.0), R[iAtom], spokeIntegralImage);

    // Finally, the sum of the area element weight should just be
    // (aEltB+aEltM)/2 * R, but again, the factor of 3
    Expression *intVolumeFinal = new TernaryProduct(
          p, new Constant(p, 1.0 / 6.0), R[iAtom],
          new BinarySum(p, AeltM[iAtom], AeltX[iBnd]));

    // Add it all up
    xObjectIntegral->AddSummand(intImageFinal);
    xVolumeIntegral->AddSummand(intVolumeFinal);
    }

  *outVolumeIntegral = xVolumeIntegral;
  *outObjectIntegral = xObjectIntegral;
}


