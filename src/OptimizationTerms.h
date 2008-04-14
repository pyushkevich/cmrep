#ifndef __OptimizationTerms_h_
#define __OptimizationTerms_h_

#include <smlmath.h>
#include "CodeTimer.h"
#include "GenericMedialModel.h"
#include "BasisFunctions2D.h"
#include "ScriptInterface.h"
#include "CoefficientMapping.h"
#include "SmoothedImageSampler.h"
#include "MedialAtomGrid.h"
#include "Registry.h"

using medialpde::FloatImage;

/******************************************************************
 * THIS EUCLIDEAN FUNCTION JUNK SHOULD GO SOMEWHERE ELSE
 * ************************************************************** */

/**
 * A simple image adapter that samples the image without applying
 * any additional processing
 */
class FloatImageEuclideanFunctionAdapter : public EuclideanFunction
{
public:
  FloatImageEuclideanFunctionAdapter(FloatImage *image)
    { this->image = image; }

  double Evaluate(const SMLVec3d &X)
    { return (double) image->Interpolate(X); }
  
  void ComputeGradient(const SMLVec3d &X, SMLVec3d &G)
    { image->InterpolateImageGradient(X, G); }
  
protected:
  FloatImage *image;
};

/** A function that computes the match as (I(x) - I0)^2, where I0 is some
 * level set in the image. Used for blurred binary images */
class FloatImageSquareValueFunctionAdapter : public EuclideanFunction
{
public:
  FloatImageSquareValueFunctionAdapter(FloatImage *image, float xLevelSet = 0.0f) 
    { this->image = image; this->xLevelSet = xLevelSet; }

  double Evaluate(const SMLVec3d &X)
    { 
    // Get the image match from the superclass
    double I = (double) image->Interpolate(X);
    double d = I - xLevelSet; 
    return d * d;
    }

  void ComputeGradient(const SMLVec3d &X, SMLVec3d &G)
    {
    // Get the gradient and match from the superclass
    double I = (double) image->Interpolate(X);
    double d = I - xLevelSet;
    image->InterpolateImageGradient(X, G);
    G *= 2.0 * d;
    }

private:
  FloatImage *image;
  float xLevelSet;
};

class ImageSmoothSamplingEuclideanFunction : public EuclideanFunction
{
public:
  ImageSmoothSamplingEuclideanFunction(
    FloatImage *image,
    double sigma, double cutoff = 3.5);

  ~ImageSmoothSamplingEuclideanFunction()
    { delete sis; }

  double ComputeFunctionAndGradient(const SMLVec3d &X, SMLVec3d &G)
    {
    SMLVec3d Ximg;
    for(size_t d = 0; d < 3; d++)
     Ximg[d] = (X[d] - xOrigin[d]) * xInvVoxSize[d];

    double f = sis->Sample(Ximg.data_block(), G.data_block());

    for(size_t d = 0; d < 3; d++)
     G[d] = G[d] * xInvVoxSize[d];

    return f;
    }

  double Evaluate(const SMLVec3d &X)
    { 
    SMLVec3d Ximg;
    for(size_t d = 0; d < 3; d++)
     Ximg[d] = (X[d] - xOrigin[d]) * xInvVoxSize[d];

    return sis->Sample(Ximg.data_block());
    }

  void ComputeGradient(const SMLVec3d &X, SMLVec3d &G)
    {
    ComputeFunctionAndGradient(X, G);
    }

private:
  FloatImage *image;
  SmoothedImageSampler *sis;

  SMLVec3d xOrigin, xInvVoxSize;

};


/** 
 * A base class for data associated with each solution as well as the 
 * data associated with the partial derivatives of solutions
 */
class SolutionDataBase
{
public:
  // Constructor. One of the parameters is the number of cuts for internal
  // point representation
  SolutionDataBase(MedialIterationContext *xGrid);
  virtual ~SolutionDataBase() {};
  
  // The medial model that specifies the grid/mesh on which atoms live
  MedialIterationContext *xAtomGrid;

  // The array of atoms in this solution. This array can either point to the 
  // xMedialModel's array or it can be pointing to another set of atoms
  MedialAtom *xAtoms;

  // Compute integration weights. This involves computing area-based weights
  // on the medial surface, as well as weights needed to integrate over the
  // interior of the model
  virtual void ComputeIntegrationWeights() = 0;

  // The array of boundary weights (area elements) and the total bnd. area
  std::vector<double> xBoundaryWeights;
  std::vector<SMLVec3d> xBoundaryAreaVector;
  std::vector<SMLVec3d> xBoundaryTriangleUnitNormal;
  std::vector<double> xBoundaryTriangleArea;
  double xBoundaryArea;

  // The array of medial weights (area elements) and the total bnd. area
  std::vector<double> xMedialWeights;
  std::vector<SMLVec3d> xMedialTriangleUnitNormal;
  double xMedialArea;

  // Array of values for computing volume element. The volume element is
  // a quadratic form a + b t + c t^2 and this array holds the ABCs as vectors
  std::vector<SMLVec3d> xInteriorVolumeElement;

  size_t nAtoms, nBndPts, nMedTri, nBndTri;
};

/** 
 * A structure to be passed in and out of the energy term functions
 */
class SolutionData : public SolutionDataBase
{
public:

  // Initialize the solution data using a grid and an array of derivative
  // atoms
  SolutionData(MedialIterationContext *xGrid, MedialAtom *xAtoms);

  void ComputeIntegrationWeights();

  // Computes the volume element at atom i, value xi, given that the sampling
  // interval in xi is delta. The last flag specifies whether xi is at the end
  // of the spoke, in which case only the volume 'on the inside' is computed
  // (SAME CODE AS USED IN VolumeIntegralEnergyTerm)
  double ComputeVolumeElement(size_t iAtom, double xi, double delta, bool edge)
    {
    // Find the right boundary atom
    size_t iSide = (xi < 0) ? 0 : 1;
    size_t iBnd = xAtomGrid->GetBoundaryPointIndex(iAtom, iSide);

    // Compute the volume element
    double xi1 = abs(xi) - 0.5 * delta, xi2 = abs(xi) + 0.5 * delta;
    SMLVec3d t1(1, xi1, xi1 * xi1), t2(1, xi2, xi2 * xi2);
    if(edge)
      return 0.5 * delta * dot_product(xInteriorVolumeElement[iBnd], t1);
    else
      return 0.5 * delta * dot_product(xInteriorVolumeElement[iBnd], t1 + t2);
    }

};

class PartialDerivativeSolutionData : public SolutionDataBase
{
public:

  // Initialize the solution data using a grid
  PartialDerivativeSolutionData(
    SolutionData *xReference, MedialAtom *dAtoms);

  void ComputeIntegrationWeights();

private:
  SolutionData *xReference;
};

/** 
 * Accumulator for statistics in energy terms
 */
class StatisticsAccumulator
{
public:
  StatisticsAccumulator() { Reset(); }

  void Reset()
    { n = 0; xMin = xMax = xSum = xSumSq = 0.0; }

  void Update(double x)
    {
    // Update Min and Max
    if(n == 0)
      { xMin = xMax = x; }
    else
      { xMin = std::min(xMin, x); xMax = std::max(xMax, x); }

    // Update the rest
    xSum += x;
    xSumSq += x * x;
    n++;
    }

  double GetMin() const { return xMin; }
  double GetMax() const { return xMax; }
  double GetMaxAbs() const { return std::max(fabs(xMin), fabs(xMax)); }
  double GetMean() const { return xSum / n; }
  double GetVariance() const { return (xSumSq - xSum * xSum / n) / (n-1); }
  double GetStdDev() const { return sqrt(GetVariance()); }
  double GetSum() const { return xSum; }
  size_t GetCount() const { return n; }

private:
  size_t n;
  double xMin, xMax, xSum, xSumSq;
};

inline ostream & operator << (ostream &sout, const StatisticsAccumulator &S)
{
  sout << "mean(sd): " << S.GetMean() << "(" << S.GetStdDev() << "); ";
  sout << "range: " << S.GetMin() << " to " << S.GetMax() << "; ";
  sout << "n: " << S.GetCount() << ";";
  return sout;
}

class EnergyTerm 
{
public:
  // Compute the energy term 
  virtual double ComputeEnergy(SolutionData *data) = 0;

  // Initialize gradient computation and return the value of the solution
  // at the current state
  virtual double BeginGradientComputation(SolutionData *SCenter)
    { return ComputeEnergy(SCenter); }

  // Compute the partial derivative of the energy function
  virtual double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS) = 0;

  // Finish gradient computation, remove all temporary data
  virtual void EndGradientComputation() {};

  // Print a verbose report
  virtual void PrintReport(ostream &sout) = 0;

  // Print a short name
  virtual string GetShortName() = 0;

  // Pass in parameters using a registry object
  virtual void SetParameters(Registry &r) {}
};

/**
 * This is a special energy term that involves integrating stuff
 * over the medial axis, ie. \int f(u,v) du dv. For such a term, we
 * want the du dv quantity, which is non-uniform over the mesh in many
 * cases (because of non-uniform grid spacing)
 */
class MedialIntegrationEnergyTerm : public EnergyTerm
{
public:
  MedialIntegrationEnergyTerm(GenericMedialModel *model);

protected:
  double xDomainArea;
  vnl_vector<double> xDomainWeights;
};

class BoundaryImageMatchTerm : public EnergyTerm
{
public:
  // Constructor
  BoundaryImageMatchTerm(GenericMedialModel *model, FloatImage *image);
    ~BoundaryImageMatchTerm();

  // Compute the image match
  double ComputeEnergy(SolutionData *data)
    { return UnifiedComputeEnergy(data, false); }

  // Initialize gradient computation and return the value of the solution
  // at the current state
  double BeginGradientComputation(SolutionData *data)
    { return UnifiedComputeEnergy(data, true); }
  
  // Compute the partial derivative (must be called in the middle of Begin and
  // End of GradientComputation.
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a verbose report
  void PrintReport(ostream &sout);

  // Print a short name
  string GetShortName() { return string("BMATCH"); }

private:
  FloatImage *xImage;

  // Common energy function
  double UnifiedComputeEnergy(SolutionData *, bool);

  // Terms used in reporting details
  double xImageMatch, xBoundaryArea, xFinalMatch;

  // Temporaries for gradient computation
  SMLVec3d *xGradI;
  double *xImageVal;
};

class LocalDistanceDifferenceEnergyTerm : public EnergyTerm
{
public:
  /** Standard initialization */
  LocalDistanceDifferenceEnergyTerm(GenericMedialModel *model);

  // Set the parameters from a registry
  void SetParameters(Registry &r);

  /** Compute the volume overlap fraction between image and m-rep */
  double ComputeEnergy(SolutionData *data)
    { return UnifiedComputeEnergy(data, false); }

  // Initialize gradient computation and return the value of the solution
  double BeginGradientComputation(SolutionData *data)
    { return UnifiedComputeEnergy(data, true); }
  
  // Finish gradient computation, remove all temporary data
  void EndGradientComputation() {};

  /** Compute directional deriavative */
  double ComputePartialDerivative(SolutionData *S, 
    PartialDerivativeSolutionData *dS);

  /** Print the report */
  void PrintReport(ostream &sout);

  /** Print a short name of the energy term */
  string GetShortName() { return string("LOCDIS"); }

private:

  /** Compute x-correlation */
  double UnifiedComputeEnergy(SolutionData *S, bool gradient_mode);

  /** Sample the reference image */
  void InitializeReferenceData(
    vector<GenericMedialModel *> &mdlReference);

  // Structure that holds intensity profile data
  struct DistanceData
    {
    double dm[3], db0[3], db1[3], r[3];
    };

  // Medial model
  GenericMedialModel *xModel;

  // Array of profiles, one for each boundary point
  std::vector<DistanceData> xRefDistData, xDistData;

  // Statistical accumulators
  StatisticsAccumulator saDist, saDistDeriv;

  // Total penalty
  double xTotalPenalty;
};

/**
   * Cross correlation energy term is initialized with a target image,
   * a medial model, a reference image, and a reference medial model.
   * The target image is the one we want to match, and the model is the
   * deformable model that will be trying to match the target image. The
   * reference model and image are the ones against which cross-correlation
   * will be computed.
   *
   * The last two parameters are nCuts (the number of sample points along 
   * each spoke, and xiMax, the extent of the spoke, which is 1.0 by default,
   * meaning that the spoke extends to the boundary. Setting xiMax to 2.0 will
   * make it extend one radius value past the boundary.
   */
class CrossCorrelationImageMatchTerm : public EnergyTerm
{
public:
  
  /** Initialize this term with the current model and target image*/
  CrossCorrelationImageMatchTerm(
    GenericMedialModel *model, FloatImage *xGrayImage);

  // Destructor
  ~CrossCorrelationImageMatchTerm();

  /** Load the parameters from file */
  void SetParameters(Registry &r);

  /** Compute the volume overlap fraction between image and m-rep */
  double ComputeEnergy(SolutionData *data)
    { return UnifiedComputeEnergy(data, false); }

  // Initialize gradient computation and return the value of the solution
  double BeginGradientComputation(SolutionData *data)
    { return UnifiedComputeEnergy(data, true); }
  
  // Finish gradient computation, remove all temporary data
  void EndGradientComputation() {};

  /** Compute directional deriavative */
  double ComputePartialDerivative(SolutionData *S, 
    PartialDerivativeSolutionData *dS);

  /** Print the report */
  void PrintReport(ostream &sout);

  /** Print a short name of the energy term */
  string GetShortName() { return string("X-CORR"); }

  EuclideanFunction *GetReferenceFunction()
    { return fReference; }

  EuclideanFunction *GetTargetFunction()
    { return fTarget; }

private:

  /** Compute x-correlation */
  double UnifiedComputeEnergy(SolutionData *S, bool gradient_mode);

  /** Sample the reference image */
  void InitializeReferenceData(
    EuclideanFunction *fReference, GenericMedialModel *mdlReference);

  // Structure that holds intensity profile data
  struct ProfileData
    {
    // Image profile and gradient
    std::vector<SMLVec3d> xImageGrad;
    std::vector<double> xImageVal;

    // Reference image profile
    std::vector<double> xRefImgVal;

    // Statistics from the reference image
    double xRefMean, xRefSD;
    double xSD, xCov, xMean, xCorr;

    ProfileData(size_t n) 
      : xImageGrad(n), xImageVal(n), xRefImgVal(n) {}
    ProfileData() {}
    };

  // Array of profiles, one for each boundary point
  std::vector<ProfileData> xProfile;

  // Array of sample points in xi
  std::vector<double> xSamples;

  // Sampling info
  size_t nSamplesPerAtom;

  // Statistical accumulators
  StatisticsAccumulator saPenalty, saCorr;

  // Total penalty
  double xTotalPenalty;

  // Model
  GenericMedialModel *xModel;

  // Target image function
  EuclideanFunction *fTarget, *fReference;
};

/**
 * This term computes the integral of some function over the volume of a 
 * cm-rep. 
 */
class VolumeIntegralEnergyTerm : public EnergyTerm
{
public:

  /** Initialize with an image and a number of sample points on each
   * medial sail vector */
  VolumeIntegralEnergyTerm(
    GenericMedialModel *model, EuclideanFunction *function, size_t nCuts);
  
  /** Compute the volume overlap fraction between image and m-rep */
  double ComputeEnergy(SolutionData *data)
    { return UnifiedComputeEnergy(data, false); }

  // Initialize gradient computation and return the value of the solution
  // at the current state
  double BeginGradientComputation(SolutionData *data)
    { return UnifiedComputeEnergy(data, true); }
  
  // Compute the partial derivative (must be called in the middle of Begin and
  // End of GradientComputation.
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Finish gradient computation, remove all temporary data
  void EndGradientComputation() {};

  // Print a verbose report
  void PrintReport(ostream &sout);

  double GetModelVolume()
    { return xVolumeIntegral; }

  // Print a short name
  string GetShortName() { return string("VOLOVL"); }
  
private:
  /** Common method for energy computation */
  double UnifiedComputeEnergy(SolutionData *S, bool gradient_mode);

  // Image object to sample
  EuclideanFunction *function;

  // Cumulative values
  double xImageIntegral, xVolumeIntegral, xObjectIntegral, xRatio;

  // Profile array sizes
  size_t nCuts, nSamplesPerAtom, iCenterSample, iLastSample, nAtoms;

  // Samples (xi coordinate)
  std::vector<double> xSamples;
  std::vector<SMLVec3d> xSampleCoeff;

  // Structure that holds intensity profile data
  struct ProfileData
    {
    std::vector<SMLVec3d> xImageGrad;
    std::vector<double> xImageVal;
    std::vector<double> xVolumeElt;
    double r_scale;
    size_t bix1, bix2;
    ProfileData(size_t n) 
      : xImageGrad(n), xImageVal(n), xVolumeElt(n) {}
    ProfileData() {};
    };
  std::vector<ProfileData> xProfile;

  friend class VolumeOverlapEnergyTerm;
};

/**
 * This term computes the amount of volume overlap between an image and the
 * m-rep. The image represents the object of interest with voxels of intensity
 * 1.0 and background with intensity 0.0. The image should have floating point
 * voxels, where values between 0 and 1 represent partial volume or
 * uncertainty
 */
class VolumeOverlapEnergyTerm : public EnergyTerm
{
public:

  /** Initialize with an image and a number of sample points on each
   * medial sail vector */
  VolumeOverlapEnergyTerm(
    GenericMedialModel *model, FloatImage *image, size_t nCuts);

  ~VolumeOverlapEnergyTerm()
    { delete worker; delete function; }
  
  /** Compute the volume overlap fraction between image and m-rep */
  double ComputeEnergy(SolutionData *data);

  // Initialize gradient computation and return the value of the solution
  // at the current state
  double BeginGradientComputation(SolutionData *data);
  
  // Compute the partial derivative (must be called in the middle of Begin and
  // End of GradientComputation.
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Finish gradient computation, remove all temporary data
  void EndGradientComputation() {};

  // Print a verbose report
  void PrintReport(ostream &sout);


  // Print a short name
  string GetShortName() { return string("VOLOVL"); }

  double GetModelVolume()
    { return worker->GetModelVolume(); }
  
private:
  // Image object to sample
  FloatImage *xImage;
  
  // Image encoded as a function
  EuclideanFunction *function;

  // Cumulative values
  double xImageIntegral, xRatio;

  // The class that does the actual work
  VolumeIntegralEnergyTerm *worker;
};

class BoundaryJacobianEnergyTerm : public EnergyTerm
{
public:
  BoundaryJacobianEnergyTerm();

  // Compute the energy
  double ComputeEnergy(SolutionData *data);

  // Print a verbose report
  void PrintReport(ostream &sout);

  // Compute the partial derivative term
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Get the worst Jacobian
  double GetMinJacobian() { return saJacobian.GetMin(); }
  
  /** Pass in parameters from a registry */
  void SetParameters(Registry &R)
    {
    xPenaltyA = R["PenaltyA"][xDefaultPenaltyA];
    xPenaltyB = R["PenaltyB"][xDefaultPenaltyB];
    }

  // Print a short name
  string GetShortName() { return string("BNDJAC"); }

private:
  // Maximum and minimum values of the Jacobian encountered during
  // optimization
  StatisticsAccumulator saJacobian, saPenalty;

  // double xMinJacobian, xMaxJacobian, xAvgJacobian, xTotalPenalty;

  // In addition to the total values, we keep track of per-quad values
  // computed in the course of calculating the penalty. This allows us
  // to compute derivatives more quickly
  struct TriangleEntry {
    SMLVec3d XU, XV, YU[2], YV[2], NX, NY[2];
    double gX2, J[2], PenA[2], PenB[2];
  };
  
  typedef vector<TriangleEntry> TriangleVector;
  TriangleVector xTriangleEntries;

  /** Penalty function applied to the squared jacobian */
  double PenaltyFunction(double x, double a, double b)
    { return exp(-a * x) + exp(x - b); }

  /** The derivative of the penalty function */
  double PenaltyFunctionDerivative(double x, double a, double b)
    { return exp(x - b) - a * exp (-a * x); }

  // Constants used in the penalty computation
  double xPenaltyA, xPenaltyB;
  const static double xDefaultPenaltyA, xDefaultPenaltyB;
};

/**
 * This heuristic penalty punishes angles below a certain minimum (e.g., 20 degrees)
 * with a quadratic penalty. It should lead to quality triangulation. However, one 
 * must be aware that limiting the triangulation quality will also limit how well a
 * model can fit data. Ideally, this will be used with a remeshing algorithm
 */
/*
class MedialMinimumTriangleAnglePenaltyTerm : public EnergyTerm
{
public:
  // Initialize the term
  MedialMinimumTriangleAnglePenaltyTerm(GenericMedialModel *model);

  // Compute the energy
  double ComputeEnergy(SolutionData *data)
    { return this->UnifiedComputeEnergy(data, false); }

  // Compute the energy and begin gradient
  double BeginGradientComputation(SolutionData *data)
    { return this->UnifiedComputeEnergy(data, true); }

  // Compute the partial derivative term
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a verbose report
  void PrintReport(ostream &sout);

  // Print a short name
  string GetShortName() { return string("MINANG"); }
  
private:
  // A structure that holds derivative data
  
  double UnifiedComputeEnergy(SolutionData *, bool); 
  StatisticsAccumulator sMinAngle, sPenalty;
};
*/


class MedialAnglesPenaltyTerm : public MedialIntegrationEnergyTerm
{
public:
  // Initialize the term
  MedialAnglesPenaltyTerm(GenericMedialModel *model);

  // Compute the energy
  double ComputeEnergy(SolutionData *data);

  // Print a verbose report
  void PrintReport(ostream &sout);

  // Compute the partial derivative term
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a short name
  string GetShortName() { return string("MEDANG"); }
  
private:
  double xTotalPenalty, xMaxPenalty;
};

/*
class CrestLaplacianEnergyTerm : public NumericalGradientEnergyTerm
{
public:
  double ComputeEnergy(SolutionData *data);
  void PrintReport(ostream &sout);
private:
  double xMaxLaplacian, xAvgLaplacian, xTotalPenalty;
  size_t nCrestAtoms, nBadSites;
  double PenaltyFunction(double x, double a, double b)
    { return exp( a * x - b ); }
};
*/

class AtomBadnessTerm : public EnergyTerm
{
public:
  double ComputeEnergy(SolutionData *data);
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Get the number of bad atoms
  size_t GetNumberOfBadAtoms() 
    { return nBadAtoms; }

  // Print a short name
  string GetShortName() { return string("ATMBAD"); }
  
private:
  vnl_vector<double> xPenalty;
  double xMinBadness, xAvgBadness, xTotalPenalty;
  size_t nBadAtoms, nAtoms;
};

/**
 * Term that penalizes |GradR| <> 1 on the boundary
 */
class BoundaryGradRPenaltyTerm : public EnergyTerm
{
public:
  // Compute the penalty
  double ComputeEnergy(SolutionData *data);

  // Describe the terms of the penalty
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a short name
  string GetShortName() { return string("GRAD-R"); }

private:
  // Accumulators for display and statistics calculation
  StatisticsAccumulator saGradR, saPenalty;

  // Scaling factor used to derive the penalty (why?)
  const static double xScale;
};

/**
 * Mesh regularization penalty term. This term computes a penalty
 * based on how much a mesh deviates from a surface parameterized
 * by Fourier harmonics. The number of harmonics in U and V is a 
 * parameter. For this to work, the u and v in each atom must be
 * specified properly
 */
class MeshRegularizationPenaltyTerm : public EnergyTerm
{
public:
  // Constructor, takes the initial model and number of components
  MeshRegularizationPenaltyTerm(
    GenericMedialModel *model, size_t nu, size_t nv);

  // Compute the penalty
  double ComputeEnergy(SolutionData *data);

  // Describe the terms of the penalty
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a short name
  string GetShortName() { return string("MSHREG"); }

private:
  // Accumulators for display and statistics calculation
  StatisticsAccumulator saPenalty, saDelta[3];

  // The matrix Q, such that the prior term has the form |x|^2 - |Qx|^2
  // matrix for a specific basis function decomposition
  vnl_matrix<double> Q, Qt, Z;

  // The vector b = 2(x - Q'Qx), useful to speed up the gradient 
  // computation
  vnl_vector<double> x[3], b[3];
};


/**
 * Term that penalizes excessive curvature of the medial axis. This
 * penalty has the form r^2 * (kappa1^2 + kappa2^2). Thus it penalizes
 * situations where one of the radii of curvature is excessively greater
 * than the radius of the medial model. Of course the model can avoid 
 * this penalty by pushing the radius to zero, so this penalty should
 * always be used in conjunction with a radius penalty.
 */
class MedialCurvaturePenalty : public EnergyTerm
{
public:
  // Compute the penalty
  double ComputeEnergy(SolutionData *data);

  // Describe the terms of the penalty
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a short name
  string GetShortName() { return string("MEDCRV"); }

private:
  // Accumulators for display and statistics calculation
  StatisticsAccumulator 
    saMeanCurv, saGaussCurv, saSumSqKappa, saRad, saFeature, saPenalty;

  // Parameters of the penalty, which is of the form pow(f/scale, power)
  const static double xScale, xPower;
};


/**
 * Term that penalizes excessive curvature of the medial axis. This
 * penalty has the form r^2 * (kappa1^2 + kappa2^2). Thus it penalizes
 * situations where one of the radii of curvature is excessively greater
 * than the radius of the medial model. Of course the model can avoid 
 * this penalty by pushing the radius to zero, so this penalty should
 * always be used in conjunction with a radius penalty.
 *
class BoundaryCurvaturePenalty : public EnergyTerm
{
public:
  // Constructr
  BoundaryCurvaturePenalty(GenericMedialModel *model);

  // Compute the penalty
  double ComputeEnergy(SolutionData *data);

  // Describe the terms of the penalty
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a short name
  string GetShortName() { return string("BNDCRV"); }

  // Get the curvature computed internally
  SMLVec3d GetCurvatureVector(size_t i)
    { return xMeanCurvVec[i]; }

private:
  // Temps
  GenericMedialModel *model;
  size_t nBnd;
  std::vector<SMLVec3d> xMeanCurvVec, dMeanCurvVec;
  double xIntegralSqrMeanCrv, xCrestCurvatureTerm, xPenalty;

  StatisticsAccumulator saCurv, saDenom;
};        */

class BoundaryCurvaturePenalty : public EnergyTerm
{
public:
  // Constructr
  BoundaryCurvaturePenalty(GenericMedialModel *model);

  // Compute the penalty
  double ComputeEnergy(SolutionData *data);

  // Describe the terms of the penalty
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a short name
  string GetShortName() { return string("BNDCRV"); }

private:
  // Temps
  GenericMedialModel *model;
  StatisticsAccumulator saSqrEdgeLen;
};

/*
class BoundaryCurvaturePenalty : public EnergyTerm
{
public:
  // Constructr
  BoundaryCurvaturePenalty(GenericMedialModel *model);

  // Compute the penalty
  double ComputeEnergy(SolutionData *data);

  // Describe the terms of the penalty
  void PrintReport(ostream &sout);

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

private:
  // Accumulators for display and statistics calculation
  StatisticsAccumulator 
    saMeanCurv, saGaussCurv, saSumSqKappa, saRad, saFeature, saPenalty;

  // Temps
  GenericMedialModel *model;
  size_t n;
  vnl_vector<double> MC, GC, dMC, dGC;

  // Parameters of the penalty, which is of the form pow(f/scale, power)
  const static double xScale, xPower;
};
*/

/**
 * This term is used to fit phi to an existing radius field. The fitting
 * is least-squares in phi, which may not be a great idea, but if we make it
 * least-squares in R, there may be unforseen problems with stability at R=0
 */
class DistanceToRadiusFieldEnergyTerm : public MedialIntegrationEnergyTerm
{
public:
  /** Constructor, takes the target radius field */
  DistanceToRadiusFieldEnergyTerm(GenericMedialModel *model, double *radius);
  
  double ComputeEnergy(SolutionData *data);
  
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  void PrintReport(ostream &sout);

  // Print a short name
  string GetShortName() { return string("DSTRAD"); }

private:

  typedef vnl_vector<double> Vec;
  Vec xTargetPhi, xLastDelta;
  double xTotalDiff, xMaxDiff, xMinDiff, xTotalMatch, xTotalArea;
};


/**
 * This term is used to fit the medial surface to an existing set of points 
 * given at u and v coordinates. This should be used in conjunction with a
 * regularization prior
 */
class DistanceToPointSetEnergyTerm : public MedialIntegrationEnergyTerm
{
public:
  /** Constructor, takes the target XYZ field */
  DistanceToPointSetEnergyTerm(
    GenericMedialModel *model, double *x, double *y, double *z);
  
  double ComputeEnergy(SolutionData *data);
  
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  void PrintReport(ostream &sout);

  // Print a short name
  string GetShortName() { return string("DSTPNT"); }

private:

  // Target points
  vector<SMLVec3d> target;

  // Match components
  double xTotalDist,  xTotalMatch, xTotalArea;
};



/**
 * Penalty for small values of the radius (phi / scale)^-8
 *
class RadiusPenaltyTerm : public EnergyTerm
{
public:
  RadiusPenaltyTerm(double xScale)
    { this->xScale = xScale; }

  double ComputeEnergy(SolutionData *data);
  
  void PrintReport(ostream &sout);

  // Print a short name
  string GetShortName() { return string("RADIUS"); }

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);
private:
  double xMinR2, xTotalPenalty, xScale;
};
*/

/**
 * A very simple radius penalty in the form
 *     { r > r_max:  (r - r_max)^2
 * p = { r < r_min:  (r - r_min)^2
 *     { 0        :  otherwise
 */
class RadiusPenaltyTerm : public EnergyTerm
{
public:
  RadiusPenaltyTerm(
    double bndLower = 0., double bndUpper = 1.e100, 
    double sclLower = 1., double sclUpper = 1.)
    {
    rMin = bndLower;
    rMax = bndUpper;
    sLower = sclLower;
    sUpper = sclUpper;
    }

  double ComputeEnergy(SolutionData *data);
  
  void PrintReport(ostream &sout);

  // Print a short name
  string GetShortName() { return string("RADIUS"); }

  // Compute the partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Set the parameters
  void SetParameters(Registry &R);


private:
  double rMin, rMax, sLower, sUpper;
  StatisticsAccumulator saRadius, saPenLow, saPenHigh;
};

class MedialBendingEnergyTerm : public MedialIntegrationEnergyTerm
{
public:
  /** Initialize with a sample solution */
  MedialBendingEnergyTerm(GenericMedialModel *model);
  
  /** Initialize the term with the template info */
  double ComputeEnergy(SolutionData *data);

  /** Compute partial derivative */
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);
  
  void PrintReport(ostream &sout);

  // Print a short name
  string GetShortName() { return string("MED-BE"); }
private:
  double xMaxBending, xTotalBending, xMeanBending;
};

class MedialRegularityTerm : public MedialIntegrationEnergyTerm
{
public:
  /** Initialize with a sample solution */
  MedialRegularityTerm(GenericMedialModel *model);
  
  /** Initialize the term with the template info */
  double ComputeEnergy(SolutionData *data);

  // Compute partial derivative
  double ComputePartialDerivative(
    SolutionData *S, PartialDerivativeSolutionData *dS);

  // Print a short name
  string GetShortName() { return string("MEDREG"); }
  
  void PrintReport(ostream &sout);
private:
  double xGradMagIntegral, xGradMagMaximum;
};

class MedialOptimizationProblem 
{
public:
  /** Initialize the problem */
  MedialOptimizationProblem(
    GenericMedialModel *xMedialModel, CoefficientMapping *xCoeff);

  ~MedialOptimizationProblem();

  /** Add an image match term with a perscribed weight */
  void AddEnergyTerm(EnergyTerm *term, double xWeight);

  /** Evaluate the function */
  double Evaluate(double *X);
  
  /** Compute the function and the gradient (comes from optima library) */
  double ComputeGradient(double *X, double *XGrad);

  /** Test gradient computation using central difference approximation */
  double TestGradientComputation(double *X, double eps);

  /** Compute the function and the partial derivative in one given direction */
  // double ComputePartialDerivative(double *xEvalPoint, double *xDirection, double &df);

  /** Print a detailed report */
  void PrintReport(ostream &sout);

  /** Get the medial model */
  GenericMedialModel *GetMedialModel()
    { return xMedialModel; }

  /** Get the list of energy terms */
  vector<EnergyTerm *> &GetEnergyTerms()
    { return xTerms; }

  /** Be quiet! */
  void QuietOn() 
    { flagQuiet = true; }

  void QuietOff() 
    { flagQuiet = false; }

  void DumpGradientMesh();

private:
  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  // Precision at which the problem is solved
  static const double xPrecision;

  // The epsilon used to compute finite difference derivatives
  static const double xEpsilon;

  // The number of coefficients
  size_t nCoeff;

  // Value of the last match
  double xLastSolutionValue;

  // Value of the solution for each of the energy terms
  Vec xLastTermValues, xLastGradEvalTermValues;

  // The value of the coefficients before optimization began
  Vec xInitialCoefficients;

  // The medial model being optimized 
  GenericMedialModel *xMedialModel;

  // The mapping from free parameters to the model's coefficients (e.g.,
  // affine transform)
  CoefficientMapping *xCoeff;

  // The weights of the different energy terms in the optimization
  vector<double> xWeights;

  // The pointers to the different energy terms
  vector<EnergyTerm *> xTerms;

  // Timers used to report code speeds
  vector<CodeTimer> xTimers, xGradTimers;
  CodeTimer xSolveTimer, xSolveGradTimer, xWeightsTimer, xWeightsGradTimer;

  // Whether the gradient is available
  bool flagLastEvalAvailable, flagPhiGuessAvailable;
  
  // Last place where the function was evaluated
  vnl_vector<double> xLastEvalPoint;

  // Last place where the gradient was evaluated and its value there
  vnl_vector<double> xLastGradPoint, xLastGradient, xLastGradHint;
  bool flagGradientComputed, flagQuiet;

  // Array of gradient vectors for each optimization term
  std::vector<vnl_vector<double>> xLastGradientPerTerm;

  // The solution at the last evaluation point
  vnl_matrix<double> xLastPhiField;

  // vnl_vector<double> xLastGradEvalPoint, xLastGrad;

  // The phi / dPhi fields at the last evaluation point
  // vnl_matrix<double> xLastPhiField, xLastPhiDerivField; 

  // This method solves the MedialPDE, potentially using the last gradient
  // evaluation as the guess
  bool SolvePDE(double *xEvalPoint);

  // Legacy central difference solver
  // void ComputeCentralDifferenceGradientPhi(double *x);

  // The representation of the variational basis
  Mat xBasis;

  // The array of derivative atoms (for each variational derivative)
  MedialAtom *dAtoms;

  // Statistics for learning scaling params
  vnl_vector<double> xGradSum, xGradSumSqr;

  // Data associated with the solution and its derivative
  SolutionData *S; 
  PartialDerivativeSolutionData *dS;

  size_t nGradCalls, nEvalCalls;
};

#endif
