#ifndef __ScriptInterface_h_
#define __ScriptInterface_h_

#include <string>
#include <vector>
#include <stdexcept>
#include <smlmath.h>

class GenericMedialModel;
class CartesianMedialModel;
class IBasisRepresentation2D;
class MedialOptimizationProblem;
class IMedialCoefficientMask;
class PrincipalComponents;
class PCACoefficientMask;
class Registry;
template <typename TPixel> class ITKImageWrapper;

void TestTriangleAreaPartialDerivative();
void TestTetrahedronVolumePartialDerivative();
void TestFDStuff();

class CoefficientMapping;
class OptimizationParameters;

namespace medialpde {

class BinaryImage;
class FloatImage;

/** This is a high level representation of a 3D non-binary image */
class FloatImage
{
public:
  typedef ITKImageWrapper<float> WrapperType;

  FloatImage();
  virtual ~FloatImage();

  // Load the image from a file
  void LoadFromFile(const char *file);

  // Load the image and the gradient from a path
  void LoadFromPath(const char *filebase, const char *ext);

  // Save the image to a file
  void SaveToFile(const char *file);

  // Load the image and the gradient from a path
  void SaveToPath(const char *filebase, const char *ext);

  // Load the image from a file
  void LoadGradientFromFile(unsigned int iComponent, const char *file);

  // Save the image to a file
  void SaveGradientToFile(unsigned int iComponent, const char *file);

  // Compute the image by processing a binary image
  void SetToGradientMagnitude(BinaryImage *imgSource, double xSigma);

  // Compute the image by blurring a binary image. The surface is represented
  // by the 0.0 level set, 1.0 is inside, -1.0 is outside of the object. The
  // optional error scale parameter applies the erf() function to the image,
  // to bring the values nearer to -1.0 or 1.0. This is useful for volume
  // overlap measures
  void SetToBlurredBinary(BinaryImage *imgSource, double xSigma, double xErrorScale = 0.0);

  // Compute the volume of the interior of the image. The image must represent
  // an follows: internal points near 1.0, external near -1.0, interface 0.0
  double ComputeObjectVolume();

  // Compute the integral of positive voxels in the image. This is the upper
  // bound for the probabilistic image match function that we use. In other
  // words, if the cm-rep model and the object in an image are perfectly
  // matched, the integral of the image inside the model will equal the output
  // of this function
  virtual double IntegratePositiveVoxels();

  // Get the voxel corresponding to a point. Return values set to -1 if point is outside
  virtual void GetVoxelIndex(const SMLVec3d &x, int &vx, int &vy, int &vz);

  // Interpolate the image at a given position
  virtual float Interpolate(const SMLVec3d &x);

  // Interpolate the image at a given position
  virtual float InterpolateNearestNeighbor(const SMLVec3d &x);

  // Interpolate the image gradient
  virtual void InterpolateImageGradient(const SMLVec3d &x, SMLVec3d &g);

  // Check whether the gradient information is available
  bool IsGradientAvailable();

  // Set the default 'outside' value (value of pixels outside of the image)
  void SetOutsideValue(float xValue);

  // Set the origin of the image, since some readers (Analyze) don't set it correctly
  // Do this after loading the image!
  void SetImageOrigin(double ox, double oy, double oz);

  // Get internal image
  WrapperType *GetInternalImage()
    { return xImage; } 

private:
  // Internally stored image
  WrapperType *xImage;

  // Optional gradient images - used to compute image match with gradient
  WrapperType *xGradient[3];

  // The outside value
  float xOutsideValue;

  // Allow the medial PDE to access this class
  friend class MedialPDE;
  friend class CartesianMPDE;
  friend class MedialPCA;
};

/** This is a binary image representation. */
class BinaryImage
{
public:
  BinaryImage();
  virtual ~BinaryImage();

  // Load the image from a file
  void LoadFromFile(const char *file);

private:
  // Internally stored image
  typedef ITKImageWrapper<unsigned char> WrapperType;
  WrapperType *xImage;

  // MedialPDE has access to private members
  friend class MedialPDE;
  friend class FloatImage;
};

/** This is the highest level class for working with the PDE application */
class MedialPDE
{
public:

  /** Public constructor. Reads the MPDE from a model file */
  MedialPDE(const char *file);

  /** Virtual destructor. Child must delete the medial model. */
  virtual ~MedialPDE();

  /**
   * Set the model. This method is virtual and child methods should check
   * that the model's type is correct
   */
  virtual void SetMedialModel(GenericMedialModel *model);

  /** Save to a parameter file */
  void SaveToParameterFile(const char *file);

  /** Load from a parameter file */
  void LoadFromParameterFile(const char *file);

  /** Compute the match between the model and a floating point image */
  double ComputeImageMatch(FloatImage *image);

  /** Compute the boundary jacobian penalty */
  double ComputeBoundaryJacobianPenalty(bool verbose = false);

  /** Fit the model to the binary image by matching moments of inertia */
  void MatchImageByMoments(FloatImage *image, unsigned int nCuts);

  /**
   * Perform an optimization step. Optimization parameters must be specified
   * in a parameter file. They have the following keys:
   *
   * Optimizer        String (ConjGrad, EvolStrat, Gradient)
   * Mapping          String (Identity, Affine, CoarseFine, PCA)
   * ImageMatch       String (VolumeOverlap, BoundaryIntegral)
   * TermWeights      Folder; each entry is a weight for a specific term
   *   BoundaryJacobian
   *   AtomBadness
   *   MedialRegularity
   *   RadiusPenalty
   *
   * Finally, the coarse-to-fine mapping can be fine-tuned. At the basic
   * level, the user can select a CTF level, which is a number between 1 and N
   * where N is the number of levels and 1 is the coarsest level. The finest
   * level can also be denoted as -1, down to N. However, for additional
   * fine-tuning, model-specific parameters can be provided.
   *
   * The idea is to eventually have a user interface for designing these
   * parameter files.
   */
  void RunOptimization(
    FloatImage *image, 
    size_t nSteps, 
    Registry &reg,
    FloatImage *imggray = NULL,
    bool do_deriv_test = false);
  
  // Same as above
  void RunOptimization(
    FloatImage *image, size_t nSteps, 
    const char *paramfile, const char *folderName = NULL);

  /** 
   * Configuration step for running optimization. Creates a mapping
   * object based on the passed in optimization parameters
   */
  CoefficientMapping *GenerateCoefficientMapping(
    OptimizationParameters &p);

  /** 
   * Configuration step for running optimization. Adds all energy terms
   * to the optimization problem with appropriate weights.
   */
  void ConfigureEnergyTerms(
    MedialOptimizationProblem &xProblem,                                     
    OptimizationParameters &p,
    FloatImage *image, 
    FloatImage *imgGray);

  /** Save the model as a BYU mesh */
  void SaveBYUMesh(const char *file);

  /** Save to a VTK file */
  void SaveVTKMesh(const char *fileMedial, const char *fileBoundary);

  /** Compute the radius function after surface/pho update */
  void Solve();

  /** Set the optimizer dump path. As the optimizer runs, it will periodically
   * dump meshes into path.number.med.vtk and path.number.bnd.vtk */
  void EnableMeshDump(const char *path, double xImprovement = 0.01)
    { strDumpPath = path; xMeshDumpImprovementPercentage = xImprovement; }

  void DisableMeshDump()
    { xMeshDumpImprovementPercentage = 0.0; }

  /** Enable the optimizer dump */
  void EnableOptimizerDump(const char *file)
    { strOptimizerDump = file; }

  void DisableOptimizerDump()
    { strOptimizerDump = ""; }

  /** Should not be here! */
  GenericMedialModel *GetMedialModel() { return xMedialModel; }

  /** Assign an intensity image to the cm-rep */
  void SetIntensityImage(FloatImage *imgIntensity);

  /** Get an intensity image from the cm-rep */
  void GetIntensityImage(FloatImage *imgIntensity);

  /** Set the PCA data for PCA-based optimization */
  void SetPCAMatrix(size_t ncu, size_t ncv, const char *fnMatrix);

  /** Store coordinates of interior (and exterior) points in a file */
  void SampleInterior(const char *file, double xStep, double xStart, double xEnd,
                      FloatImage *fim = NULL);


protected:

  /** Constructor. Takes the medial model from child class. */
  MedialPDE();

  // The solver
  GenericMedialModel *xMedialModel;

  // A file where the mesh info is dumped
  std::string strDumpPath, strOptimizerDump;
  double xMeshDumpImprovementPercentage;

  // An image containing grayscale intensities for this m-rep
  FloatImage imgIntensity;

  // Whether the intensities have been loaded
  bool flagIntensityPresent;

  // Fit the model to a radius field, used for fitting discrete m-reps, etc
  void FitPDEModelToRadius(double *rfield);
  void FitPDEModelToPointData(double *x, double *y, double *z);

  void ExportIterationToVTK(unsigned int iIter);
  void ConjugateGradientOptimization(MedialOptimizationProblem *xProblem,
    vnl_vector<double> &xSolution, unsigned int nSteps, double xStep);
  void ConjugateGradientOptimizationTOMS(MedialOptimizationProblem *xProblem,
    vnl_vector<double> &xSolution, unsigned int nSteps, double xStep);

  // Friend functions
  friend class MedialPCA;

  // PCA matrix (used to create PCA coefficient masks)
  vnl_matrix<double> mPCAMatrix;

  // Number of coefficients used in the PCA matrix
  size_t ncuPCA, ncvPCA, nPCAModes;
};

/**
 * This is a cm-rep on the Cartesian grid (basically unit box). It uses the
 * Fourier basis (cosine basis) to define the skelton and the rho function.
 */
class CartesianMPDE : public MedialPDE
{
public:

  /**
   * Constructor. Takes the number of basis functions in u and v and the
   * specification of the sampling grid
   */
  CartesianMPDE(
    unsigned int nBasesU, unsigned int nBasesV, unsigned int xResU, unsigned int xResV,
    double xFineScale = 0.0, unsigned int xFineU = 0, unsigned int xFineV = 0);

  /**
   * File-IO constructor
   */
  CartesianMPDE(const char *file);

  /** Access the medial surface */
  IBasisRepresentation2D *GetMedialSurface();

  /**
   * Set the medial model. This will ensure that the model is of the correct
   * type, and throw a ModelIOException if it is not.
   */
  virtual void SetMedialModel(GenericMedialModel *model);

  /** Set the number of Fourier coefficients */
  void SetNumberOfCoefficients(unsigned int m, unsigned int n);

  /** Set the size of the evaluation grid */
  void SetGridSize(size_t nu, size_t nv, bool useHint = false);

  /** Set the size of the evaluation grid with special sampling along edges */
  void SetGridSize(size_t nu, size_t nv,
    size_t eu, size_t ev, double eFactor = 0.5,
    bool useHint = false);

  /**
   * This method takes an image and an m-rep and samples the image using the
   * m-reps coordinate system, with linear interpolation. The result is stored in
   * a separate image.
   */
  void SampleImage(FloatImage *imgInput, FloatImage *imgOutput, size_t zSamples);

  /**
   * This method is the opposite of SampleImage. Given an image in the medial coordinate
   * system (it must match the sampling rate of the cm-rep), this with map this image back
   * into the ambient space
   */
  void SampleReferenceFrameImage(FloatImage *imgInput, FloatImage *imgOutput, size_t zSamples);

  /** Initialize the surface using a discrete m-rep. If the initial rho is set
   * to zero the method will perform an optimization problem, where we look for
   * the rho that gives the closest approximation of the radius function in the
   * m-rep */
  void LoadFromDiscreteMRep(const char *file, double xInitRho = 0);

  /** Import the model from a point list file. This file contains has the
   * structure [u v xi x y z rho r] */
  void ImportFromPointFile(const char *file, double xConstRho,
    bool flagFixRegularity, bool flagFitToRadius);

  /** Some default initialization */
  void GenerateSampleModel();

private:

  // Structure to represent importable data
  struct ImportAtom {
    double x, y, z;
    double u, v, xi, rho, r;
  };

  // Method to import a model from a set of sample points
  void ImportFromPointData(
    const std::vector<ImportAtom> &atoms,
    bool flagFixRegularity,
    bool flagFitToRadius);

  // The cartesian medial model
  CartesianMedialModel *xCartesianMedialModel;
};


/**
 * This is a cm-rep on Subdivision surfaces
 */
class SubdivisionMPDE : public MedialPDE
{
public:

  /**
   * File-IO constructor. Reads the surface from file.
   */
  SubdivisionMPDE(const char *file) : MedialPDE()
    {
    this->LoadFromParameterFile(file);
    }

  /**
   * Subdivide the cm-rep. The subdivision surface cm-rep is defined by a
   * coefficient mesh and an atom mesh, which are typically separated by a
   * level of atoms. This method will subdivide either of the meshes.
   */
  void SubdivideMeshes(size_t iCoeffSub, size_t iAtomSub);

  /**
   * Remesh coefficient and atom-level meshes
   */
  void Remesh();

  /**
   * Convert a Brute Force cm-rep to a PDE one. The radius function is used
   * as the initial guess for the PDE solution, and the LBO of the radius 
   * function from the brute force model is used to initialize rho.
   */
  void BruteForceToPDE();

private:

};

class MedialPCA
{
public:
  MedialPCA();
  ~MedialPCA();

  // Add a sample to the PCA
  void AddSample(MedialPDE *pde);

  // Compute the PCA
  void ComputePCA(MedialPDE *mpde);

  // Move current sample location to the mean
  void SetFSLocationToMean();

  // Move along a given mode a certain number of S.D.
  void SetFSLocation(unsigned int iMode, double xSigma);

  // Generate a sample at the current location
  void GetShapeAtFSLocation(MedialPDE *target);

  // Perform leave-one-out analysis using PCA
  // void LeaveOneOutAnalysis();

  // Export the PCA data for use in Mathematica, etc
  void ExportShapeMatrix(const char *filename);

private:

  typedef vnl_vector<double> Vec;
  typedef vnl_matrix<double> Mat;

  std::vector< FloatImage* > xAppearance;
  std::vector< Vec > xModelCoefficients;
  vnl_vector<double> xPCALocation;
  vnl_vector<double> xAppearancePCALocation;

  vnl_matrix<double> xDataShape, xDataAppearance;

  // Hint arrays for all the subjects
  std::vector<Vec> xHints;

  PrincipalComponents *xPCA, *xAppearancePCA;
};

} // Namespace medial PDE!

#endif
