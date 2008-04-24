#include "ScriptInterface.h"
#include "CartesianMedialModel.h"
#include "BasisFunctions2D.h"
#include "ITKImageWrapper.h"
#include "OptimizationTerms.h"
#include "Procrustes.h"
#include "PrincipalComponents.h"
#include "Registry.h"
#include "MedialAtomGrid.h"
#include "MedialModelIO.h"
#include "MedialException.h"
#include "PDESubdivisionMedialModel.h"
#include "BruteForceSubdivisionMedialModel.h"
#include "DiffeomorphicEnergyTerm.h"
#include "JacobianDistortionPenaltyTerm.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkBSplineDeformableTransform.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "TestSolver.h"

#include "vtkUnstructuredGrid.h"
#include "vtkPointLocator.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkCellLocator.h"


#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_random.h>


#include <iostream>
#include <algorithm>
#include <fstream>






// References to FORTRAN code
extern "C" {
  void deflt_(int &alg, int *iv, int &liv, int &lv, double *v);

  void sumsl_(
    int &n, double *d, double *x,
    void (*calcf)(int &, double *, int &, double &, int *, double *, void *),
    void (*calcg)(int &, double *, int &, double *, int *, double *, void *),
    int *iv, int &liv, int &lv, double *v,
    int *uiparm, double *urparm, void *ufparm);
}

// Reference to constants in Fortran code
const int mxiter_ = 18, mxfcal_ = 17, solprt_ = 22;

// VTK Export Methods
vtkUnstructuredGrid * ExportVolumeMeshToVTK(
  GenericMedialModel *xModel, size_t nSamples);

void ExportMedialMeshToVTK(
  GenericMedialModel *xModel, ITKImageWrapper<float> *xImage, const char *file);

void ExportBoundaryMeshToVTK(
  GenericMedialModel *xModel, ITKImageWrapper<float> *xImage, const char *file);

using namespace std;

void WriteMatrixFile(const vnl_matrix<double> &mat, const char *file)
{
  ofstream ofs(file, ios_base::out);
  ofs << "# MATRIX " << mat.rows() << " " << mat.columns() << endl;
  for(size_t r = 0; r < mat.rows(); r++)
    for(size_t c = 0; c < mat.columns(); c++)
      {
      ofs << mat[r][c];
      if(c == mat.columns() - 1) ofs << endl; else ofs << " ";
      }
  ofs.close();
}

/*
#include "gtb/graphics/indexed_triangle_set.hpp"
#include "gtb/graphics/triangle_mesh.h"
typedef gtb::tTriangleMesh<float> AFrontMesh;
extern int do_tri_mesh(int csubdiv, AFrontMesh &mesh);
*/


struct DiscreteAtom
{
  double u, v, x, y, z, r;
  unsigned int iu, iv;
};

/***************************************************************************
 * MedialPDE Code
 * ----------------
 *  This is the main application code
 **************************************************************************/
MedialPDE::MedialPDE()
{
  // Set default values of attributes
  xMeshDumpImprovementPercentage = 0.0;
  flagIntensityPresent = false;
  xMedialModel = NULL;
}

MedialPDE::MedialPDE(const char *file)
{
  // Set default values of attributes
  xMeshDumpImprovementPercentage = 0.0;
  flagIntensityPresent = false;
  xMedialModel = NULL;

  // Load the model
  LoadFromParameterFile(file);
}

MedialPDE::~MedialPDE()
{
  if(xMedialModel != NULL) delete xMedialModel;
}

void MedialPDE::SetMedialModel(GenericMedialModel *model)
{
  if(xMedialModel != model)
    {
    if(xMedialModel != NULL) delete xMedialModel;
    xMedialModel = model;
    }
}

void MedialPDE::SaveToParameterFile(const char *file)
{
  // Write the model to the file
  MedialModelIO::WriteModel(xMedialModel, file);
}

void MedialPDE::LoadFromParameterFile(const char *file)
{
  try
    {
    GenericMedialModel *model = MedialModelIO::ReadModel(file);
    this->SetMedialModel(model);
    }
  catch(ModelIOException &exc)
    {
    cerr << "Error reading model: " << exc.what() << endl;
    }
}

void MedialPDE::Solve()
{
  xMedialModel->ComputeAtoms();
}

/*
double ComputeImageMatchGradient(
  CartesianMedialModel *xMedialModel, PDEFourierWrapper *xModel, FourierSurface *xSurface,
  FloatImage *image, vnl_vector<double> &xGrad)
{
  unsigned int i, j, k;

  // Get the atom grid
  MedialAtomGrid *xAtomGrid = xMedialModel->GetAtomGrid();
  MedialAtom *xAtoms = xMedialModel->GetAtomArray();

  // Create a solution data object representing current solution
  SolutionData *S0 = new SolutionData(xMedialModel, true);

  // Solve the medial PDE for a set of finite difference offsets
  size_t nCoeff = xSurface->GetNumberOfRawCoefficients();
  double *xCoeff = xSurface->GetRawCoefficientArray();

  // Create an array of solution data objects
  SolutionData **SGrad = new SolutionData *[nCoeff];

  // Use current solution as the guess
  xMedialModel->SetSolutionAsInitialGuess();

  // Compute the time it takes to compute the derivative m-reps
  double tStart = clock();

  // Repeat for each coefficient
  MedialBoundaryPointIterator *itBnd = xAtomGrid->NewBoundaryPointIterator();
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    {
    // Increment the coefficient by epsilon
    double xOldValue = xCoeff[iCoeff];
    xCoeff[iCoeff] = xOldValue + xEpsilon;

    // Solve the PDE for the new value
    xMedialModel->Solve(xModel, 1.0e-13); cout << "." << flush;

    // Compute the solution data
    SGrad[iCoeff] = new SolutionData(xMedialModel, true);
    }

  // Restore the solution (should be immediate)
  xMedialModel->Solve(xModel); cout << endl;

  // See how much time elapsed for the medial computation
  double tMedial = clock() - tStart;

  // At this point, we have computed the solution at current X and at a set of
  // X + eps positions. We can use this to compute the gradient of each of the
  // terms involved in the optimization
  vnl_vector<double> xGradImage(nCoeff), xGradJacobian(nCoeff);
  double xSolution, xImageTerm, xJacobianTerm;

  // Compute the gradient of the image match term
  BoundaryImageMatchTerm termImage(image);
  xImageTerm = termImage.ComputeGradient(S0, SGrad, xEpsilon, xGradJacobian);

  // Compute the gradient of the jacobian penalty term
  BoundaryJacobianEnergyTerm termJacobian();
  xJacobianTerm = termJacobian.ComputeGradient(S0, SGrad, xEpsilon, xGradJacobian);

  // Finish timing
  double tGradient = clock() - tMedial;
  cout << "PDE Time: " << (tMedial / CLOCKS_PER_SEC) <<
   "; Gradient Time: " << (tGradient / CLOCKS_PER_SEC) << endl;

  // Clean up everything
  for(size_t iCoeff = 0; iCoeff < nCoeff; iCoeff++)
    delete SGrad[iCoeff];
  delete S0; delete SGrad;

  // Compute the gradient and the solution
  xGrad = xGradImage + 0.01 * xGradJacobian;
  return xImageTerm + 0.01 * xJacobianTerm;
}
*/

double MedialPDE::ComputeImageMatch(FloatImage *image)
{
  // Create a solution object (reference atoms in the solver)
  SolutionData S(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());

  // Create an image match term
  BoundaryImageMatchTerm termImage(xMedialModel, image);

  // Compute the energy
  double xMatch = termImage.ComputeEnergy(&S);

  // Print a report
  cout << "REPORT: " << endl;
  termImage.PrintReport(cout);

  return xMatch;
}

double MedialPDE::ComputeBoundaryJacobianPenalty(bool verbose)
{
  // Create a solution object
  SolutionData S(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());

  // Compute the penalty term
  BoundaryJacobianEnergyTerm bjet;
  double x = bjet.ComputeEnergy(&S);

  // Optionally, print a report
  if(verbose)
    {
    cout << "REPORT: " << endl;
    bjet.PrintReport(cout);
    }

  return x;
}

void GradientDescentOptimization(MedialOptimizationProblem *xProblem,
  vnl_vector<double> &xSolution, unsigned int nSteps, double xStep)
{
  // Create the initial solution
  size_t nCoeff = xSolution.size();
  vnl_vector<double> xGradient(nCoeff, 0.0);

  // Print report information
  // ofstream fdump("conjgrad.txt",ios_base::out);
  // fdump << "CONJUGATE GRADIENT OPTIMIZER DUMP" << endl;

  // Compute the 'last' value
  double xLastValue = xProblem->Evaluate(xSolution.data_block());

  // Iterate, moving in the gradient direction
  for(unsigned int p = 0; p < nSteps; p++)
    {
    // Compute the gradient and the image match
    double xMatch =
      xProblem->ComputeGradient( xSolution.data_block(), xGradient.data_block());

    // Take a step in the gradient direction
    vnl_vector<double> xTentativeSolution = xSolution - xStep * xGradient;

    // Compute the new solution
    double val = xProblem->Evaluate(xTentativeSolution.data_block());

    // Check if it got worse
    if(val > xLastValue)
      {
      while(val > xLastValue && xStep > 0.0001)
        {
        xStep = 0.5 * xStep;
        cout << "Step adjusted to " << xStep << endl;
        cout << "Gradient Magnitude " << xGradient.magnitude() << endl;
        xTentativeSolution = xSolution - xStep * xGradient;
        val = xProblem->Evaluate(xTentativeSolution.data_block());
        }
      if(val > xLastValue)
        break;
      }
    else
      {
      cout << "Step adjusted to " << xStep << endl;
      xStep *= 2;
      }

    // Store the last value
    xLastValue = val;
    xSolution -= xStep * xGradient;

    cout << "STEP " << p << "\t match " << xMatch << endl;


    // Report the current state
    // fdump << "STEP " << p << endl;
    // xProblem->PrintReport(fdump);

    // cout << "Step " << p << ", best value: " << xMatch << endl;
    }
}


void MedialPDE::ExportIterationToVTK(unsigned int iter)
{
  ostringstream oss;
  oss << strDumpPath << ".";
  oss.width(5);
  oss.fill('0');
  oss << iter;
  string fMedial = oss.str() + ".med.vtk";
  string fBoundary = oss.str() + ".bnd.vtk";

  ExportMedialMeshToVTK(xMedialModel, NULL, fMedial.c_str());
  ExportBoundaryMeshToVTK(xMedialModel, NULL, fBoundary.c_str());
}

void MedialPDE::SaveVTKMesh(const char *fMedial, const char *fBnd)
{
  ExportMedialMeshToVTK(xMedialModel, imgIntensity.xImage, fMedial);
  ExportBoundaryMeshToVTK(xMedialModel, imgIntensity.xImage, fBnd);
}

template<class TProblem>
class GSLProblemWrapper
{
public:
  double Evaluate(double *x)
    { return problem->Evaluate(x); }

  double ComputeGradient(double *x, double *df)
    {
    // if(flagFirstCall)
    //  {
      double f = problem->ComputeGradient(x, df);
      xLast = vnl_vector<double>(x, n);
      flagFirstCall = false;
      return f;
    //  }
    /* else
      {
      size_t i;
      double dfdv;
      vnl_vector<double> v = vnl_vector<double>(x, n) - xLast;
      cout << "DD !!" << endl;
      double f = problem->ComputePartialDerivative(x, v.data_block(), dfdv);

      for(i = 0; i < n; i++) df[i] = 0.0;
      for(size_t i = 0; i < n; i++)
        if(v[i] != 0.0)
          { df[i] = dfdv / v[i]; break; }

      return f;
      } */
    }

    void OnNextIteration()
      { flagFirstCall = true; }

    GSLProblemWrapper(TProblem *problem, size_t n)
      {
      this->problem = problem;
      this->n = n;
      flagFirstCall = true;
      }

private:
    TProblem *problem;
    size_t n;
    vnl_vector<double> xLast;
    bool flagFirstCall;
};



/** CALLBACK FUNCTIONS FOR TOMS611 ROUTINE */
void my_calcf(int &n, double *x, int &nf, double &f, int *dummy1, double *dummy2, void *info)
{
  static int it = 0;

  // Get a pointer to the object on which to evaluate
  MedialOptimizationProblem *mop = static_cast<MedialOptimizationProblem *>(info);

  // Check if x is ok
  vnl_vector<double> vx(x, n);
  if(!vx.is_finite())
    {
    nf = 0;
    return;
    }

  // Try to evaluate at this point, and if we fail return 0 in nf (request
  // smaller step from the solver, this is a nice feature of toms611)
  try { f = mop->Evaluate(x); }
  catch(MedialModelException &exc) 
    {
    cout << " [BAD SOLN (calcf) " << exc.what() << "]" << endl;
    nf = 0; 
    }

  // Try to flush stdout, otherwise output is garbled
  cout << flush;

  // Save a mesh
  // ExportMedialMeshToVTK(mop->GetMedialModel(), NULL, "iter.med.vtk");
  // ExportBoundaryMeshToVTK(mop->GetMedialModel(), NULL, "iter.bnd.vtk");

}

void my_calcg(int &n, double *x, int &nf, double *g, int *, double *, void *info)
{
  // Get a pointer to the object on which to evaluate
  MedialOptimizationProblem *mop = static_cast<MedialOptimizationProblem *>(info);

  // Try to evaluate at this point, and if we fail return 0 in nf (request
  // smaller step from the solver, this is a nice feature of toms611)
  try 
    { 
    mop->ComputeGradient(x, g); 
    }
  catch(MedialModelException &exc) 
    {
    cout << " [BAD SOLN (calcg) " << exc.what() << "]" << endl;
    nf = 0; 
    }
}

void MedialPDE::ConjugateGradientOptimizationTOMS(
  MedialOptimizationProblem *xProblem,
  vnl_vector<double> &xSolution, unsigned int nSteps, double xStep)
{
  int nCoeff = (int) xSolution.size();

  // For some reason there is a problem if we pass in VNL.data_blocks to the
  // fortran routines. So, instead, we will pass hand-created C arrays
  double *scaling = (double *) malloc(nCoeff * sizeof(double));
  double *x = (double *) malloc(nCoeff * sizeof(double));

  // Initialize the scaline and x arrays
  std::fill_n(scaling, nCoeff, 1.0);
  std::copy(xSolution.data_block(), xSolution.data_block() + nCoeff, x);

  // Specify the parameters to the sumsl_ routine
  int liv = 60, lv = 71+ nCoeff * (nCoeff+15) / 2;
  int *iv = (int *) malloc(liv * sizeof(int));
  double *v = (double *) malloc(lv * sizeof(double));

  std::fill_n(iv, liv, 0);
  std::fill_n(v, lv, 0.0);

  // Load the defaults
  int xAlg = 2;

  // Initialize the parameters of the method
  deflt_(xAlg, iv, liv, lv, v);
  iv[mxiter_ - 1] = nSteps;
  iv[mxfcal_ - 1] = 10 * nSteps;
  iv[19 - 1] = 0;
  iv[22 - 1] = 0;
  iv[24 - 1] = 0;

  // Execute the routine
  sumsl_(
    nCoeff, scaling, x,
    &my_calcf, &my_calcg,
    iv, liv, lv, v,
    NULL, NULL, (void *) xProblem);

  // Copy the solution value into xSolution
  xSolution.copy_in(x);

  // Evaluate the problem at the best solution
  xProblem->Evaluate(xSolution.data_block());
  xProblem->PrintReport(std::cout);

  // Free the data
  free(iv); free(v); free(scaling); free(x);
}


void MedialPDE::ConjugateGradientOptimization(
  MedialOptimizationProblem *xProblem,
  vnl_vector<double> &xSolution, unsigned int nSteps, double xStep)
{
  this->ConjugateGradientOptimizationTOMS(xProblem, xSolution, nSteps, xStep);
}

void EvolutionaryOptimization(
  MedialOptimizationProblem *xProblem,
  vnl_vector<double> &xSolution,
  unsigned int nSteps)
{
  assert(0);
}

/*

void MedialPDE
::RunOptimization(
  FloatImage *image, size_t nSteps,
  const char *paramfile, const char *folderName)
{
  Registry reg;
  reg.ReadFromFile(paramfile);
  if(folderName)
    {
    Registry fld = reg.Folder(folderName);
    this->RunOptimization(image,nSteps,fld);
    }
  else
    {
    this->RunOptimization(image,nSteps,reg);
    }
}

*/

CoefficientMapping *
MedialPDE::GenerateCoefficientMapping(OptimizationParameters &p)
{
  // Create the coefficient mapping
  CoefficientMapping *xMapping = NULL;

  // We might need a PCA as well
  PrincipalComponents *pca = NULL;

  // Create an appropriate optimization mapping
  if(p.xMapping == OptimizationParameters::AFFINE)
    {
    xMapping = new AffineTransformCoefficientMapping(xMedialModel);
    }
  else if(p.xMapping == OptimizationParameters::IDENTITY)
    {
    xMapping = new IdentityCoefficientMapping(xMedialModel->GetNumberOfCoefficients());
    }
  else if(p.xMapping == OptimizationParameters::PCA)
    {
    try
      {
      // Read principal components from the file
      vnl_matrix<double> pcaMatrix;
      ReadMatrixFile(pcaMatrix, p.xPCAFileName.c_str());

      // TODO: this must be deleted somewhere!!!
      pca = new PrincipalComponents(pcaMatrix);

      // Create the PCA/Affine optimizer
      xMapping = new PCAPlusAffineCoefficientMapping(xMedialModel, pca, p.nPCAModes);
      }
    catch(...)
      {
      throw ModelIOException("Unable to read PCA matrix file in optimization parameters");
      }
    }
  else if(p.xMapping == OptimizationParameters::COARSE_TO_FINE)
    {
    // Try to use the coarse-to-fine settings. If unavailable, use whole model
    if(p.xCTFSettings == NULL)
      throw ModelIOException("No coarse-to-fine settings specified!");

    // Create the coarse to fine mask and mapping
    const CoarseToFineMappingDescriptor *ctfDesc
      = xMedialModel->GetCoarseToFineMappingDescriptor();
    xMapping = new SubsetCoefficientMapping(ctfDesc->GetMask(p.xCTFSettings));
    }
  else if(p.xMapping == OptimizationParameters::RADIUS_SUBSET)
    {
    // Create a coarse to fine mask for the radius
    vnl_vector<size_t> mask = xMedialModel->GetRadialCoefficientMask();

    // Create a mapping
    xMapping = new SubsetCoefficientMapping(mask);
    }
  else if(p.xMapping == OptimizationParameters::POSITION_SUBSET)
    {
    // Create a coarse to fine mask for the radius
    vnl_vector<size_t> mask = xMedialModel->GetSpatialCoefficientMask();

    // Create a mapping
    xMapping = new SubsetCoefficientMapping(mask);
    }
  else if(p.xMapping == OptimizationParameters::REFLECTION)
    {
    try 
      {
      SubdivisionMedialModel *smm = 
        reinterpret_cast<SubdivisionMedialModel *>(xMedialModel);
      xMapping = 
        new ReflectionCoefficientMapping(smm, p.xReflectionPlane, p.xReflectionIntercept);
      }
    catch(...)
      {
      throw ModelIOException("Reflection coefficient mapping works only for subdivision models");
      }
    }
  else
    {
    throw ModelIOException("Unknown coefficient mapping specified");
    }
  
  return xMapping;
}

void MedialPDE::ConfigureEnergyTerms(
  MedialOptimizationProblem &xProblem,                                     
  OptimizationParameters &p,
  FloatImage *image, 
  FloatImage *imgGray)
{
  // Configure the image match term (TODO: these should be treated same as
  // penalty terms
  EnergyTerm *xTermImage = NULL;
  switch(p.xImageMatch)
    {
    case OptimizationParameters::VOLUME:
      xTermImage = new VolumeOverlapEnergyTerm(xMedialModel, image, 8);
      break;

    case OptimizationParameters::BOUNDARY:
      xTermImage = new BoundaryImageMatchTerm(xMedialModel, image);
      break;

    case OptimizationParameters::RADIUS_VALUES:
      {
      vnl_vector<double> xRadius(xMedialModel->GetNumberOfAtoms());
      for(size_t i = 0; i < xRadius.size(); i++)
        xRadius[i] = xMedialModel->GetAtomArray()[i].R;
      xTermImage = new DistanceToRadiusFieldEnergyTerm(xMedialModel, xRadius.data_block());
      }
      break;

    default:
      throw MedialModelException("Unknown energy term specified");
    }

  // Add the energy term to the problem
  xProblem.AddEnergyTerm(xTermImage, 1.0);

  // Add all the penalty terms
  for(size_t iTerm = 0; iTerm < OptimizationParameters::N_PENALTY_TERMS; iTerm++)
    {
    // Ignore all terms with zero weights
    OptimizationParameters::PenaltyTerm pt = 
      static_cast<OptimizationParameters::PenaltyTerm>(iTerm);

    if(p.xTermWeights[pt] == 0.0)
      continue;

    // Create the energy term    
    EnergyTerm *xTermPenalty;
    switch(pt)
      {
      case OptimizationParameters::BOUNDARY_JACOBIAN:
        xTermPenalty = new BoundaryJacobianEnergyTerm(); break;
      case OptimizationParameters::ATOM_BADNESS:
        xTermPenalty = new AtomBadnessTerm(); break;
      case OptimizationParameters::BOUNDARY_GRAD_R:
        xTermPenalty = new BoundaryGradRPenaltyTerm(); break;
      case OptimizationParameters::BND_JACOBIAN_DISTORTION:
        xTermPenalty = new BoundaryJacobianDistortionPenaltyTerm(xMedialModel); break;
      case OptimizationParameters::MED_JACOBIAN_DISTORTION:
        xTermPenalty = new MedialJacobianDistortionPenaltyTerm(xMedialModel); break;
      case OptimizationParameters::MEDIAL_REGULARITY:
        xTermPenalty = new MedialRegularityTerm(xMedialModel); break;
      case OptimizationParameters::MEDIAL_ANGLES:
        xTermPenalty = new MedialAnglesPenaltyTerm(xMedialModel); break;
      case OptimizationParameters::MEDIAL_CURVATURE:
        xTermPenalty = new MedialCurvaturePenalty(); break;
      case OptimizationParameters::BND_CURVATURE:
        xTermPenalty = new BoundaryCurvaturePenalty(xMedialModel); break;
      case OptimizationParameters::RADIUS:
        xTermPenalty = new RadiusPenaltyTerm(); break;
      case OptimizationParameters::DIFFEOMORPHIC:
        xTermPenalty = new DiffeomorphicEnergyTerm(xMedialModel); break;
      case OptimizationParameters::CROSS_CORRELATION:
        xTermPenalty = new CrossCorrelationImageMatchTerm(xMedialModel, imgGray); break;
      case OptimizationParameters::LOCAL_DISTANCE:
        xTermPenalty = new LocalDistanceDifferenceEnergyTerm(xMedialModel); break;
      default:
        throw MedialModelException("Unknown Penalty Term");
      }

    // Configure the energy term's parameters
    xTermPenalty->SetParameters(p.xTermParameters[pt]);

    // Add the penalty term to the optimization problem
    xProblem.AddEnergyTerm(xTermPenalty, p.xTermWeights[pt]); 
    }
}


void MedialPDE
::RunOptimization(
  FloatImage *image, size_t nSteps, Registry &folder, OptimizationFlags flags, FloatImage *imgGray)
{
  // Read the optimization parameters from the file
  OptimizationParameters p; p.ReadFromRegistry(folder);

  // Create the coefficient mapping
  CoefficientMapping *xMapping = GenerateCoefficientMapping(p);

  // Create the initial solution (all zeros)
  vnl_vector<double> xSolution(xMapping->GetNumberOfParameters(), 0.0);

  // Save the initial values of the coefficients
  vnl_vector<double> xInitialCoeff = xMedialModel->GetCoefficientArray();

  // Create the optimization problem
  MedialOptimizationProblem xProblem(xMedialModel, xMapping);

  // Set up debugging features
  if(flags.flagDumpGradientMesh)
    xProblem.DumpGradientMeshOn();

  // Add all terms to the optimization problem
  ConfigureEnergyTerms(xProblem, p, image, imgGray);

  // Initial solution report
  cout << "INITIAL SOLUTION REPORT: " << endl;
  xProblem.Evaluate(xSolution.data_block());
  xProblem.PrintReport(cout);

  // Test the term-by-term gradient computation
  if(flags.flagTestGradient)
    {
    printf("TESTING GRADIENT COMPUTATION\n");
    printf("%12s  %12s : %9s  %9s  %9s  %9s  %9s  P/F\n",
      "ENERGY TERM", "COEFFICIENT MAPPING",
      "ENERGY VAL", 
      "ABS ERROR", "REL ERROR", "T(AGRAD)", "T(NGRAD)");

    vector<EnergyTerm *> &terms = xProblem.GetEnergyTerms();
    for(size_t i = 0; i < terms.size(); i++)
      {
      MedialOptimizationProblem mop(xMedialModel, xMapping);
      mop.AddEnergyTerm(terms[i], 1.0);
      TestOptimizerGradientComputation(
        mop, *xMapping, xMedialModel, 0.0, 
        terms[i]->GetShortName().c_str(), "X");
      }
    }

  // Get the step size for optimization
  double xStepSize = folder["StepSize"][10.0];

  // At this point, split depending on the method
  if(p.xOptimizer == OptimizationParameters::CONJGRAD)
    ConjugateGradientOptimizationTOMS(&xProblem, xSolution, nSteps, xStepSize);
  else if(p.xOptimizer == OptimizationParameters::GRADIENT)
    GradientDescentOptimization(&xProblem, xSolution, nSteps, xStepSize);
  else if(p.xOptimizer == OptimizationParameters::EVOLUTION)
    EvolutionaryOptimization(&xProblem, xSolution, nSteps);
  else throw ModelIOException("Unknown optimization technique");

  // Test the gradient computation again
  // cout << "GRADIENT COMPUTATION ABS MAX ERROR (eps = 0.001): " << 
  //   xProblem.TestGradientComputation(xSolution.data_block(), 0.001);

  // After optimization, apply the best result
  xMedialModel->SetCoefficientArray(xMapping->Apply(xInitialCoeff, xSolution));
  xMedialModel->ComputeAtoms();

  // Delete the mapping
  delete xMapping;
  
  // Delete all the energy terms
  for(size_t i = 0; i < xProblem.GetEnergyTerms().size(); i++)
    delete xProblem.GetEnergyTerms()[i];
}


/******************************************************************************
 * PROBLEM: FIT PDE-LIKE problem to a given radius field
 *****************************************************************************/
void MedialPDE
::FitPDEModelToRadius(double *rfield)
{
  // The coefficient mapping selects only the rho-related components of the 
  // coefficient vector - the manifold does not deform
  SubsetCoefficientMapping xCoeffMapping(
    xMedialModel->GetRadialCoefficientMask());

  // Create the initial solution (all zeros)
  vnl_vector<double> xSolution(xCoeffMapping.GetNumberOfParameters(), 0.0);

  // Save the initial values of the coefficients
  vnl_vector<double> xInitialCoeff = xMedialModel->GetCoefficientArray();

  // Create an optimization problem (same problem as always)
  MedialOptimizationProblem xProblem(xMedialModel, &xCoeffMapping);

  // Create the match term based on the radius field
  DistanceToRadiusFieldEnergyTerm xTermRadius(xMedialModel, rfield);
  xProblem.AddEnergyTerm(&xTermRadius, 1.0);

  // Add a penalty for bad atoms
  AtomBadnessTerm xTermBadness;
  xProblem.AddEnergyTerm(&xTermBadness, 1.0);

  // Add a penalty for boundary jacobian
  BoundaryJacobianEnergyTerm xTermJacobian;
  xProblem.AddEnergyTerm(&xTermJacobian, 1.0);

  // Print the initial solution
  cout << "INITIAL SOLUTION REPORT: " << endl;
  xProblem.Evaluate(xSolution.data_block());
  xProblem.PrintReport(cout);

  // Run the optimization
  ConjugateGradientOptimizationTOMS(&xProblem, xSolution, 200, 0.001);

  // After optimization, apply the best result
  xMedialModel->SetCoefficientArray(xCoeffMapping.Apply(xInitialCoeff, xSolution));
  xMedialModel->ComputeAtoms();
}



/******************************************************************************
 * PROBLEM: Interpolate data under regularity constraints
 *****************************************************************************/
void MedialPDE
::FitPDEModelToPointData(double *x, double *y, double *z)
{
  // The coefficient mapping selects only the rho-related components of the 
  // coefficient vector - the manifold does not deform
  SubsetCoefficientMapping xCoeffMapping(
    xMedialModel->GetSpatialCoefficientMask());

  // Create the initial solution (all zeros)
  vnl_vector<double> xSolution(xCoeffMapping.GetNumberOfParameters(), 0.0);

  // Save the initial values of the coefficients
  vnl_vector<double> xInitialCoeff = xMedialModel->GetCoefficientArray();

  // Create an optimization problem (same problem as always)
  MedialOptimizationProblem xProblem(xMedialModel, &xCoeffMapping);

  // Create the match term based on the radius field
  DistanceToPointSetEnergyTerm xTermDistance(xMedialModel, x, y, z);
  xProblem.AddEnergyTerm(&xTermDistance, 1.0);

  // Add a penalty for bad atoms
  MedialRegularityTerm xTermRegularity(xMedialModel);
  xProblem.AddEnergyTerm(&xTermRegularity, 0.01);

  // Print the initial solution
  cout << "INITIAL SOLUTION REPORT: " << endl;
  xProblem.Evaluate(xSolution.data_block());
  xProblem.PrintReport(cout);

  // Run the optimization
  ConjugateGradientOptimizationTOMS(&xProblem, xSolution, 100, 0.1);

  // After optimization, apply the best result
  xMedialModel->SetCoefficientArray(xCoeffMapping.Apply(xInitialCoeff, xSolution));
  xMedialModel->ComputeAtoms();
}



class FirstMomentComputer : public EuclideanFunction
{
public:
  FirstMomentComputer(unsigned int i)
    { this->i = i; }

  double Evaluate(const SMLVec3d &X)
    { return X[i]; }

private:
  unsigned int i;
};

class SecondMomentComputer : public EuclideanFunction
{
public:
  SecondMomentComputer(unsigned int i, unsigned int j)
    { this->i = i; this->j = j;}

  double Evaluate(const SMLVec3d &X)
    { return X[i] * X[j]; }

private:
  unsigned int i, j;
};

template<typename TPixel>
class ImageVolumeMatch : public EuclideanFunction
{
public:
  ImageVolumeMatch(ITKImageWrapper<TPixel> *xWrapper)
    { this->xWrapper = xWrapper; }

  double Evaluate(const SMLVec3d &X)
    { return xWrapper->Interpolate(X[0], X[1], X[2], 0.0f); }

private:
  ITKImageWrapper<TPixel> *xWrapper;
};

/** Match an mrep to a floating point (level set) image using moments. */
void MedialPDE::MatchImageByMoments(FloatImage *image, unsigned int nCuts)
{
  // Typedefs
  typedef vnl_matrix_fixed<double,3,3> Mat;
  typedef SMLVec3d Vec;

  // Compute the mean and covariance matrix of the non-zero voxels in the
  // image.
  Mat xCov(0.0f); Vec xMean(0.0f); double xVolume = 0.0;
  size_t i, j, f;

  // Iterate over all voxels in the image
  typedef FloatImage::WrapperType::ImageType ImageType;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> IteratorType;
  ImageType::Pointer xImage = image->xImage->GetInternalImage();
  IteratorType it(xImage, xImage->GetBufferedRegion());
  for( ; !it.IsAtEnd(); ++it)
    {
    // Get the volume of this voxel
    if(it.Get() > 0.0)
      {
      // Get the spatial position of the point
      itk::ContinuousIndex<double, 3> iRaw(it.GetIndex());
      itk::Point<double, 3> ptPosition;
      xImage->TransformContinuousIndexToPhysicalPoint(iRaw, ptPosition);
      SMLVec3d xPosition(ptPosition.GetDataPointer());

      // Add to the mean and 'covariance'
      xVolume += 1.0;
      xMean += xPosition;
      xCov += outer_product(xPosition, xPosition);
      }
    }

  // Compute the actual volume of the image
  double xPhysicalVolume = xVolume *
    xImage->GetSpacing()[0] * xImage->GetSpacing()[1] * xImage->GetSpacing()[2];

  // Scale the mean and covariance by accumulated weights
  xMean = xMean / xVolume;
  xCov = (xCov - xVolume *  outer_product(xMean, xMean)) / xVolume;

  // Compute the mean and covariance of the image
  cout << "--- MATCHING BY MOMENTS ---" << endl;
  cout << "Image Volume: " << xPhysicalVolume << endl;
  cout << "Image Mean: " << xMean << endl;
  cout << "Image Covariance: " << endl << xCov << endl;

  // Now compute the same statistics for the m-rep (compute coordinates)
  Mat yCov(0.0); Vec yMean(0.0); double yVolume = 0.0;

  // Create a solution data object to speed up computations
  SolutionData S0(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());
  S0.ComputeIntegrationWeights();

  for(i = 0; i < 3; i++)
    {
    // Compute the first moment in i-th direction
    FirstMomentComputer fmc(i);
    VolumeIntegralEnergyTerm viet1(xMedialModel, &fmc, nCuts);
    yMean[i] = viet1.ComputeEnergy(&S0);
    yVolume = viet1.GetModelVolume();

    for(j = 0; j <= i; j++)
      {
      // Compute the second moment in i-th and j-th directions
      SecondMomentComputer smc(i,j);
      VolumeIntegralEnergyTerm viet2(xMedialModel, &smc, nCuts);
      yCov[i][j] = yCov[j][i] = viet2.ComputeEnergy(&S0);
      }
    }

  // Compute the actual mean and covariance
  yMean /= yVolume;
  yCov = (yCov - yVolume * outer_product(yMean, yMean)) / yVolume;

  cout << "Model Volume: " << yVolume << endl;
  cout << "Model Mean: " << yMean << endl;
  cout << "Model Covariance: " << endl << yCov << endl;

  // Now, compute the alignment that will line up the moments
  // Decompose xCov and yCov into eigenvalues
  vnl_vector<double> Dx(3), Dy(3);
  vnl_matrix<double> Vx(3,3), Vy(3,3);
  vnl_symmetric_eigensystem_compute(xCov, Vx, Dx);
  vnl_symmetric_eigensystem_compute(yCov, Vy, Dy);

  // Compute the scale factor
  double s = sqrt(dot_product(Dx,Dy) / dot_product(Dy,Dy));

  // The best set of coefficients and the associated match value
  size_t iBestSurface; double xBestMatch;

  // Get the numbers of coefficients
  size_t nCoeff = xMedialModel->GetNumberOfCoefficients();
  vnl_vector<double> xRotatedCoeff[8], xInitCoeff = xMedialModel->GetCoefficientArray();

  // Create an affine transform for the coefficients
  const AffineTransformDescriptor *affDesc = xMedialModel->GetAffineTransformDescriptor();

  // Create a volume match object
  VolumeOverlapEnergyTerm tVolumeMatch(xMedialModel, image, nCuts);

  // Compute the best image match over all possible flips of the eigenvalues
  // (there are 8 possible matches, including mirroring)
  for(f = 0; f < 8; f++)
    {
    // Set the flip/scale matrix
    vnl_matrix<double> F(3, 3, 0.0);
    // F(0,0) = (f & 1) ? -sqrt(Dx(0) / Dy(0)) : sqrt(Dx(0) / Dy(0));
    // F(1,1) = (f & 2) ? -sqrt(Dx(1) / Dy(1)) : sqrt(Dx(1) / Dy(1));
    // F(2,2) = (f & 4) ? -sqrt(Dx(2) / Dy(2)) : sqrt(Dx(2) / Dy(2));

    // Flips only, no scaling, to avoid breaking model
    F(0,0) = (f & 1) ? -1.0 : 1.0;
    F(1,1) = (f & 2) ? -1.0 : 1.0;
    F(2,2) = (f & 4) ? -1.0 : 1.0;

    // Compute the rotation+flip matrix
    vnl_matrix<double> R = Vx * F * Vy.transpose();

    // Rotate the surface by matrix R and shift to yMean
    xRotatedCoeff[f] =
      affDesc->ApplyAffineTransform( xInitCoeff, R, xMean - yMean, yMean);
    xMedialModel->SetCoefficientArray(xRotatedCoeff[f]);

    // Compute the boundary
    xMedialModel->ComputeAtoms();

    // Create a solution object and a volume match term
    SolutionData SRot(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());
    SRot.ComputeIntegrationWeights();
    double xMatch = tVolumeMatch.ComputeEnergy(&SRot);

    // Record the best match ever
    if(f == 0 || xMatch < xBestMatch)
      { xBestMatch = xMatch; iBestSurface = f; }

    // Report the match quality
    cout << "ROTATION " << f << " :" << endl;
    tVolumeMatch.PrintReport(cout);
    }

  // Use the best surface as the new surface
  xMedialModel->SetCoefficientArray(xRotatedCoeff[iBestSurface]);
  xMedialModel->ComputeAtoms();

  // Test the results
  SolutionData S1(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());
  S1.ComputeIntegrationWeights();
  for(i = 0; i < 3; i++)
    {
    // Compute the first moment in i-th direction
    FirstMomentComputer fmc(i);
    VolumeIntegralEnergyTerm viet1(xMedialModel, &fmc, nCuts);
    yMean[i] = viet1.ComputeEnergy(&S1); 
    yVolume = viet1.GetModelVolume();

    for(j = 0; j <= i; j++)
      {
      // Compute the second moment in i-th and j-th directions
      SecondMomentComputer smc(i,j);
      VolumeIntegralEnergyTerm viet2(xMedialModel, &smc, nCuts);
      yCov[i][j] = yCov[j][i] = viet2.ComputeEnergy(&S1);
      }
    }

  // Compute the actual mean and covariance
  yMean /= yVolume;
  yCov = (yCov - yVolume * outer_product(yMean, yMean)) / yVolume;

  cout << "Model Volume: " << yVolume << endl;
  cout << "Model Mean: " << yMean << endl;
  cout << "Model Covariance: " << endl << yCov << endl;
}

void MedialPDE::SaveBYUMesh(const char *file)
{
  // Export the model as a BYU mesh
  fstream fout(file, ios_base::out);

  // Get the atom grid
  MedialIterationContext *grid = xMedialModel->GetIterationContext();
  MedialAtom *atoms = xMedialModel->GetAtomArray();

  // Write the number of vertices, edges, parts
  unsigned int nVertices = grid->GetNumberOfBoundaryPoints();
  unsigned int nFaces = grid->GetNumberOfBoundaryTriangles();
  unsigned int nEdges = 2 * (nFaces + nVertices - 2);
  fout << 1 << "\t" << nVertices << "\t" << nFaces << "\t" << nEdges << "\t" << endl;

  // Write the number of faces in the first part
  fout << 1 << "\t" << nFaces << endl;

  // Write all the vertices
  for(MedialBoundaryPointIterator itv(grid); !itv.IsAtEnd(); ++itv)
    fout << GetBoundaryPoint(itv, atoms).X << endl;

  // Write all the faces
  for(MedialBoundaryTriangleIterator itt(grid); !itt.IsAtEnd(); ++itt)
    {
    int i0 = 1 + itt.GetBoundaryIndex(0);
    int i1 = 1 + itt.GetBoundaryIndex(1);
    int i2 = 1 + itt.GetBoundaryIndex(2);

    fout << i0 << "\t" << i1 << "\t" << -i2 << endl;
    }

  // Close the file
  fout.close();
}

// This method estimates the parameters of a transform to a set of point pairs
// using least square fitting
void LeastSquaresFit(
  itk::BSplineDeformableTransform<double, 3, 3> *t,
  size_t nx, SMLVec3d *x, SMLVec3d *y)
{
  // Get the number of parameters
  size_t np = t->GetNumberOfParameters(), p, i;

  // Store the original parameters
  typedef itk::BSplineDeformableTransform<double, 3, 3> TTransform;
  TTransform::ParametersType p0(t->GetNumberOfParameters());

  // Allocate three matrices
  vnl_matrix<double> Z1(np, nx), Z2(np, nx), Z3(np, nx);

  // Compute the right hand side vectors
  vnl_vector<double> b1(nx), b2(nx), b3(nx);


  for(p = 0; p < np; p++)
    {
    // Create a basis function
    p0.fill(0.0); p0[p] = 1.0;
    t->SetParameters(p0);

    // Compute values for all x
    for(i = 0; i < nx; i++)
      {
      TTransform::InputPointType xin;
      TTransform::OutputPointType xout;
      xin[0] = x[i][0]; xin[1] = x[i][1]; xin[2] = x[i][2];
      xout = t->TransformPoint(xin);
      Z1[p][i] = xout[0]; Z2[p][i] = xout[1]; Z3[p][i] = xout[2];

      b1[i] = y[i][0]; b2[i] = y[i][1]; b3[i] = y[i][2];
      }
    }

  // Print stuff
  cout << "Fitting BSpline parameters" << endl;
  cout << "  Coefficients    : " << np << endl;
  cout << "  Landmark Points : " << nx << endl;

  // Compute the square matrix A
  vnl_matrix<double> A =
    Z1 * Z1.transpose() + Z2 * Z2.transpose() + Z3 * Z3.transpose();
  vnl_vector<double> B = Z1 * b1 + Z2 * b2 + Z3 * b3;

  cout << "  Matrix A        : " << A.rows() << " x " << A.columns() << endl;
  cout << "  Vector B        : " << B.size() << endl;

  // Solve the system Ax = b (LU decomposition)
	vnl_qr<double> qr(A);
	vnl_vector<double> Y = qr.solve(B);

  // Set the parameters
  for(p = 0; p < np; p++)
    p0[p] = Y[p];

  t->SetParametersByValue(p0);

  // Double check that the mean squared difference is reasonable
  double distMax = 0.0;
  for(i = 0; i < nx; i++)
    {
    TTransform::InputPointType xin;
    TTransform::OutputPointType xout;
    xin[0] = x[i][0]; xin[1] = x[i][1]; xin[2] = x[i][2];
    xout = t->TransformPoint(xin);
    SMLVec3d z(xout[0], xout[1], xout[2]);
    double dist = (z - y[i]).two_norm();
    if(dist > distMax) distMax = dist;
    }

  cout << "  Max Sqr. Err   : " << sqrt(distMax) << endl;
}

// Update region by a point
void UpdateRegion(itk::ImageRegion<3> &R, const itk::Index<3> &idx, bool first)
{
  if(first)
    {
    R.SetIndex(idx);
    R.SetSize(0, 1); R.SetSize(1, 1); R.SetSize(2, 1);
    return;
    }

  if(R.IsInside(idx))
    { return; }

  // Update the index
  for(size_t i = 0; i < 3; i++)
    {
    if(R.GetIndex(i) < idx[i])
      R.SetIndex(i, idx[i]);
    if(R.GetSize(i) + R.GetIndex(i) <= idx[i])
      R.SetSize(i, 1 + idx[i] - R.GetIndex(i));
    }
}



void MedialPDE::SetIntensityImage(FloatImage *imgSource)
{
  // Copy the internal image
  imgIntensity.xImage->SetInternalImage(imgSource->xImage->GetInternalImage());

  // Set the flag
  flagIntensityPresent = true;
}

void MedialPDE::GetIntensityImage(FloatImage *imgTarget)
{
  imgTarget->xImage->SetInternalImage(imgIntensity.xImage->GetInternalImage());
}

void MedialPDE::SetPCAMatrix(size_t ncu, size_t ncv, const char *fname)
{
  // Store the dimensions
  ncuPCA = ncu;
  ncvPCA = ncv;

  // Read the matrix
  ReadMatrixFile(mPCAMatrix, fname);
  cout << "READ PCA MATRIX " << mPCAMatrix.rows()
    << " x " << mPCAMatrix.columns() << endl;
}

void MedialPDE::SampleInterior(
  const char *file,
  double xStep,
  double xStart,
  double xEnd,
  FloatImage *fim)
{
  // Create a file for writing the points
  ofstream fout(file, ios_base::out);

  // Create a solution data object for volume element computation
  SolutionData sd(xMedialModel->GetIterationContext(), xMedialModel->GetAtomArray());
  sd.ComputeIntegrationWeights();

  // We can not simply iterate over the interior points because there may also
  // be points outside of the range -1, 1 in xi, and we don't have an iterator
  // that handles them. So we just go atom by atom, outputting the coordinates
  MedialAtomIterator it = xMedialModel->GetAtomIterator();
  for( ; !it.IsAtEnd(); ++it)
    {
    // Get the atom out
    MedialAtom a = xMedialModel->GetAtomArray()[it.GetIndex()];

    // Write out the points for this atom
    for(double xi = xStart; xi < xEnd + 0.5 * xStep; xi += xStep)
      {
      SMLVec3d b = (xi < 0) ? a.xBnd[0].X : a.xBnd[1].X;
      SMLVec3d x = (1.0 - fabs(xi)) * a.X + fabs(xi) * b;
      fout << a.u << " " << a.v << " " << xi << " ";
      fout << x[0] << " " << x[1] << " " << x[2];

      // Write rho and r values
      fout << " " << a.xLapR;
      fout << " " << a.R;

      // Compute the volume element
      double velt = sd.ComputeVolumeElement(
        it.GetIndex(), xi, xStep, (xi == xStart || xi > xEnd - 0.5 * xStep));
      fout << " " << velt << endl;

      // Write the image value as well
      if(fim) fout << " " << fim->InterpolateNearestNeighbor(x);
      fout << endl;
      }
    }

  // Close the output stream
  fout.close();
}

/*
IMedialCoefficientMask *MedialPDE::CreatePCACoefficientMask(size_t nModes)
{
  size_t ncu, ncv;

  // Create a PCA mask
  IMedialCoefficientMask *xMask = NULL;

  // See if the numbers of coefficients match between current surface and matrix
  xSurface->GetNumberOfCoefficientsUV(ncu, ncv);
  if(ncu == ncuPCA && ncv == ncvPCA)
    {
    // Simple, use the matrix that was passed in
    xMask = new AffineAndPCACoefficientMask(mPCAMatrix, xSurface, nModes);
    }
  else
    {
    // We need to extract the coefficients from the source matrix to match the
    // current coefficients
    FourierSurface xSource(ncuPCA, ncvPCA), xTarget(ncu, ncv);
    vnl_matrix<double> mPCAChunk(
      mPCAMatrix.rows(), xSurface->GetNumberOfCoefficients());

    // For each subject in the matrix, remap using current coefficients
    for(size_t i = 0; i < mPCAMatrix.rows(); i++)
      {
      xSource.SetCoefficientArray(mPCAMatrix.get_row(i).data_block());
      for(size_t icu = 0; icu < ncu; icu++) for(size_t icv = 0; icv < ncv; icv++)
        for(size_t iComp = 0; iComp < xSurface->GetNumberOfDimensions(); iComp++)
          {
          double xVal = (icu >= ncuPCA || icv >= ncvPCA) ? 0.0 :
            xSource.GetCoefficient(icu, icv, iComp);
          xTarget.SetCoefficient(icu, icv, iComp, xVal);
          }
      mPCAChunk.set_row(i, xTarget.GetCoefficientArray());
      }

    // Pass this matrix to the PCA
    xMask = new Affine3DAndPCACoefficientMask(mPCAChunk, xSurface, nModes);
    }

  return xMask;
}
*/

/**
void MedialPDE::ReleasePCACoefficientMask(IMedialCoefficientMask *xMask)
{
  delete xMask;
}
*/

/***************************************************************************
 * SubdivisionMPDE Code
 * ----------------
 * This is the code for the subdivision surface based cm-rep
 **************************************************************************/
void SubdivisionMPDE::SubdivideMeshes(size_t iCoeffSub, size_t iAtomSub)
{
  // Get the subdivision surface medial model that is currently available
  SubdivisionMedialModel *smm = dynamic_cast<SubdivisionMedialModel *>(xMedialModel);
  if(!smm)
    throw MedialModelException(
      "SubdivisionMPDE given a cmrep not based on subdivision surfaces");

  // Get the coefficient-level mesh
  SubdivisionSurface::MeshLevel mlCoeffOld = smm->GetCoefficientMesh();
  vnl_vector<double> xCoeffOld = smm->GetCoefficientArray();
  vnl_vector<double> uCoeffOld = smm->GetCoefficientU();
  vnl_vector<double> vCoeffOld = smm->GetCoefficientV();

  SubdivisionSurface::MeshLevel mlCoeffNew;
  vnl_vector<double> xCoeffNew, uCoeffNew, vCoeffNew;

  // Select the coefficient mesh
  if(iCoeffSub > 0)
    {
    // Subdivide the coefficient-level mesh as requested
    SubdivisionSurface::RecursiveSubdivide(&mlCoeffOld, &mlCoeffNew, iCoeffSub);

    // Compute the coefficients for the new cm-rep
    xCoeffNew.set_size(mlCoeffNew.nVertices * 4);
    SubdivisionSurface::ApplySubdivision(
      xCoeffOld.data_block(), xCoeffNew.data_block(), 4, mlCoeffNew);

    // Interpolate u and v arrays
    uCoeffNew.set_size(mlCoeffNew.nVertices);
    vCoeffNew.set_size(mlCoeffNew.nVertices);
    SubdivisionSurface::ApplySubdivision(
      uCoeffOld.data_block(), uCoeffNew.data_block(), 1, mlCoeffNew);
    SubdivisionSurface::ApplySubdivision(
      vCoeffOld.data_block(), vCoeffNew.data_block(), 1, mlCoeffNew);

    // Set the subdivided mesh at the root
    mlCoeffNew.SetAsRoot();
    }
  else
    {
    // Use the same mesh as before
    xCoeffNew = xCoeffOld;
    mlCoeffNew = mlCoeffOld;
    uCoeffNew = uCoeffOld;
    vCoeffNew = vCoeffOld;
    }

  // Create a new mesh where the subdivided data will be stored
  SubdivisionMedialModel *smmNew = NULL;
  if(dynamic_cast<PDESubdivisionMedialModel *>(smm))
    smmNew = new PDESubdivisionMedialModel();
  else if (dynamic_cast<BruteForceSubdivisionMedialModel *>(smm))
    smmNew = new BruteForceSubdivisionMedialModel();

  // Set the mesh topology for the new mesh
  smmNew->SetMesh(
    mlCoeffNew, xCoeffNew, uCoeffNew, vCoeffNew, 
    smm->GetSubdivisionLevel() + iAtomSub - iCoeffSub, 0);

  // Now, we must solve the model. The way we do this depends on whether the
  // atom level has been changed or not. If it has, we will interpolate the
  // phi from the previous level to generate an initialization. Hopefully
  // this will help solve the PDE better...
  if(iAtomSub > 0)
    {
    // Get the mesh currently stored in the medial model
    SubdivisionSurface::MeshLevel mCurrentAtomLevel = smm->GetAtomMesh();

    // Set this mesh as root for subdivision
    mCurrentAtomLevel.SetAsRoot();

    // Subdivide to get the new mesh level
    SubdivisionSurface::MeshLevel mNewAtomLevel;
    SubdivisionSurface::RecursiveSubdivide(&mCurrentAtomLevel, &mNewAtomLevel, iAtomSub);

    // Get the array of phi values for the source
    vnl_vector<double> phiCurrent = smm->GetPhi();

    // Allocate the phi vector for the target
    vnl_vector<double> phiNew(mNewAtomLevel.nVertices);

    // Subdivide the source into the target
    SubdivisionSurface::ApplySubdivision(
      phiCurrent.data_block(), phiNew.data_block(), 1, mNewAtomLevel);

    // Put the new phi into the subdivision surface model
    smmNew->SetPhi(phiNew);

    // Solve the PDE without hints
    smmNew->ComputeAtoms();
    }
  else
    {
    // The solution is already there in the source-level mesh
    vnl_vector<double> xHint = smm->GetHintArray();
    smmNew->ComputeAtoms(xHint.data_block());
    }

  // Set the new mesh to replace the old mesh
  this->SetMedialModel(smmNew);
}


void SubdivisionMPDE::Remesh()
  {
  // Get the subdivision surface medial model that is currently available
  SubdivisionMedialModel *smm = dynamic_cast<SubdivisionMedialModel *>(xMedialModel);
  if(!smm)
    throw MedialModelException(
      "SubdivisionMPDE given a cmrep not based on subdivision surfaces");

  // Get the coefficient-level mesh
  SubdivisionSurface::MeshLevel mlCoeffOld = smm->GetCoefficientMesh(); 
  SubdivisionSurface::MeshLevel mlAtomsOld = smm->GetAtomMesh(); 


  // We need to get a list of coordinates for remeshing
  typedef vnl_vector_fixed<double, 3> Vec;
  Vec *X = new Vec[mlCoeffOld.nVertices];
  for(size_t i = 0; i < mlCoeffOld.nVertices; i++)
    for(size_t k = 0; k < 3; k++)
      X[i][k] = smm->GetCoefficient(i * 4 + k);








  // --- AFRONT CODE ---
  // Create a triangle mesh for Afront
  /*
  AFrontMesh amesh;
  for(size_t i = 0; i < mlCoeffOld.nVertices; i++)
    {
    gtb::tTriangleMeshVertex<float> v;
    for(size_t k = 0; k < 3; k++)
      v.point[k] = smm->GetCoefficient(i * 4 + k);
    amesh.verts.push_back(v);
    }
  for(size_t i = 0; i < mlCoeffOld.triangles.size(); i++)
    {
    gtb::TriangleMeshFace f;
    f.verts[0] = mlCoeffOld.triangles[i].vertices[0];
    f.verts[1] = mlCoeffOld.triangles[i].vertices[1];
    f.verts[2] = mlCoeffOld.triangles[i].vertices[2];
    amesh.faces.push_back(f);
    }
  std::vector<int> facemap, vertmap;
	amesh.IdentityMap(facemap, mlCoeffOld.triangles.size());
	amesh.IdentityMap(vertmap, mlCoeffOld.nVertices);
  amesh.build_structures(facemap, vertmap);

  // Run AFRONT
  do_tri_mesh(smm->GetSubdivisionLevel(), amesh); 
  */

  // Apply remeshing to the coefficient mesh
  mlCoeffOld.MakeDelaunay(X);

  // Second run (to be sure)
  cout << "RUN 2" << endl;
  mlCoeffOld.MakeDelaunay(X);

  // Check consistency
  SubdivisionSurface::CheckMeshLevel(mlCoeffOld);

  // Get rid of bad off-edge triangles
  SubdivisionSurface::MeshLevel mlCoeffNew;
  TriangleMeshGenerator tmg(&mlCoeffNew, mlCoeffOld.nVertices);
  for(size_t p = 0; p < mlCoeffOld.triangles.size(); p++)
    {
    size_t *v = mlCoeffOld.triangles[p].vertices;
    MedialAtom *A = smm->GetAtomArray();
    if(A[v[0]].flagCrest && A[v[1]].flagCrest && A[v[2]].flagCrest)
      {
      cerr << "Dropping triangle " << v[0] << " " << v[1] << " " << v[2] << endl;
      }
    else
      {
      tmg.AddTriangle(v[0], v[1], v[2]);
      }
    }
  tmg.GenerateMesh();

  // Find all unused vertices and create a vertex replacement rule
  vnl_vector<double> C(smm->GetCoefficientArray().data_block(), mlCoeffOld.nVertices * 4);
  vnl_vector<double> Cu(smm->GetCoefficientU().data_block(), mlCoeffOld.nVertices);
  vnl_vector<double> Cv(smm->GetCoefficientV().data_block(), mlCoeffOld.nVertices);

  std::vector<size_t> vcnt(mlCoeffOld.nVertices,0), vmap(mlCoeffOld.nVertices, 0);
  for(size_t i = 0; i < mlCoeffNew.triangles.size(); i++)
    for(size_t q = 0; q < 3; q++)
      vcnt[mlCoeffNew.triangles[i].vertices[q]] |= 1;
  size_t idx = 0;
  for(size_t j = 0; j < mlCoeffOld.nVertices; j++)
    {
    vmap[j] = vcnt[j] ? idx : NOID;    
    for(size_t c = 0; c < 4; c++)
      C[idx * 4 + c] = C[j * 4 + c];
    Cu[idx] = Cu[j];
    Cv[idx] = Cv[j];
    idx += vcnt[j];
    }

  // Now that we have a map from vertices to coefficients, we can shift up 
  // the vertex arrays
  for(size_t i = 0; i < mlCoeffNew.triangles.size(); i++)
    for(size_t q = 0; q < 3; q++)
      mlCoeffNew.triangles[i].vertices[q] = vmap[
        mlCoeffNew.triangles[i].vertices[q]];
  mlCoeffNew.nVertices = idx;
  mlCoeffNew.ComputeWalks();

  // Apply subdivision to the refined mesh
  try 
    {
    smm->SetMesh(mlCoeffNew, 
      C.extract(4*idx), 
      Cu.extract(idx), 
      Cv.extract(idx), 
      smm->GetSubdivisionLevel(), 0);
    }
  catch(...)
    {
    cerr << "Bad triangles encountered" << endl;
    }

  // Clean up
  delete X;
  }

void SubdivisionMPDE::BruteForceToPDE()
{
  // Get the brute force model, or throw exception
  BruteForceSubdivisionMedialModel *brute = 
    dynamic_cast<BruteForceSubdivisionMedialModel *>(xMedialModel);
  if(!brute)
    throw MedialModelException(
      "SubdivisionMPDE given a cmrep not based on subdivision surfaces");

  // Get coefficient arrays from the brute force model
  vnl_vector<double> C = brute->GetCoefficientArray();
  vnl_vector<double> u = brute->GetCoefficientU();
  vnl_vector<double> v = brute->GetCoefficientV();

  // Subdivide to get the atom-level coefficients
  vnl_vector<double> Ca(4 * brute->GetNumberOfAtoms());
  vnl_vector<double> ua(brute->GetNumberOfAtoms());
  vnl_vector<double> va(brute->GetNumberOfAtoms());
  for(size_t a = 0; a < brute->GetNumberOfAtoms(); a++)
    {
    Ca[4 * a] = brute->GetAtomArray()[a].X(0);
    Ca[4 * a + 1] = brute->GetAtomArray()[a].X(1);
    Ca[4 * a + 2] = brute->GetAtomArray()[a].X(2);
    Ca[4 * a + 3] = 0.0;
    ua[a] = brute->GetAtomArray()[a].u;
    va[a] = brute->GetAtomArray()[a].v;
    }
  
  // Clear the rho values
  for(size_t i = 0; i < C.size(); i+=4)
    C[i + 3] = 0.0;

  // Create a PDE model with the same mesh levels
  PDESubdivisionMedialModel *pde = new PDESubdivisionMedialModel();
  // pde->SetMesh(brute->GetCoefficientMesh(), C, u, v, 
  //  brute->GetSubdivisionLevel(), 0);
  //  cout << "START HERE" << endl;
  pde->SetMesh(brute->GetAtomMesh(), Ca, ua, va, 0, 0);
  // cout << "GOT HERE" << endl;

  // vnl_vector<double> hint(brute->GetAtomMesh().nVertices, 0.0);
  // try { pde->ComputeAtoms(hint.data_block()); } catch(...) {}

  // Compute the Laplace-Beltrami operator, using the current phi values
  vnl_vector<double> rho_atom = 
    pde->ComputeLBO(brute->GetPhi().data_block());
  vnl_vector<double> rho_ctl;

  // If there is subdivision in the model, we have to do a least square fit
  // to determine the rho values in the control-level mesh
  if(pde->GetSubdivisionLevel() > 0)
    {
    // Perform a least squares fit in order to get the control point rho
    // this involves solving the system W'Wx = W'y
    vnl_matrix<double> W = brute->GetAtomMesh().weights.GetDenseMatrix();
    vnl_matrix<double> WT = W.transpose();

    // Now we are going to set the columns in WT that correspond to boundary 
    // vertices to zeros. This corresponds to solving the least square problem
    // |JWx - Jy|^2, where J is a selection matrix (an identity matrix missing
    // the rows we don't want to fit). The solution is W'J'JWx = W'J'Jy, so
    // the matrix WT below is actually W'J'J. 
    for(size_t i = 0; i < WT.cols(); i++)
      if(!pde->GetAtomMesh().IsVertexInternal(i))
        WT.set_column(i, 0.0);

    // Solve the system for rho
    rho_ctl = vnl_qr<double>(WT * W).solve(WT * rho_atom);
    }
  else
    {
    // If no subdivision, the control-level rho is the same as the atom-level
    rho_ctl = rho_atom;
    cout << "RHO = " << rho_ctl << endl;
    }
  
  /*
  // Perform a least squares fit in order to get the control point rho
  // this involves solving the system W'Wx = W'y
  const ImmutableSparseMatrix<double> &W = brute->GetAtomMesh().weights;

  // Perform the ATA computation
  ImmutableSparseMatrix<double> WTW;
  ImmutableSparseMatrix<double>::InitializeATA(WTW, W);
  ImmutableSparseMatrix<double>::ComputeATA(WTW, W);

  // Compute W'y
  vnl_vector<double> WTy = W.MultiplyTransposeByVector(rho);

  // Solution vector
  vnl_vector<double> rho_ctl(brute->GetCoefficientMesh().nVertices, 0.0);

  // Solve the sparse system using Pardiso
  SymmetricPositiveDefiniteRealPARDISO solver;
  solver.SymbolicFactorization(WTW);
  solver.NumericFactorization(WTW);
  solver.Solve(WTy.data_block(), rho_ctl.data_block());
  */

  // Assign the rho values to the PDE model
  for(size_t i = 0; i < rho_ctl.size(); i++)
    pde->SetCoefficient(4*i+3, rho_ctl[i]);

  // Try solving the PDE problem
  try
    { pde->ComputeAtoms(brute->GetPhi().data_block()); }
  catch(MedialModelException &exc)
    { cout << "PDE failed: " << exc.what() << endl; }

  // Compare the rho in the atoms to actual rho
  double diff = 0.0; size_t ndiff = 0;
  for(size_t q = 0; q < pde->GetNumberOfAtoms(); q++)
    {
    if(pde->GetAtomMesh().IsVertexInternal(q))
      {
      double d = pde->GetAtomArray()[q].xLapR - rho_atom[q];
      diff += d * d;
      ndiff++;
      }
    }
  cout << "Mean squared difference in rho = " << diff / ndiff << endl;

  // Replace the brute force problem with the pde problem
  delete xMedialModel;
  xMedialModel = pde;
}

/***************************************************************************
 * CartesianMPDE Code
 * ----------------
 * This is the code for the rectangular grid-based cm-rep
 **************************************************************************/
CartesianMPDE::CartesianMPDE(unsigned int nBasesU, unsigned int nBasesV,
  unsigned int xResU, unsigned int xResV, double xFineScale,
  unsigned int xFineU, unsigned int xFineV)
{
  // Create a new Cartesian medial model (for now)
  this->xCartesianMedialModel =
    new CartesianMedialModel(xResU, xResV, xFineScale, xFineU, xFineV);

  // Associate the model with a Fourier surface
  FourierSurface *xSurface = new FourierSurface(nBasesU, nBasesV);
  xCartesianMedialModel->AdoptMedialSurface(xSurface);

  // Set the internal surface representation
  this->xMedialModel = xCartesianMedialModel;
}


/** Set the size of the evaluation grid */
void CartesianMPDE::SetGridSize(size_t nu, size_t nv, bool useHint)
{
  this->SetGridSize(nu, nv, 0, 0, 0.0, useHint);
}

/** Set the size of the evaluation grid with special sampling along edges */
void CartesianMPDE::SetGridSize(
  size_t nu, size_t nv, size_t eu, size_t ev, double eFactor, bool useHint)
{
  // Create a new model
  CartesianMedialModel *cmmnew = new CartesianMedialModel(nu, nv, eFactor, eu, ev);

  // Create a new surface
  size_t ncu, ncv;
  FourierSurface *surfold = dynamic_cast<FourierSurface *>(
    xCartesianMedialModel->GetMedialSurface());
  surfold->GetNumberOfCoefficientsUV(ncu, ncv);
  FourierSurface *surfnew = new FourierSurface(ncu, ncv);
  cmmnew->AdoptMedialSurface(surfnew);

  // Set the coefficients of the new surface
  cmmnew->SetCoefficientArray(xCartesianMedialModel->GetCoefficientArray());

  // Solve, potentially using the hint array
  if(useHint &&
    xCartesianMedialModel->GetGridU().size() == cmmnew->GetGridU().size() &&
    xCartesianMedialModel->GetGridV().size() == cmmnew->GetGridV().size())
    {
    CartesianMedialModel::Vec xHint = xCartesianMedialModel->GetHintArray();
    cmmnew->ComputeAtoms(xHint.data_block());
    }
  else
    cmmnew->ComputeAtoms();

  // Get rid of the old surface and model
  this->SetMedialModel(cmmnew);
  this->xCartesianMedialModel = cmmnew;
}


CartesianMPDE::CartesianMPDE(const char *file) :
  MedialPDE()
{
  // Load the model
  this->LoadFromParameterFile(file);
}

void CartesianMPDE::SetMedialModel(GenericMedialModel *model)
{
  // Cast the model to Cartesian model
  CartesianMedialModel *mc = dynamic_cast<CartesianMedialModel *>(model);
  if(mc == NULL)
    throw ModelIOException("Type mismatch in CartesianMPDE::SetMedialModel");

  // Call the parent's method
  MedialPDE::SetMedialModel(model);

  // Set my own pointer
  xCartesianMedialModel = mc;
}

IBasisRepresentation2D *
CartesianMPDE::GetMedialSurface()
{
  return xCartesianMedialModel->GetMedialSurface();
}

void CartesianMPDE::SetNumberOfCoefficients(unsigned int m, unsigned int n)
{
  // TODO: This is bad C++
  FourierSurface *xSurface = dynamic_cast<FourierSurface *>(
    xCartesianMedialModel->GetMedialSurface());
  xSurface->SetNumberOfCoefficients(m, n);
  xCartesianMedialModel->AdoptMedialSurface(xSurface);
  xCartesianMedialModel->ComputeAtoms();
}

vnl_vector<double> GetDoubleVectorInStides(const void *p, size_t n, size_t sz)
{
  const char *off = reinterpret_cast<const char *>(p);
  vnl_vector<double> vec(n, 0.0);
  for(size_t i = 0; i < n; i++, off+=sz)
    {
    vec[i] = *(reinterpret_cast<const double *>(off));
    }
  return vec;
}

void CartesianMPDE::ImportFromPointFile(
  const char *file, double xConstRho,
  bool flagFixRegularity, bool flagFitToRadius)
{
  vector<ImportAtom> atoms;

  ifstream fin(file, ios_base::in);
  while(!fin.eof())
    {
    ImportAtom atom;

    // Read a line from the file
    string line;
    getline(fin, line);

    // Drop comments
    if(line[0] == '#') continue;

    // Parse the line
    istringstream iss(line);

    // Read the atom
    iss >> atom.u; iss >> atom.v; iss >> atom.xi;
    iss >> atom.x; iss >> atom.y; iss >> atom.z;
    iss >> atom.rho; iss >> atom.r;

    // Allow rho override
    if(xConstRho != 0.0)
      atom.rho = xConstRho;

    if(!iss.fail())
      atoms.push_back(atom);
    }
  fin.close();

  // Do the actual import
  ImportFromPointData(atoms, flagFixRegularity, flagFitToRadius);
}


void CartesianMPDE::ImportFromPointData(
  const std::vector<ImportAtom> &atoms,
  bool flagFixRegularity,
  bool flagFitToRadius)
{
  // Get the number of import atoms
  size_t k = atoms.size();
  size_t n = xMedialModel->GetNumberOfAtoms();

  // Mapping from floats to ints, tolerance
  double scale = 100000000;

  // First, we need to see if the import data matches 1-for-1 the data in
  // the current model. 
  bool flagSameGrid = true;
  
  // Ok, now we will fit the x, y, z data to get model coefficients
  typedef vnl_vector<double> Vec;
  Vec ul = GetDoubleVectorInStides(&atoms[0].u, k, sizeof(ImportAtom));
  Vec vl = GetDoubleVectorInStides(&atoms[0].v, k, sizeof(ImportAtom));
  Vec xl = GetDoubleVectorInStides(&atoms[0].x, k, sizeof(ImportAtom));
  Vec yl = GetDoubleVectorInStides(&atoms[0].y, k, sizeof(ImportAtom));
  Vec zl = GetDoubleVectorInStides(&atoms[0].z, k, sizeof(ImportAtom));
  Vec wl = GetDoubleVectorInStides(&atoms[0].rho, k, sizeof(ImportAtom));
  Vec rl = GetDoubleVectorInStides(&atoms[0].r, k, sizeof(ImportAtom));

  // Get the underlying surface
  IBasisRepresentation2D *xSurface = xCartesianMedialModel->GetMedialSurface();

  // Perform the fitting on x, y, z and rho
  xSurface->FitToData(k, 0, ul.data_block(), vl.data_block(), xl.data_block());
  xSurface->FitToData(k, 1, ul.data_block(), vl.data_block(), yl.data_block());
  xSurface->FitToData(k, 2, ul.data_block(), vl.data_block(), zl.data_block());

  // The rho will only be fitted if there is no clean-up (regularization),
  // otherwise, we start with a very small rho and proceed to regularization
  if(flagFixRegularity)
    {
    Vec dummy(k, -0.01);
    xSurface->FitToData(k, 3, 
      ul.data_block(), vl.data_block(), dummy.data_block());
    }
  else
    {
    xSurface->FitToData(k, 3, ul.data_block(), vl.data_block(), wl.data_block());
    }

  // Compute the atoms, this may produce a bogus model, but that is OK
  try 
    {
    xMedialModel->ComputeAtoms();
    } 
  catch(MedialModelException &)
    {
    cerr << "Exception caught importing medial data" << endl;
    cerr << "   *** try using smaller rho ***" << endl;
    return;
    }

  // Create a mapping from u,v to atom index
  typedef std::pair<long, long> PairUV;
  std::map<PairUV, size_t> mapuv;
  for(size_t i = 0; i < k; i++)
    {
    PairUV key;
    key.first = static_cast<long>(atoms[i].u * scale + 0.5);
    key.second = static_cast<long>(atoms[i].v * scale + 0.5);
    mapuv[key] = i;
    }

  // Parse over all the atoms in the model to see if they all have a 
  // corresponding import atom
  MedialIterationContext *xGrid = xCartesianMedialModel->GetIterationContext();
  MedialAtom *xModelAtoms = xCartesianMedialModel->GetAtomArray();
  for(MedialAtomIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    MedialAtom &a = xModelAtoms[it.GetIndex()];
    PairUV key;
    key.first = static_cast<long>(a.u * scale + 0.5);
    key.second = static_cast<long>(a.v * scale + 0.5);
    if(mapuv.find(key) == mapuv.end()) 
      flagSameGrid = false;
    }

  // After the trivial fitting, we may want to do a closer fit to the input
  // data under various regularity constraints, so that we get a decent model
  if(flagFixRegularity)
    {
    // First thing that we need here is the XYZ value for every input point. If
    // the grid is the same as in the input, we can use the input directly, if 
    // not, we use the interpolation
    Vec xt(n, 0.0), yt(n, 0.0), zt(n, 0.0);
    for(size_t j = 0; j < n; j++)
      {
      MedialAtom &a = xModelAtoms[j];
      if(flagSameGrid)
        {
        PairUV key;
        key.first = static_cast<long>(a.u * scale + 0.5);
        key.second = static_cast<long>(a.v * scale + 0.5);
        size_t q = mapuv[key];
        xt[j] = atoms[q].x; yt[j] = atoms[q].y; zt[j] = atoms[q].z;
        }
      else
        {
        xSurface->EvaluateDerivative(a.u, a.v, 0, 0, 0, 1, &xt[j]);
        xSurface->EvaluateDerivative(a.u, a.v, 0, 0, 1, 1, &yt[j]);
        xSurface->EvaluateDerivative(a.u, a.v, 0, 0, 2, 1, &zt[j]);
        }
      }

    // Fit the model to the points
    FitPDEModelToPointData(xt.data_block(), yt.data_block(), zt.data_block());

    // Now, use the actual rho supplied to us
    xSurface->FitToData(k, 3, ul.data_block(), vl.data_block(), wl.data_block());
    xCartesianMedialModel->ComputeAtoms();
    }

  // Now, the model's XYZ geometry should be pretty reasonable, and we can
  // focus on getting a nice estimate of the radius field. We need to estimate
  // an initial radius field
  if(flagFitToRadius)
    {
    // Iniitalize the target radius field
    Vec rt(n, 0.0);

    if(flagSameGrid)
      {
      for(size_t j = 0; j < n; j++)
        {
        MedialAtom &a = xModelAtoms[j];
        rt[j] = atoms[mapuv[make_pair(
          static_cast<long>(a.u * scale + 0.5),
          static_cast<long>(a.v * scale + 0.5))]].r;
        }
      }
    else
      {
      // We need to fit to the radius field
      xSurface->FitToData(k, 3, ul.data_block(), vl.data_block(), rl.data_block());

      // Sample the interpolated radius values
      for(size_t j = 0; j < n; j++)
        {
        MedialAtom &a = xModelAtoms[j];
        xSurface->EvaluateDerivative(a.u, a.v, 0, 0, 3, 1, &rt[j]);
        }

      // Restore the surface to the rho fit
      xSurface->FitToData(k, 3, ul.data_block(), vl.data_block(), wl.data_block());
      }

    // Now solve the optimization problem
    FitPDEModelToRadius(rt.data_block());
    }
}

void CartesianMPDE::LoadFromDiscreteMRep(const char *file, double xRhoInit)
{
  ifstream fin(file, ios_base::in);
  int iend;

  // Vector to store atoms
  typedef vector<DiscreteAtom> AList;
  AList xAtoms;

  // Get the max u and v values
  unsigned int uMax = 0, vMax = 0;

  // Read atoms
  bool done = false;
  while(!fin.eof())
    {
    DiscreteAtom atom;

    // Read a line from the file
    string line;
    getline(fin, line);

    // Parse the line
    istringstream iss(line);

    // Read the atom
    iss >> atom.iu; iss >> atom.iv; iss >> iend;
    iss >> atom.x; iss >> atom.y; iss >> atom.z;
    iss >> atom.r; 

    if(iss.good())
      {
      xAtoms.push_back(atom);
      if(uMax < atom.iu) uMax = atom.iu;
      if(vMax < atom.iv) vMax = atom.iv;
      }
    }
  fin.close();

  // Scale by u and v to unit square
  for(AList::iterator it = xAtoms.begin(); it!=xAtoms.end(); ++it)
    { it->u = it->iu * 1.0 / uMax; it->v = it->iv * 1.0 / vMax; }

  // Create double arrays
  double *xx = new double[xAtoms.size()];
  double *yy = new double[xAtoms.size()];
  double *zz = new double[xAtoms.size()];
  double *uu = new double[xAtoms.size()];
  double *vv = new double[xAtoms.size()];
  double *rr = new double[xAtoms.size()];
  double *rd = new double[xAtoms.size()];

  for(unsigned int i = 0; i < xAtoms.size(); i++)
    {
    xx[i] = xAtoms[i].x;
    yy[i] = xAtoms[i].y;
    zz[i] = xAtoms[i].z;
    uu[i] = xAtoms[i].u;
    vv[i] = xAtoms[i].v;
    rd[i] = xAtoms[i].r;
    rr[i] = xRhoInit;
    }

  // Get the underlying surface
  IBasisRepresentation2D *xSurface = xCartesianMedialModel->GetMedialSurface();

  // Perform the fitting on x, y and z
  xSurface->FitToData(xAtoms.size(), 0, uu, vv, xx);
  xSurface->FitToData(xAtoms.size(), 1, uu, vv, yy);
  xSurface->FitToData(xAtoms.size(), 2, uu, vv, zz);

  // Fit the rho function to a constant or by default, zero
  xSurface->FitToData(xAtoms.size(), 3, uu, vv, rr);

  // Compute the atoms
  try {
    xCartesianMedialModel->ComputeAtoms();
    this->SaveToParameterFile("import01.cmrep");
  }
  catch(...) {
    cerr << "initially failed, trying to go on" << endl;
  }

  // Now, there may be all kinds of folding problems in the data. We want to
  // correct those problems by fitting the data under regularity constraints
  size_t n = xCartesianMedialModel->GetNumberOfAtoms();
  vnl_vector<double> target_x(n), target_y(n), target_z(n);
  for(size_t i = 0; i < n; i++)
    {
    target_x[i] = xCartesianMedialModel->GetAtomArray()[i].X[0];
    target_y[i] = xCartesianMedialModel->GetAtomArray()[i].X[1];
    target_z[i] = xCartesianMedialModel->GetAtomArray()[i].X[2];
    }
  FitPDEModelToPointData( 
    target_x.data_block(), target_y.data_block(), target_z.data_block());

  this->SaveToParameterFile("import02.cmrep");

  // If the rho is 0, set it for now to a small constant
  if(xRhoInit == 0.0)
    {
    for(unsigned int i = 0; i < xAtoms.size(); i++)
      rr[i] = -0.05;

    // Fit the rho function to a constant or by default, zero
    xSurface->FitToData(xAtoms.size(), 3, uu, vv, rr);
    xCartesianMedialModel->ComputeAtoms();

    // Fit to the radius values to get the radius field (temporary)
    xSurface->FitToData(xAtoms.size(), 3, uu, vv, rd);
    vnl_vector<double> xTargetRadius(xCartesianMedialModel->GetNumberOfAtoms());
    for(size_t i = 0; i < xTargetRadius.size(); i++)
      xSurface->EvaluateDerivative(
        xCartesianMedialModel->GetAtomArray()[i].u,
        xCartesianMedialModel->GetAtomArray()[i].v,
        0, 0, 3, 1, &xTargetRadius[i]);

    // Fit back to the rho value
    xSurface->FitToData(xAtoms.size(), 3, uu, vv, rr);

    // Now solve the optimization problem
    FitPDEModelToRadius(xTargetRadius.data_block());
    }

  // Clean up
  delete xx; delete yy; delete zz; delete uu; delete vv; delete rr; delete rd;
}




void CartesianMPDE::GenerateSampleModel()
{
  // Decide how many points to interpolate
  const unsigned int nSide = 21;
  const unsigned int nPoints = nSide * nSide;
  double uStep = 1.0 / (nSide - 1);

  // Allocate arrays of points and coordinates
  double xPoints[nPoints], yPoints[nPoints], zPoints[nPoints], rhoPoints[nPoints];
  double uPoints[nPoints], vPoints[nPoints];

  // Create an array of points
  unsigned int i = 0;
  for(unsigned int u = 0; u < nSide; u++)
    for(unsigned int v = 0; v < nSide; v++)
      {
      double uu = u * uStep, vv = v * uStep;

      uPoints[i] = uu;
      vPoints[i] = vv;
      xPoints[i] = 0.5 * uu + 0.25;
      yPoints[i] = vv;
      zPoints[i] = ((uu - 0.5) * (uu - 0.5) + (vv - 0.5) * (vv - 0.5)) * 0.25;
      rhoPoints[i] = -0.35;
      ++i;
      }

  // Get the underlying surface
  IBasisRepresentation2D *xSurface = xCartesianMedialModel->GetMedialSurface();

  // Peform the fit
  xSurface->FitToData(nPoints, 0, uPoints, vPoints, xPoints);
  xSurface->FitToData(nPoints, 1, uPoints, vPoints, yPoints);
  xSurface->FitToData(nPoints, 2, uPoints, vPoints, zPoints);
  xSurface->FitToData(nPoints, 3, uPoints, vPoints, rhoPoints);

  // Generate the atoms
  xCartesianMedialModel->ComputeAtoms();
}


// This is just one approach to doing this. Let's hope that it works
void CartesianMPDE::
SampleReferenceFrameImage(FloatImage *imgInput, FloatImage *imgOutput, size_t zSamples)
{
  // Get the internal images
  typedef FloatImage::WrapperType::ImageType ImageType;
  ImageType::Pointer iInput = imgInput->xImage->GetInternalImage();
  ImageType::Pointer iOutput = imgOutput->xImage->GetInternalImage();

  // Get the image dimensions
  size_t m = xCartesianMedialModel->GetNumberOfUPoints();
  size_t n = xCartesianMedialModel->GetNumberOfVPoints();

  // The image dimensions must match the grid size
  itk::Size<3> szInput = iInput->GetBufferedRegion().GetSize();
  if ( szInput[0] != m || szInput[1] != n || szInput[2] != zSamples * 2 + 1)
    {
    cerr << "Input image dimensions do not match mpde size" << endl;
    return;
    }

  // Create a flat image matched up with the internal iterator
  size_t npix = xMedialModel->GetNumberOfInternalPoints(zSamples);
  vnl_vector<float> xpix(npix, 0.0f);

  MedialInternalPointIterator it(xMedialModel->GetIterationContext(), zSamples);
  for( ; !it.IsAtEnd(); ++it)
    {
    MedialAtom &xAtom = xMedialModel->GetAtomArray()[it.GetAtomIndex()];
    itk::Index<3> idx;
    idx[0] = xAtom.uIndex;
    idx[1] = xAtom.vIndex;
    idx[2] = it.GetBoundarySide() ? zSamples - it.GetDepth() : zSamples + it.GetDepth();
    xpix[it.GetIndex()] = iInput->GetPixel(idx);
    }

  // Generate a VTK cell grid
  vtkUnstructuredGrid *cells = ExportVolumeMeshToVTK(xMedialModel, zSamples);

  // Create a VTk locator for these cells
  vtkPointLocator *loc = vtkPointLocator::New();
  loc->SetDataSet(cells);
  loc->BuildLocator();

  // Iterate over the voxels in the output image. The image should be initialized
  // to a mask, so that we don't waste time on outside pixels
  itk::ImageRegionIteratorWithIndex<ImageType> itOut(
    iOutput, iOutput->GetBufferedRegion());
  for( ; !itOut.IsAtEnd(); ++itOut)
    {
    // Only do something at masked points
    if(itOut.Get() != 0.0f)
      {
      // Convert this voxel into a point
      ImageType::PointType pVoxel;
      iOutput->TransformIndexToPhysicalPoint(itOut.GetIndex(), pVoxel);

      // Define the extent
      double x0 = pVoxel[0] - 0.5 * iOutput->GetSpacing()[0];
      double x1 = pVoxel[0] + 0.5 * iOutput->GetSpacing()[0];
      double y0 = pVoxel[1] - 0.5 * iOutput->GetSpacing()[1];
      double y1 = pVoxel[1] + 0.5 * iOutput->GetSpacing()[1];
      double z0 = pVoxel[2] - 0.5 * iOutput->GetSpacing()[2];
      double z1 = pVoxel[2] + 0.5 * iOutput->GetSpacing()[2];

      // Find all samples inside that pixel
      double xMax = 0.0;
      MedialInternalPointIterator it(xMedialModel->GetIterationContext(), zSamples);
      for( ; !it.IsAtEnd(); ++it)
        {
        SMLVec3d x = GetInternalPoint(it, xMedialModel->GetAtomArray());
        if(x[0] <= x1 && x[0] > x0 && x[1] <= y1 && x[1] > y0 && x[2] <= z1 && x[2] > z0)
          {
          xMax = xMax < xpix[it.GetIndex()] ? xpix[it.GetIndex()] : xMax;
          }
        }

      // Locate the cell that includes this point
      vnl_vector<double> vox = pVoxel.GetVnlVector();
      size_t iClosest = (size_t)
        loc->FindClosestPoint(vox(0), vox(1), vox(2));

      // Convert this to a pixel index
      // itOut.Set(xpix[iClosest]);
      itOut.Set(xMax);

      }
    }

  loc->Delete();
  cells->Delete();
/*
  // Compute the image region that includes the mask
  itk::ImageRegion<3> xMaskRegion;

  // Define an ITK transform for mapping the points
  typedef itk::BSplineDeformableTransform<double, 3, 3> TransformType;
  TransformType::Pointer tps = TransformType::New();

  // Define the lanmark arrays
  size_t nLMCuts = 2;
  size_t nLandmarks = xMedialModel->GetAtomGrid()->GetNumberOfInternalPoints(nLMCuts);

  SMLVec3d *lSource = new SMLVec3d[nLandmarks];
  SMLVec3d *lTarget = new SMLVec3d[nLandmarks];

  // Define the landmarks
  MedialInternalPointIterator *it =
    xMedialModel->GetAtomGrid()->NewInternalPointIterator(nLMCuts);
  for( size_t i=0 ; !it->IsAtEnd(); ++(*it), ++i)
    {
    // The source are the positions in patient space
    lSource[i] = GetInternalPoint(it, xMedialModel->GetAtomArray());

    // Convert the source point to an image index
    ImageType::PointType pt; ImageType::IndexType idx;
    pt[0] = lSource[i][0]; pt[1] = lSource[i][1]; pt[2] = lSource[i][2];
    iOutput->TransformPhysicalPointToIndex(pt, idx);
    UpdateRegion(xMaskRegion, idx, i==0);

    // The target are the positions in cm-rep coordinate system
    MedialAtom &xAtom = xMedialModel->GetAtomArray()[it->GetAtomIndex()];
    lTarget[i][0] = round(xAtom.u * (m - 1));
    lTarget[i][1] = round(xAtom.v * (n - 1));
    lTarget[i][2] = it->GetBoundarySide() ?
      zSamples - it->GetDepth() : zSamples + it->GetDepth();
    }

  // Pass the ladmarks to the transform
  xMaskRegion.SetIndex(0, 0);
  xMaskRegion.SetIndex(1, 0);
  xMaskRegion.SetIndex(2, 0);
  xMaskRegion.SetSize(0, 6);
  xMaskRegion.SetSize(1, 10);
  xMaskRegion.SetSize(2, 6);

  tps->SetGridRegion(xMaskRegion);
  // tps->SetGridSpacing(iOutput->GetSpacing() / 10.0);
  // tps->SetGridOrigin(iOutput->GetOrigin());

  cout << "Selected Region: " << xMaskRegion << endl;
  cout << "Number of params: " << tps->GetNumberOfParameters() << endl;

  LeastSquaresFit(tps, nLandmarks, lSource, lTarget);

  // Iterate over the voxels in the output image. The image should be initialized
  // to a mask, so that we don't waste time on outside pixels
  itk::ImageRegionIteratorWithIndex<ImageType> itOut(
    iOutput, iOutput->GetBufferedRegion());
  for( ; !itOut.IsAtEnd(); ++itOut)
    {
    // Only do something at masked points
    if(itOut.Get() != 0.0f)
      {
      // Convert this voxel into a point
      ImageType::PointType pVoxel;
      iOutput->TransformIndexToPhysicalPoint(itOut.GetIndex(), pVoxel);

      // Map this point into the other image space
      TransformType::OutputPointType pReference =
        tps->TransformPoint(pVoxel);

      // Interpolate the second image at this point
      itOut.Set(imgInput->Interpolate(
          SMLVec3d(pReference[0], pReference[1], pReference[2])));
      }
    }
*/
}

void CartesianMPDE::
SampleImage(FloatImage *fiInput, FloatImage *fiOut, size_t zSamples)
{
  // Allocate the flattened intensity output image
  typedef FloatImage::WrapperType::ImageType ImageType;
  ImageType::Pointer imgOutSmall = ImageType::New();
  ImageType::RegionType regOutput;

  // Get the image dimensions
  size_t m = xCartesianMedialModel->GetNumberOfUPoints();
  size_t n = xCartesianMedialModel->GetNumberOfVPoints();

  // Initialize the output region size
  regOutput.SetSize(0, m);
  regOutput.SetSize(1, n);
  regOutput.SetSize(2, zSamples * 2 + 1);

  // Continue creating the image
  imgOutSmall->SetRegions(regOutput);
  imgOutSmall->Allocate();
  imgOutSmall->FillBuffer(0.0f);

  // Get the input image
  ImageType::Pointer imgInput = fiInput->xImage->GetInternalImage();

  // Create the output images for U, V, Tau coordinates
  ImageType::Pointer imgCoord[3];
  for(size_t d = 0; d < 3; d++)
    {
    imgCoord[d] = ImageType::New();
    imgCoord[d]->SetRegions(imgInput->GetBufferedRegion());
    imgCoord[d]->SetSpacing(imgInput->GetSpacing());
    imgCoord[d]->SetOrigin(imgInput->GetOrigin());
    imgCoord[d]->Allocate();
    imgCoord[d]->FillBuffer(-2.0);
    }

  // Create a nearest neighbor interpolator
  //typedef itk::NearestNeighborInterpolateImageFunction<
  //  ImageType, double> InterpolatorType;
  // typedef itk::BSplineInterpolateImageFunction<
  //  ImageType, double> InterpolatorType;
  typedef itk::LinearInterpolateImageFunction<
    ImageType, double> InterpolatorType;
  InterpolatorType::Pointer fInterp = InterpolatorType::New();
  fInterp->SetInputImage(imgInput);

  // Get an internal point iterator
  MedialInternalPointIterator itPoint(xMedialModel->GetIterationContext(), zSamples-1);
  for(; !itPoint.IsAtEnd(); ++itPoint)
    {
    // Interpolate the position in 3-space
    SMLVec3d xPoint = GetInternalPoint(itPoint, xMedialModel->GetAtomArray());

    // Get the corresponding medial atom
    MedialAtom &xAtom = xMedialModel->GetAtomArray()[itPoint.GetAtomIndex()];

    // Get the i and j coordinates of the atom (using cartesian logic)
    size_t i = xAtom.uIndex;
    size_t j = xAtom.vIndex;
    size_t k = itPoint.GetBoundarySide() ?
      zSamples - itPoint.GetDepth() :
       zSamples + itPoint.GetDepth();

    // Create an index into the output image
    ImageType::IndexType idxTarget;
    idxTarget.SetElement(0, i);
    idxTarget.SetElement(1, j);
    idxTarget.SetElement(2, k);

    // Get the point at which to sample the image
    double tau = itPoint.GetDepth() * 1.0 / itPoint.GetMaxDepth();
    SMLVec3d Z =
      xAtom.X + tau * xAtom.R * xAtom.xBnd[itPoint.GetBoundarySide()].N;

    // Create a point and a continuous index
    itk::Point<double, 3> ptZ(Z.data_block());
    itk::ContinuousIndex<double, 3> idxZ;
    imgInput->TransformPhysicalPointToContinuousIndex(ptZ, idxZ);
    float f = fInterp->EvaluateAtContinuousIndex(idxZ);

    // Print some random statistics
    if(rand() % 4000 == 0)
      {
      cout << "MCoord : " << xAtom.u << ", " << xAtom.v << ", " << tau << "; ";
      cout << "ZCoord : " << Z << "; ";
      cout << "ICoord : " << (int)(idxZ[0]+0.5) << ", " << (int)(idxZ[1] + 0.5) << ", " << (int)(idxZ[2]+0.5) << "; ";
      cout << "PCoord : " << ptZ[0] << ", " << ptZ[1] << ", " << ptZ[2] << "; ";
      cout << "IVal : " << f << endl;
      }

    // Sample the input image
    // float f = imgInput->xImage->Interpolate(Z[0], Z[1], Z[2], 0.0f);
    imgOutSmall->SetPixel(idxTarget, f);

    // This is junk!
    SMLVec3d Z1 = Z + 0.5;
    ImageType::IndexType idxCoord;
    itk::Point<double, 3> ptZ1(Z1.data_block());
    imgCoord[0]->TransformPhysicalPointToIndex(ptZ1, idxCoord);
    imgCoord[0]->SetPixel(idxCoord, xAtom.u);
    imgCoord[1]->SetPixel(idxCoord, xAtom.v);
    imgCoord[2]->SetPixel(idxCoord, itPoint.GetBoundarySide() ? tau : -tau);

    //if(rand() % 80 == 0)
    //  {
    //  cout << "Sample " << Z << " value " << f << endl;
    //  }

    // If the index is at the edge, there are two output pixels
    if(itPoint.IsEdgeAtom())
      {
      idxTarget.SetElement(2, zSamples * 2 - k);
      imgOutSmall->SetPixel(idxTarget, f);
      }
    }

  // Store the output image
  fiOut->xImage->SetInternalImage(imgOutSmall);
  fiOut->xGradient[0]->SetInternalImage(imgCoord[0]);
  fiOut->xGradient[1]->SetInternalImage(imgCoord[1]);
  fiOut->xGradient[2]->SetInternalImage(imgCoord[2]);
}


/**
 * MEDIAL PCA CODE
 */
void MedialPCA::AddSample(MedialPDE *pde)
{
  // Make sure that the number of coefficients matches
  if(xModelCoefficients.size())
    assert(xModelCoefficients.back().size() == pde->xMedialModel->GetNumberOfCoefficients());

  // Get the coefficients from this medial PDE
  xModelCoefficients.push_back(pde->xMedialModel->GetCoefficientArray());

  // Get the hint array from the medial PDE
  xHints.push_back(pde->xMedialModel->GetHintArray());

  // Make a copy of the intensity values associated with this PDE
  if(pde->flagIntensityPresent)
    {
    // Make a copy of the image
    FloatImage *img = new FloatImage();
    img->xImage->SetInternalImage(
      pde->imgIntensity.xImage->GetInternalImage());
    xAppearance.push_back(img);
    }
}

MedialPCA::MedialPCA()
{
  xPCA = NULL;
  xAppearancePCA = NULL;
}

MedialPCA::~MedialPCA()
{
  if(xPCA) delete xPCA;
  if(xAppearancePCA) delete xAppearancePCA;
}

void MedialPCA::ComputePCA(MedialPDE *mpde)
{
  // The non-shape information must be cleaned from the data using
  // something like the generalized Procrustes method. We begin by
  // interpolating each mpde
  size_t nSamples = xModelCoefficients.size();

  // These are the parameters to the GPA method
  Mat *A = new Mat[nSamples];
  Mat *R = new Mat[nSamples];
  Vec *t = new Vec[nSamples];
  double *s = new double[nSamples];
  size_t i, j;

  // Generate same parameters for the boundary-based alignment and PCA
  Mat *Abnd = new Mat[nSamples];
  Mat *Rbnd = new Mat[nSamples];
  Vec *tbnd = new Vec[nSamples];
  double *sbnd = new double[nSamples];

  // Populate the input matrices
  GenericMedialModel *xMedialModel = mpde->xMedialModel;
  MedialIterationContext *xGrid = xMedialModel->GetIterationContext();

  for(i = 0; i < nSamples; i++)
    {
    // Solve for this surface using these coefficients
    cout << "Compute subject " << i << endl;
    xMedialModel->SetCoefficientArray(xModelCoefficients[i]);
    xMedialModel->ComputeAtoms(xHints[i].data_block());

    // Initialize the A matrix
    A[i].set_size(xMedialModel->GetNumberOfAtoms(), 3);
    R[i].set_size(3,3);
    t[i].set_size(3);

    // Parse over all the boundary sites
    for(MedialAtomIterator it(xGrid); !it.IsAtEnd(); ++it)
      A[i].set_row(it.GetIndex(), xMedialModel->GetAtomArray()[it.GetIndex()].X);

    // Same for boundary atoms
    Abnd[i].set_size(xMedialModel->GetNumberOfBoundaryPoints(), 3);
    Rbnd[i].set_size(3,3);
    tbnd[i].set_size(3);

    for(MedialBoundaryPointIterator bit(xGrid); !bit.IsAtEnd(); ++bit)
      Abnd[i].set_row(bit.GetIndex(), GetBoundaryPoint(bit, xMedialModel->GetAtomArray()).X);
    }

  // cout << "Testing Procrustes" << endl;
  // TestProcrustes(A[0]);

  // Run the procrustes method
  cout << "Procrustenating..." << endl;
  GeneralizedProcrustesAnalysis(nSamples, A, R, t, s);

  // Create an affine mapper
  const AffineTransformDescriptor *affDesc = xMedialModel->GetAffineTransformDescriptor();

  // Rotate each of the models into the common coordinate frame
  for(i = 0; i < nSamples; i++)
    {
    // Transform the coefficients under the affine transform
    Vec xAlignCoeff = affDesc->ApplyAffineTransform(
      xModelCoefficients[i], s[i] * R[i].transpose(), t[i], Vec(3, 0.0));

    // Create a new MedialPDE (?)
    // CartesianMPDE xJunk(8,12,32,80);
    // xJunk.xMedialModel->SetCoefficientArray(xAlignCoeff);
    // xJunk.Solve();
    cout << "Compute aligned " << i << endl;
    xMedialModel->SetCoefficientArray(xAlignCoeff);
    xMedialModel->ComputeAtoms(xHints[i].data_block());
    }

  // Create a principal components object

  // Compute the mean shape and the covariance matrix on the fourier
  // parameters. Since the Fourier basis is orthonormal, doing PCA on the
  // surface and on the Fourier components is identical
  size_t m = xModelCoefficients[0].size();
  size_t n = nSamples;

  // Populate the data matrix
  xDataShape.set_size(n, m);
  for(i = 0; i < n; i++) for(j = 0; j < m; j++)
    xDataShape[i][j] = xModelCoefficients[i][j];

  // Compute the principal components of xDataShape
  if(xPCA) delete xPCA;
  xPCA = new PrincipalComponents(xDataShape);
  xPCALocation.set_size(xPCA->GetNumberOfModes());
  xPCALocation.fill(0.0);

  // If there is appearance information, also compute the appearance PCA
  if(xAppearance.size() == n)
    {
    // Create a matrix of appearance values
    typedef FloatImage::WrapperType::ImageType ImageType;
    ImageType::Pointer imgFirst =
      xAppearance.front()->xImage->GetInternalImage();
    size_t mpix = imgFirst->GetBufferedRegion().GetNumberOfPixels();
    xDataAppearance.set_size(n, mpix);

    // Populate the matrix from the image
    for(i = 0; i < n; i++)
      {
      float *pix = xAppearance[i]->xImage->GetInternalImage()->GetBufferPointer();
      for(j = 0; j < mpix; j++)
        xDataAppearance(i, j) = pix[j];
      }

    // Compute the PCA from the image array
    if(xAppearancePCA) delete xAppearancePCA;
    xAppearancePCA = new PrincipalComponents(xDataAppearance);
    xAppearancePCALocation.set_size(xAppearancePCA->GetNumberOfModes());
    xAppearancePCALocation.fill(0.0);
    }

  // Report the principal components
  cout << "EIGENVALUES: " << endl;
  for(size_t l = 0; l < xPCA->GetNumberOfModes(); l++)
    { cout << xPCA->GetEigenvalue(l) << " "; }
  cout << endl;

  /*
  // Compute the boundary PCA (just to compare)
  GeneralizedProcrustesAnalysis(xSurfaces.size(), Abnd, Rbnd, tbnd, sbnd);
  Mat xBndMat(n, Abnd[0].rows() * 3);
  for(i = 0; i < n; i++) for(j = 0; j < Abnd[0].rows(); j++)
    {
    size_t k = j * 3;
    Vec p = tbnd[i] + sbnd[i] * Rbnd[i].transpose() * Abnd[i].get_row(j);
    xBndMat[i][k] = p[0];
    xBndMat[i][k+1] = p[1];
    xBndMat[i][k+2] = p[2];
    if(j == 0)
      cout << p << endl;
    }

  PrincipalComponents xPCABnd(xBndMat);
  cout << "BND EIGENVALUES: " << endl;
  for(size_t l = 0; l < xPCABnd.GetNumberOfModes(); l++)
    { cout << xPCABnd.GetEigenvalue(l) << " "; }
  cout << endl;
  */

  // Project everything into the PCA
  for(i = 0; i < n; i++)
    {
    Vec y = xPCA->MapToCoefficientSpace(xDataShape.get_row(i), 8);
    cout << "Inp. Mah. Dist.: " << i << " " << y.two_norm() << " " << y << endl;
    }

  /*
  // Generate shapes along the first four eigenvectors
  for(i = 0; i < 4; i++)
    {
    for(int y = 0; y < 60; y++)
      {
      // A shape!
      double c = 0.1 * (y - 30);
      Vec z = xPCA->MapToFeatureSpace(i, c);

      cout << "MODE " << i << ", COEFF " << c << ", LAMBDA " <<
        xPCA->GetEigenvalue(i) << endl;

      MedialPDE xJunk(8,12,32,80);
      xJunk.xSurface = new FourierSurface(*xSurfaces[0]);
      xJunk.xMedialModel->SetMedialSurface(xJunk.xSurface);
      xJunk.xSurface->SetRawCoefficientArray(z.data_block());
      xJunk.Solve();

      ostringstream sJunk1, sJunk2;
      sJunk1 << "/tmp/pcamed_m" << i << "_" << (y / 10) << (y % 10) << ".vtk";
      sJunk2 << "/tmp/pcabnd_m" << i << "_" << (y / 10) << (y % 10) << ".vtk";
      xJunk.SaveVTKMesh(sJunk1.str().c_str(), sJunk2.str().c_str());
      }
    }
  */
}

// Move along a given mode a certain number of S.D.
void MedialPCA::SetFSLocationToMean()
{
  xPCALocation.fill(0.0);
  if(xAppearancePCALocation.size())
    xAppearancePCALocation.fill(0.0);
}

// Export shape matrix
void MedialPCA::ExportShapeMatrix(const char *filename)
{
  WriteMatrixFile(xDataShape, filename);
}


// Move along a given mode a certain number of S.D.
void MedialPCA::SetFSLocation(unsigned int iMode, double xSigma)
{
  xPCALocation(iMode) = xSigma;
  if(xAppearancePCALocation.size())
    xAppearancePCALocation(iMode) = xSigma;
}

// Generate a sample at the current location
void MedialPCA::GetShapeAtFSLocation(MedialPDE *target)
{
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;

  // Compute the shape space point
  Vec z = xPCA->MapToFeatureSpace(xPCALocation);

  // Compute the shape
  target->xMedialModel->SetCoefficientArray(z);
  target->Solve();

  // Compute the feature image
  if(xAppearance.size() == xModelCoefficients.size())
    {
    FloatImage *imgJunk = new FloatImage();
    FloatImage::WrapperType::ImageType *img = imgJunk->xImage->GetInternalImage();
    img->SetRegions(
      xAppearance[0]->xImage->GetInternalImage()->GetBufferedRegion());
    img->Allocate();

    Vec z1 = xAppearancePCA->MapToFeatureSpace(xAppearancePCALocation);
    for(size_t i = 0; i < z1.size(); i++)
      img->GetBufferPointer()[i] = z1[i];

    target->SetIntensityImage(imgJunk);
    delete imgJunk;
    }
}

/*
void RenderMedialPDE(MedialPDE *model)
{
  // Create a renderer
  PDESplineRenderer *rndPDESpline = new PDESplineRenderer(model->xMedialModel);

  // Initialize the GL environment
  char *argv[] = {"test",NULL};
  GLDisplayDriver::init(1, argv);

  // Add the important renderers
  GLDisplayDriver::addRenderer(new DefaultLightRenderer(), GLDisplayDriver::EYE);
  GLDisplayDriver::addRenderer(new StarfieldRenderer(), GLDisplayDriver::UNIT);

  // Add our spline renderer
  GLDisplayDriver::addRenderer(rndPDESpline);

  // Get the affine descriptor to compute the center of rotation
  SMLVec3d xCenter = model->GetMedialModel()->GetCenterOfRotation();

  // Put the display into Talairach coordinates
  GLDisplayDriver::center[0] = xCenter[0];
  GLDisplayDriver::center[1] = xCenter[1];
  GLDisplayDriver::center[2] = xCenter[2];
  GLDisplayDriver::scale = 20;

  // Start GLUT
  glutMainLoop();
}
*/

