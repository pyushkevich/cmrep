#include "ScriptInterface.h"
#include "CartesianMedialModel.h"
#include "BasisFunctions2D.h"
#include "ITKImageWrapper.h"
#include "MedialAtomGrid.h"
#include "OptimizationTerms.h"
#include "Procrustes.h"
#include "Registry.h"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm>

using namespace std;

// An inline function to update max counters
inline void UpdateMax(double xValue, double &xMaxValue)
{
  if(xValue > xMaxValue)
    xMaxValue = xValue;
}

class DerivativeTestQuantities {
public:
  typedef std::map<std::string, double> MapType;
  typedef MapType::iterator MapIter;
  typedef MapType::const_iterator ConstMapIter;

  DerivativeTestQuantities()
    { slenmax = 0; }

  void Update(const std::string &s, MapType &xMap, double v)
    {
    MapIter it = xMap.find(s);
    if(it == xMap.end())
      {
      xMap.insert(std::make_pair(s, v));
      slenmax = std::max(slenmax, s.size());
      }
    else
      it->second = std::max(it->second, v);
    }

  void Update(const std::string &s, double dv, double v1, double v2, double eps)
    {
    double eabs = fabs(dv - (v1 - v2) * 0.5 / eps);
    double rms = sqrt((dv * dv + ((v1 - v2) * 0.5 / eps) * ((v1 - v2) * 0.5 / eps)) / 2);
    double erel = eabs / (rms + eps);
    this->Update(s, xAbsError, eabs); 
    this->Update(s, xRelError, erel); 
    }

  void Update(const std::string &s, double dv, double v1, double v2, double v3, double v4, double eps)
    {
    // Estimate the derivative up to the fourth order of epsilon
    double fd = -(v1 + 27 * (v3 - v2) - v4) / (24 * eps);    
    double eabs = fabs(dv - fd);
    double mean = 0.5 * (fabs(dv) + fabs(fd)) + eps;
    double erel = eabs / mean;
    this->Update(s, xAbsError, eabs);
    this->Update(s, xRelError, erel);
    }

  void Update(const std::string &s, SMLVec3d dv, SMLVec3d v1, SMLVec3d v2, double eps)
    {
    // Absolute error
    double eabs = (dv - (v1 - v2) * 0.5 / eps).two_norm();

    // Get the RMS of analytical and numerical gradients
    double rms = 
      sqrt((((v1 - v2) * 0.5 / eps).squared_magnitude() + dv.squared_magnitude()) / 6);

    // Get the relative error
    double erel = eabs / (rms + eps);
    
    // The way to compute e-rel below is wrong! comparing derivative
    // to the function value
    // double erel = eabs / (0.5 * (v1 + v2).two_norm() + eps);
    this->Update(s, xAbsError, eabs); 
    this->Update(s, xRelError, erel); 
    }

  void Update(const std::string &s, SMLVec3d dv, 
    SMLVec3d v1, SMLVec3d v2, SMLVec3d v3, SMLVec3d v4, double eps)
    {
    SMLVec3d fd = -(v1 + 27. * (v3 - v2) - v4) / (24 * eps); 
    double eabs = (dv - fd).two_norm();
    double mean = 0.5 * (dv.two_norm() + fd.two_norm()) + eps;
    double erel = eabs / mean;
    this->Update(s, xAbsError, eabs); 
    this->Update(s, xRelError, erel); 
    }

  void Update(const DerivativeTestQuantities &src)
    {
    for(ConstMapIter it = src.xAbsError.begin(); it != src.xAbsError.end(); ++it)
      this->Update(it->first, xAbsError, it->second);
    for(ConstMapIter it = src.xRelError.begin(); it != src.xRelError.end(); ++it)
      this->Update(it->first, xRelError, it->second);
    }

  bool Test(double xMaxAbsError, double xMaxRelError)
    {
    for(ConstMapIter it = xAbsError.begin(); it != xAbsError.end(); ++it)
      if(it->second > xMaxAbsError) return false;
    for(ConstMapIter it = xRelError.begin(); it != xRelError.end(); ++it)
      if(it->second > xMaxRelError) return false;
    return true;
    }

  void PrintReport() 
    {
    for(ConstMapIter it = xAbsError.begin(); it != xAbsError.end(); ++it)
      cout << std::setw(slenmax) << it->first 
        << " :  abs = " << std::setw(12) << it->second 
        << ";  rel = "  << std::setw(12) << xRelError[it->first] << endl;
    }

private:
  MapType xAbsError, xRelError;

  size_t slenmax;
};

int TestGradientComputation(
  GenericMedialModel *xSolver, CoefficientMapping *xMapping, 
  vnl_vector<double> P0, int nRandVar)
{
  size_t iVar;

  // Timers for different routines
  CodeTimer tCentral, tAnalytic, tInitialSolver;
  
  // The epsilon for central difference tests
  double master_eps = 1.0e-5;

  // The number of atoms in the solver
  size_t nAtoms = xSolver->GetNumberOfAtoms();

  // The number of parameters in the source coefficient space
  size_t nParams = xMapping->GetNumberOfParameters();
  size_t nCoeffs = xMapping->GetNumberOfCoefficients();
  size_t nVar = nParams + nRandVar;
  
  // Holder for all the test quantities. This object will report the maximum
  // value for each quantity over all the tests
  DerivativeTestQuantities dtqGlobal;

  cout << "*****************************************" << endl;
  cout << "Comparing Numeric Derivatives to Analytic" << endl;
  cout << "*****************************************" << endl;

  // Create a list of variational surfaces: some components, some random
  // variations. There are two basis sets: one in parameter space (PS) 
  // another in coefficient space (CS)
  vnl_matrix<double> xVariationalBasisPS(nParams + nRandVar, nParams, 0.0);
  vnl_matrix<double> xVariationalBasisCS(nParams + nRandVar, nCoeffs, 0.0);
  vector<string> xVariationNames;

  // The initial coefficient vector is C0
  vnl_vector<double> C0 = xSolver->GetCoefficientArray();

  // Create mixed variations
  for(size_t i = 0; i < nRandVar; i++)
    {    
    vnl_vector<double> dx(nParams);
    for(size_t j = 0; j < nParams; j++)
      dx[j] = rand() * 2.0 / RAND_MAX - 1.0;

    xVariationNames.push_back(" Random Variation ");
    xVariationalBasisPS.set_row(i, dx);
   }

  // Create the per-component variations
  size_t iComp;
  for( iComp = 0; iComp < nParams; iComp++)
    {
    ostringstream oss;
    oss << " Coefficient " << iComp;
    
    vnl_vector<double> xParamVariation(nParams, 0.0); xParamVariation[iComp] = 1.0;
    
    xVariationNames.push_back(oss.str());
    xVariationalBasisPS.set_row(nRandVar + iComp, xParamVariation);
    }


  // Compute the variations in coefficient space
  for(iVar = 0; iVar < nVar; iVar++)
    {
    // Transform the variation in parameters to a variation in coefficients
    xVariationalBasisCS.set_row(iVar,
      xMapping->ApplyJacobianInParameters(C0, P0, xVariationalBasisPS.get_row(iVar)));
    }

  // Pass in the variations to the solver
  xSolver->SetVariationalBasis(xVariationalBasisCS);

  // For each variation, allocate a derivative array
  MedialAtom * dAtoms = new MedialAtom[nAtoms];

  // Set the coefficients for the initial solution
  xSolver->SetCoefficientArray(xMapping->Apply(C0, P0));

  // Solve the medial PDE for all these variations
  tInitialSolver.Start();
  xSolver->ComputeAtoms();
  tInitialSolver.Stop();

  // Make three copies of the atom array
  MedialAtom *A0 = new MedialAtom[nAtoms];
  MedialAtom *A1 = new MedialAtom[nAtoms];
  MedialAtom *A2 = new MedialAtom[nAtoms];
  MedialAtom *A3 = new MedialAtom[nAtoms];
  MedialAtom *A4 = new MedialAtom[nAtoms];

  // Create a hint vector
  GenericMedialModel::Vec xHint = xSolver->GetHintArray();

  // Copy the current solution into the first array
  std::copy(xSolver->GetAtomArray(), xSolver->GetAtomArray() + nAtoms, A0);

  // Begin the gradient computation
  tAnalytic.Start();
  xSolver->BeginGradientComputation();
  tAnalytic.Stop();
  
  // Create curvature arrays
  size_t nb = xSolver->GetNumberOfBoundaryPoints();

  // Compute all the other derivatives
  for(iVar = 0; iVar < nVar; iVar++) 
    {
    cout << "----------------------------------------" << endl;
    cout << xVariationNames[iVar] << endl;
    cout << "----------------------------------------" << endl;

    // Compute the directional derivative
    tAnalytic.Start();
    xSolver->ComputeAtomVariationalDerivative(iVar, dAtoms);
    tAnalytic.Stop();

    // Get the current variation in parameter space
    vnl_vector<double> xVar = xVariationalBasisPS.get_row(iVar);

    // Too big of an epsilon can cause us to break the atoms, so we
    // play to stay in legal range
    double eps = master_eps;

    while(eps > 1.0e-9)
      {
      try 
        {
        // Solve for the forward difference
        xSolver->SetCoefficientArray(xMapping->Apply(C0, P0 + 1.5 * eps * xVar));
        tCentral.Start(); xSolver->ComputeAtoms(xHint.data_block()); tCentral.Stop();
        std::copy(xSolver->GetAtomArray(), xSolver->GetAtomArray() + nAtoms, A1);

        // Solve for the forward difference
        xSolver->SetCoefficientArray(xMapping->Apply(C0, P0 + 0.5 * eps * xVar));
        tCentral.Start(); xSolver->ComputeAtoms(xHint.data_block()); tCentral.Stop();
        std::copy(xSolver->GetAtomArray(), xSolver->GetAtomArray() + nAtoms, A2);

        // Solve for the backward difference
        xSolver->SetCoefficientArray(xMapping->Apply(C0, P0 - 0.5 * eps * xVar));
        tCentral.Start(); xSolver->ComputeAtoms(xHint.data_block()); tCentral.Stop();
        std::copy(xSolver->GetAtomArray(), xSolver->GetAtomArray() + nAtoms, A3);

        // Solve for the backward difference
        xSolver->SetCoefficientArray(xMapping->Apply(C0, P0 - 1.5 * eps * xVar));
        tCentral.Start(); xSolver->ComputeAtoms(xHint.data_block()); tCentral.Stop();
        std::copy(xSolver->GetAtomArray(), xSolver->GetAtomArray() + nAtoms, A4);

        // Revert to original values
        xSolver->SetCoefficientArray(xMapping->Apply(C0, P0));
        xSolver->ComputeAtoms(xHint.data_block());

        // Exit loop 
        break;
        }
      catch(MedialModelException &)
        {
        // Epsilon is too large
        eps = 0.5 * eps;

        // Print a message
        cout << "Scaling eps to " << eps << endl;
        }
      }

    // Define an accumulator for all atom-wise errors
    DerivativeTestQuantities dtq;

    // Compute atom-wise errors
    for(unsigned int i = 0; i < nAtoms; i++)
      {

      // Point to the atoms
      MedialAtom &a1 = A1[i];
      MedialAtom &a2 = A2[i];
      MedialAtom &a3 = A3[i];
      MedialAtom &a4 = A4[i];
      MedialAtom &a0 = dAtoms[i];

      // Look at such a simple thing as a.X
      dtq.Update("Atom's X", a0.X, a1.X, a2.X, a3.X, a4.X, eps);
      dtq.Update("Atom's Xu", a0.Xu, a1.Xu, a2.Xu, a3.Xu, a4.Xu, eps);
      dtq.Update("Atom's Xv", a0.Xv, a1.Xv, a2.Xv, a3.Xv, a4.Xv, eps);
      dtq.Update("Atom's Xuu", a0.Xuu, a1.Xuu, a2.Xuu, a3.Xuu, a4.Xuu, eps);
      dtq.Update("Atom's Xuv", a0.Xuv, a1.Xuv, a2.Xuv, a3.Xuv, a4.Xuv, eps);
      dtq.Update("Atom's Xvv", a0.Xvv, a1.Xvv, a2.Xvv, a3.Xvv, a4.Xvv, eps);
      dtq.Update("Atom's R", a0.R, a1.R, a2.R, a3.R, a4.R, eps);
      dtq.Update("Atom's Ru", a0.Ru, a1.Ru, a2.Ru, a3.Ru, a4.Ru, eps);
      dtq.Update("Atom's Rv", a0.Rv, a1.Rv, a2.Rv, a3.Rv, a4.Rv, eps);
      dtq.Update("Atom's Fu", a0.Fu, a1.Fu, a2.Fu, a3.Fu, a4.Fu, eps);
      dtq.Update("Atom's Fv", a0.Fv, a1.Fv, a2.Fv, a3.Fv, a4.Fv, eps);
      dtq.Update("Atom's N", a0.N, a1.N, a2.N, a3.N, a4.N, eps);

      // Look at the difference in the specific terms of metric tensor
      for(size_t u = 0; u < 2; u++) for(size_t v = 0; v < 2; v++)
        {
        dtq.Update("Metric Tensor (Covariant)", 
          a0.G.xCovariantTensor[u][v],
          a1.G.xCovariantTensor[u][v],
          a2.G.xCovariantTensor[u][v],
          a3.G.xCovariantTensor[u][v],
          a4.G.xCovariantTensor[u][v], eps);

        dtq.Update("Metric Tensor (Contra)", 
          a0.G.xContravariantTensor[u][v],
          a1.G.xContravariantTensor[u][v],
          a2.G.xContravariantTensor[u][v],
          a3.G.xContravariantTensor[u][v],
          a4.G.xContravariantTensor[u][v], eps);

        for(size_t w = 0; w < 2; w++)
          {
          dtq.Update("ChristoffelSymols(1)",
            a0.G.xChristoffelFirst[u][v][w],
            a1.G.xChristoffelFirst[u][v][w],
            a2.G.xChristoffelFirst[u][v][w],
            a3.G.xChristoffelFirst[u][v][w],
            a4.G.xChristoffelFirst[u][v][w], eps);

          dtq.Update("ChristoffelSymols(2)",
            a0.G.xChristoffelSecond[u][v][w],
            a1.G.xChristoffelSecond[u][v][w],
            a2.G.xChristoffelSecond[u][v][w],
            a3.G.xChristoffelSecond[u][v][w],
            a4.G.xChristoffelSecond[u][v][w], eps);
          }
        }

      // Take the largest difference in Phi
      dtq.Update("Phi", a0.F, a1.F, a2.F, a3.F, a4.F, eps);

      // Also look at the differences in grad R mag sqr.
      dtq.Update("Sqr. Grad. Mag. of R", 
        a0.xGradRMagSqr, 
        a1.xGradRMagSqr, 
        a2.xGradRMagSqr,
        a3.xGradRMagSqr, 
        a4.xGradRMagSqr, eps);

      // Look at the gradR difference
      dtq.Update("Grad R", a0.xGradR, a1.xGradR, a2.xGradR, a3.xGradR, a4.xGradR, eps);

      // Mean and Gaussian curvature
      dtq.Update("Medial mean curvature", 
        a0.xMeanCurv, 
        a1.xMeanCurv, 
        a2.xMeanCurv,
        a3.xMeanCurv, 
        a4.xMeanCurv, eps);
      dtq.Update("Medial Gaussian curvature", 
        a0.xGaussCurv, 
        a1.xGaussCurv, 
        a2.xGaussCurv,
        a3.xGaussCurv, 
        a4.xGaussCurv, eps);

      // Take the largest difference in boundary derivs
      // for(unsigned int k = 0; k < 2; k++)
      //   {
      //   dtq.Update("Boundary Node Positions", a0.xBnd[k].X, a1.xBnd[k].X, a2.xBnd[k].X, eps);
      //   dtq.Update("Boundary Node Normals", a0.xBnd[k].N, a1.xBnd[k].N, a2.xBnd[k].N, eps);
      //   }
      }

    // Now, test the effect on SolutionData objects
    MedialIterationContext *xGrid = xSolver->GetIterationContext();
    SolutionData sd0(xGrid, A0);
    SolutionData sd1(xGrid, A1);
    SolutionData sd2(xGrid, A2);
    SolutionData sd3(xGrid, A3);
    SolutionData sd4(xGrid, A4);
    PartialDerivativeSolutionData pdsd(&sd0, dAtoms);

    // Compute the boundary weights derivatives
    sd0.ComputeIntegrationWeights();
    sd1.ComputeIntegrationWeights();
    sd2.ComputeIntegrationWeights();
    sd3.ComputeIntegrationWeights();
    sd4.ComputeIntegrationWeights();
    pdsd.ComputeIntegrationWeights();

    for(size_t q = 0; q < xGrid->GetNumberOfAtoms(); q++)
      {
      dtq.Update("Weights (medial)",
        pdsd.xMedialWeights[q],
        sd1.xMedialWeights[q],
        sd2.xMedialWeights[q],
        sd3.xMedialWeights[q],
        sd4.xMedialWeights[q], eps);
      }

    // more stuff
    for(MedialBoundaryTriangleIterator tit(xGrid); !tit.IsAtEnd(); ++tit)
      {
      dtq.Update("Weights (check)",
        pdsd.xBoundaryTriangleArea[tit.GetIndex()],
        sd1.xBoundaryTriangleArea[tit.GetIndex()],
        sd2.xBoundaryTriangleArea[tit.GetIndex()],
        sd3.xBoundaryTriangleArea[tit.GetIndex()],
        sd4.xBoundaryTriangleArea[tit.GetIndex()], eps);
      }

    // Compute the area weight difference
    double xMaxBndAreaDiff = 0.0;
    for(MedialBoundaryPointIterator it(xGrid); !it.IsAtEnd(); ++it)
      {
      dtq.Update("Weights (boundary)",
        pdsd.xBoundaryWeights[it.GetIndex()],
        sd1.xBoundaryWeights[it.GetIndex()],
        sd2.xBoundaryWeights[it.GetIndex()],
        sd3.xBoundaryWeights[it.GetIndex()],
        sd4.xBoundaryWeights[it.GetIndex()], eps);

      dtq.Update("Boundary Node Positions",
        GetBoundaryPoint(it, dAtoms).X,
        GetBoundaryPoint(it, A1).X,
        GetBoundaryPoint(it, A2).X,
        GetBoundaryPoint(it, A3).X,
        GetBoundaryPoint(it, A4).X, eps);

      dtq.Update("Boundary Node Normals",
        GetBoundaryPoint(it, dAtoms).N,
        GetBoundaryPoint(it, A1).N,
        GetBoundaryPoint(it, A2).N,
        GetBoundaryPoint(it, A3).N,
        GetBoundaryPoint(it, A4).N, eps);

      dtq.Update("Boundary Mean Curv",
        GetBoundaryPoint(it, dAtoms).curv_mean,
        GetBoundaryPoint(it, A1).curv_mean,
        GetBoundaryPoint(it, A2).curv_mean,
        GetBoundaryPoint(it, A3).curv_mean,
        GetBoundaryPoint(it, A4).curv_mean, eps);

      dtq.Update("Boundary Gauss Curv",
        GetBoundaryPoint(it, dAtoms).curv_gauss,
        GetBoundaryPoint(it, A1).curv_gauss,
        GetBoundaryPoint(it, A2).curv_gauss,
        GetBoundaryPoint(it, A3).curv_gauss,
        GetBoundaryPoint(it, A4).curv_gauss, eps);

      /* UNCOMMENT LATER! 
      dtq.Update("Integration Weights",
        pdsd.xInteriorVolumeElement[it.GetIndex()],
        sd1.xInteriorVolumeElement[it.GetIndex()],
        sd2.xInteriorVolumeElement[it.GetIndex()], eps);
      */
      }

    // Print the report
    cout << "Maximal differences between numeric and analytic derivatives:" << endl;
    dtq.PrintReport();

    // Update the global maximum counters
    dtqGlobal.Update(dtq);
    }

  // Delete the atoms
  delete dAtoms;

  // Delete the atom arrays
  delete A0; delete A1; delete A2; delete A3; delete A4;

  // Finalize the gradient computation
  tAnalytic.Start();
  xSolver->EndGradientComputation();
  tAnalytic.Stop();
  
    
  // Print the final report
  cout << "========================================" << endl;
  cout << "             FINAL REPORT               " << endl;
  cout << "----------------------------------------" << endl;

  cout << "Maximal differences between numeric and analytic derivatives:" << endl;
  dtqGlobal.PrintReport();

  cout << "----------------------------------------" << endl;
  cout << "Numeric Gradient CPU Secs. :" << tCentral.Read() << endl;
  cout << "Analytic Gradient CPU Secs.:" << tAnalytic.Read() << endl; 
  cout << "Initial Solution CPU Secs.:" << tInitialSolver.Read() << endl; 
  cout << "----------------------------------------" << endl;

  // Check if all the errors are reasonable
  bool flagPass = dtqGlobal.Test(1e100, 1);

  // Describe what happened
  if(flagPass)
    cout << "TEST PASSED! (Relative errors < 1.0)" << endl;
  else
    cout << "TEST FAILED! (Relative errors > 1.0)" << endl;
  cout << "========================================" << endl;

  // Return value
  return flagPass ? 0 : -1;
}

int TestOptimizerGradientComputation(
  MedialOptimizationProblem &mop, 
  CoefficientMapping &xMapping,
  GenericMedialModel *xSolver,
  double xStepSize, const char *nm_term, const char *nm_map)
{
  // Set up timers
  CodeTimer tAnalytic, tNumeric;

  // Get the number of parameters that are being optimized over
  size_t nParams = xMapping.GetNumberOfParameters();

  // Create the initial parameter vector and the vector to hold the gradient
  vnl_vector<double> P(nParams, 0.0), xGradient(nParams, 0.0);

  // We do not want to only perform the test at P = 0 because that may not
  // detect problems for other P values. So instead, we can move along the
  // gradient a certain direction.
  // mop.ComputeGradient(P.data_block(), xGradient.data_block());
  // cout << "GRADIENT AT 0 " << xGradient << endl;
  // P += xStepSize * xGradient;

  // Now we are ready to perform the test.
  // mop.PrintReport(cout);

  // Compute the analytic gradient
  tAnalytic.Start();
  mop.ComputeGradient(P.data_block(), xGradient.data_block());
  tAnalytic.Stop();
  
  // Keep track of the maximum error value
  double xMaxError = 0.0, xMaxRelError = 0.0;
  double fVal = mop.Evaluate(P.data_block());
  
  // Compute the numeric gradient of the objective function
  double eps = 0.0001;
  for(size_t i = 0; i < nParams; i++)
    {
    // Get the current value of the coefficient
    double ci = P[i];

    // Compute the objective function at point 1
    P[i] = ci + eps;
    tNumeric.Start();
    double f1 = mop.Evaluate(P.data_block());
    tNumeric.Stop();

    // Compute the objective at point 2
    P[i] = ci - eps;
    tNumeric.Start();
    double f2 = mop.Evaluate(P.data_block());
    tNumeric.Stop();

    // Restore the coefficient
    P[i] = ci;
    
    // Print the derivative
    double dNumeric = 0.5 * (f1 - f2) / eps;
    double xDifference = fabs(xGradient[i] - dNumeric);
    double xRelError = xDifference / (0.5 * (fabs(dNumeric) + fabs(xGradient[i])) + eps);
    // printf("D[x_%04d](f):  AN = %+6E   CD = %+6E   AE = %+6E   RE = %+6E\n",
    //   i, dNumeric, xGradient[i], xDifference, xRelError);

    // Update the max error tracker
    UpdateMax(xDifference, xMaxError);
    UpdateMax(xRelError, xMaxRelError);
    }

  // Report difference in time
  printf("%12s  %12s :", nm_term, nm_map);
  printf("   %+4.2le   %+4.2le   %+4.2le   %7.2f   %7.2f  %s\n",
    fVal, xMaxError, xMaxRelError, tAnalytic.Read(), tNumeric.Read(),
    (xMaxRelError > 1 ? -1 : 0) ? "FAIL" : "pass");

  return xMaxError > eps ? -1 : 0;
}
