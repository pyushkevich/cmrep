#include "ScriptInterface.h"
#include "BasisFunctions2D.h"
#include "MedialAtom.h"
#include "CartesianMedialModel.h"
#include "OptimizationTerms.h"
#include "DiffeomorphicEnergyTerm.h"
#include "CoefficientMapping.h"
#include "MedialAtomGrid.h"
#include "PrincipalComponents.h"
#include "ITKImageWrapper.h"
#include "SmoothedImageSampler.h"
#include "itkImage.h"
#include "TestSolver.h"
#include "vnl/vnl_erf.h"
#include "vnl/vnl_random.h"

#include "vtkOBJReader.h"
#include "vtkBYUWriter.h"
#include "vtkPolyData.h"

#include <string>
#include <iostream>

using namespace std;
using namespace medialpde;

void TestTetGen(GenericMedialModel *model);

string dirWork = "/home/pauly/data2005/Stanley/data/";
// string dirWork = "/mnt/data2/PUBLIC/Data/Input/StanleySchizophrenia/";

/**
 * This is a test class that is a sphere that pretends to be a floating
 * point image. The idea is to use this ellipsoid for testing image-based
 * derivative code, as the derivatives computed here are exact rather than
 * numerical.
 */
class TestFloatImage : public FloatImage
{
public:
  TestFloatImage(const SMLVec3d &C, double R, double sigma)
    {
    this->C = C;
    this->R = R;
    this->sigma = sigma;
    } 

  // Interpolate the image at a given position
  float Interpolate(const SMLVec3d &x)
    { 
    double u = ( R - (x - C).two_norm() ) / sigma;
    return vnl_erf(u);
    }
    
  // Interpolate the image gradient
  void InterpolateImageGradient(const SMLVec3d &x, SMLVec3f &g)
    {
    SMLVec3d g1;
    InterpolateImageGradient(x, g1);
    g[0] = (float) g1[0];
    g[1] = (float) g1[1];
    g[2] = (float) g1[2];
    }

  // Interpolate the image gradient
  void InterpolateImageGradient(const SMLVec3d &x, SMLVec3d &g)
    {
    double z = (x - C).two_norm();
    double u = ( R - z ) / sigma;
    double a = 2.0 * exp(-u * u) / sqrt(M_PI);
    
    g[0] = a * (- (x[0] - C[0]) / z) / sigma;
    g[1] = a * (- (x[1] - C[1]) / z) / sigma;
    g[2] = a * (- (x[2] - C[2]) / z) / sigma;
    }

  // Integrate positive voxels in the image volume
  double IntegratePositiveVoxels()
    {
    // This is a hack!
    return 4.0 * M_PI * R * R * R / 3.0;
    }

private:
  SMLVec3d C;
  double R, sigma;
};



class TestFunction01 : public EuclideanFunction
{
public:
  double Evaluate(const SMLVec3d &x)
    { return 1; }
};

void ExportMedialMeshToVTK(
  GenericMedialModel *xModel, ITKImageWrapper<float> *xImage, const char *file);

void ExportBoundaryMeshToVTK(
  GenericMedialModel *xModel, ITKImageWrapper<float> *xImage, const char *file);

/**
 * This test ensures that all the volume elements in a model are computed
 * correctly.
 */
int TestVolumeComputation(const char *file)
{
  // The number of cuts at which to test
  size_t i,j,k,nCuts = 5;
 
  // Success code
  bool flagSuccess = true;

  // Read the model
  MedialPDE mp(file);
  GenericMedialModel *model = mp.GetMedialModel();

  size_t nInternalPoints = model->GetNumberOfInternalPoints(nCuts);

  // Use wedge-based system to compute volume. This is an old way that is
  // fairly slow but seemingly reliable
  double *xVolWedgeWeights = new double[model->GetNumberOfInternalPoints(nCuts)];
  SMLVec3d *xInternalPoints = new SMLVec3d[model->GetNumberOfInternalPoints(nCuts)];
  ComputeMedialInternalPoints(
    model->GetIterationContext(),
    model->GetAtomArray(),
    nCuts,
    xInternalPoints);
  double volWedgeBased = 
    ComputeMedialInternalVolumeWeights(
      model->GetIterationContext(), 
      model->GetAtomArray(),
      xInternalPoints,
      nCuts, 
      xVolWedgeWeights, NULL);

  // Use the newer method
  double *xMedSurfaceAreaWeights = new double[model->GetNumberOfAtoms()];
  double *xBndSurfaceAreaWeights = new double[model->GetNumberOfBoundaryPoints()];
  double *xFastVolumeWeights = new double[model->GetNumberOfInternalPoints(nCuts)];
  SMLVec3d *xMedNormals = new SMLVec3d[model->GetNumberOfAtoms()]; 
  SMLVec3d *xBndNormals = new SMLVec3d[model->GetNumberOfBoundaryPoints()];
  double xBndSurfaceArea = 
    ComputeMedialBoundaryAreaWeights(
      model->GetIterationContext(), 
      model->GetAtomArray(),
      xBndNormals,
      xBndSurfaceAreaWeights);
  double xMedSurfaceArea = 
    ComputeMedialSurfaceAreaWeights(
      model->GetIterationContext(), 
      model->GetAtomArray(),
      xMedNormals,
      xMedSurfaceAreaWeights);
  double volFastMethod = 
    FastComputeMedialInternalVolumeWeights(
      model->GetIterationContext(), 
      model->GetAtomArray(),
      xMedSurfaceAreaWeights, 
      xBndSurfaceAreaWeights,
      nCuts, xFastVolumeWeights);

  // Report what's computed by the two methods
  cout << "Wedge-Based Volume: " << volWedgeBased << endl;
  cout << "Fast Method Volume: " << volFastMethod << endl;

  // Compute the volume via VolumeOverlap class
  TestFloatImage testimg(model->GetCenterOfRotation(), 100.0, 1.0);
  VolumeOverlapEnergyTerm voet(model, &testimg, nCuts);
  SolutionData sdat(model->GetIterationContext(), model->GetAtomArray());
  sdat.ComputeIntegrationWeights();
  voet.ComputeEnergy(&sdat);
  cout << "Vol Ovl ETerm Vol : " << voet.GetModelVolume() << endl;

  // Compare the volumes
  StatisticsAccumulator volRatio, volFast, volWedge;
  for(size_t i = 0; i < model->GetNumberOfInternalPoints(nCuts); i++)
    {
    volFast.Update(xFastVolumeWeights[i]);
    volWedge.Update(xVolWedgeWeights[i]);
    volRatio.Update(xFastVolumeWeights[i] / xVolWedgeWeights[i]);
    }
  cout << "Mean F/W Ratio    : " << volRatio.GetMean() << endl;
  cout << "SDev F/W Ratio    : " << volRatio.GetStdDev() << endl;
  cout << "Min/Max F/W Ratio : " << volRatio.GetMin()  
    << ", " << volRatio.GetMax() << endl;

  // Create a solution data based on this model
  MedialIterationContext *grid = model->GetIterationContext();

  // Test 1. Ensure that the volume of every wedge in the model is positive
  cout << "Max wedge velement: " << volWedge.GetMax() << endl;
  cout << "Min wedge velement: " << volWedge.GetMin() << endl;
  cout << "Avg wedge velement: " << volWedge.GetMean() << endl;

  if(volWedge.GetMin() < 0.0)
    flagSuccess = false;

  // Test 2. Ensure that the boundary area elements are also positive
  StatisticsAccumulator saArea;
  for(i = 0; i < model->GetNumberOfBoundaryPoints(); i++)
    saArea.Update(xBndSurfaceAreaWeights[i]);

  cout << "Max area element: " << saArea.GetMax() << endl;
  cout << "Min area element: " << saArea.GetMin() << endl;
  cout << "Avg area element: " << saArea.GetMean() << endl;
  cout << "Sum area element: " << saArea.GetSum() << endl;
  cout << "Surface area: " << xBndSurfaceArea << endl;

  if(saArea.GetMin() < 0.0 || saArea.GetSum() != xBndSurfaceArea)
    flagSuccess = false;

  // Test 3. Compare volume computed by integration to volume estimated by
  // Green's theorem
  double xVolGreen = 0.0;
  for(MedialBoundaryPointIterator bit(grid); !bit.IsAtEnd(); ++bit)
    {
    BoundaryAtom bat = GetBoundaryPoint(bit, model->GetAtomArray());
    xVolGreen += dot_product(bat.X, bat.N) * xBndSurfaceAreaWeights[bit.GetIndex()];
    }
  xVolGreen /= 3.0;

  cout << "Green Thm. Integral : " << xVolGreen << endl;

  // Test 4. Another way to compute Green's integral
  double xVolGreen2 = 0.0;
  for(MedialBoundaryTriangleIterator btt(grid); !btt.IsAtEnd(); ++btt)
    {
    // Get the triangle
    SMLVec3d A = GetBoundaryPoint(btt, model->GetAtomArray(), 0).X;
    SMLVec3d B = GetBoundaryPoint(btt, model->GetAtomArray(), 1).X;
    SMLVec3d C = GetBoundaryPoint(btt, model->GetAtomArray(), 2).X;

    // Get the area-scaled normal vector
    SMLVec3d N = vnl_cross_3d(B-A,C-A) / 2.0;
    SMLVec3d X = (A + B + C) / 3.0;
    xVolGreen2 += dot_product(X, N) / 3.0;
    }
  cout << "Green Thm. Integral 2 : " << xVolGreen2 << endl;

  // Export model for external verification
  ExportBoundaryMeshToVTK(model, NULL, "volume_test.vtk");

  // Return the test score
  return flagSuccess ? 0 : -1;
}

void Test01()
{
  CartesianMPDE *mp = new CartesianMPDE(3, 5, 20, 40);
  //mp->LoadFromParameterFile(fMrep.c_str());
  
  mp->LoadFromDiscreteMRep("/tmp/surf01.txt",-0.3);
  // mp->GenerateSampleModel();
  mp->Solve();
  mp->SaveBYUMesh("temp.byu");

  // Make sure areas and volumes add up
  
  // Load the image and gradients
  // FloatImage img;
  // img.LoadFromFile((dirWork + "avg/average_hippo_blurred_hi.mha").c_str());

  // Match the volume to the image
  // mp->MatchImageByMoments(&img, 5);

  
  
  // RenderMedialPDE(mp);
}

/**
void Test02()
{
  // Decide how many points to interpolate
  unsigned int nSide = 11, nPoints = nSide * nSide;
  double uStep = 1.0 / (nSide - 1);

  // Allocate arrays of points and coordinates
  double xPoints[nPoints], yPoints[nPoints], zPoints[nPoints];
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
      ++i;
      }

  // Peform the fit
  FourierSurfaceOld s1(5, 5, 3);
  s1.FitData(0, nPoints, uPoints, 1, vPoints, 1, xPoints, 1);
  s1.FitData(1, nPoints, uPoints, 1, vPoints, 1, yPoints, 1);
  s1.FitData(2, nPoints, uPoints, 1, vPoints, 1, zPoints, 1);

  FourierSurface s2(5, 5);
  s2.FitToData(nPoints, 0, uPoints, vPoints, xPoints);
  s2.FitToData(nPoints, 1, uPoints, vPoints, yPoints);
  s2.FitToData(nPoints, 2, uPoints, vPoints, zPoints);

  SMLVec3d T1, T2;
  for(double u = 0.0; u <= 1.0; u+=0.234)
    for(double v = 0.0; v <= 1.0; v+=0.1234)
      {
      T1[0] = s1.Evaluate(u, v, 1, 0, 0);
      T1[1] = s1.Evaluate(u, v, 1, 0, 1);
      T1[2] = s1.Evaluate(u, v, 1, 0, 2);
      s2.EvaluateDerivative(u, v, 1, 0, 0, 3, T2.data_block());
      cout << T1 - T2 << endl;
      }
  
}
*/

void Test03()
{
  // Create a medial PDE object
  CartesianMPDE *mp = new CartesianMPDE(5, 5, 50, 100);
  mp->LoadFromParameterFile((dirWork + "avg/average_mrepL_01.mpde").c_str());
  mp->Solve();

  // Convert the binary image into a floating point and compute match
  FloatImage *img = new FloatImage();
  img->LoadFromFile((dirWork + "avg/average_hippo_blurred.mha").c_str());
  img->LoadGradientFromFile(0, (dirWork + "avg/average_hippo_blurred.dx.mha").c_str());
  img->LoadGradientFromFile(1, (dirWork + "avg/average_hippo_blurred.dy.mha").c_str());
  img->LoadGradientFromFile(2, (dirWork + "avg/average_hippo_blurred.dz.mha").c_str());
  img->SetOutsideValue(-1.0);

  // Compute the gradient of the match
  cout << "Image Match = " << mp->ComputeImageMatch(img) << endl;

  // TODO: Load optimization parameters from a file
  mp->RunOptimization(img, 30, "file.txt");
  mp->SaveToParameterFile((dirWork + "avg/average_mrepL_02.mpde").c_str());
  mp->SaveBYUMesh((dirWork + "avg/average_mrepL_02.byu").c_str());
}

int TestWedgeVolume()
{
  // Define an object with a volume that we know
  double vol1 = WedgeVolume(
    SMLVec3d(0, 0, 0), SMLVec3d(1, 0, 0), SMLVec3d(1, 1, 0),
    SMLVec3d(0, 0, 1), SMLVec3d(1, 0, 1), SMLVec3d(1, 1, 1));

  cout << "Volume 1 = " << vol1 << " and should be 0.5" << endl;
  
  // Compare to the right amount
  return (vol1 == 0.5) ? 0 : -1;
}

void TestCellVolume()
{
  // Compute the volume of a unit cube
  cout << "Unit Cube Volume: " << CellVolume(
    SMLVec3d(0,0,0),SMLVec3d(0,0,1),SMLVec3d(0,1,0),SMLVec3d(0,1,1),
    SMLVec3d(1,0,0),SMLVec3d(1,0,1),SMLVec3d(1,1,0),SMLVec3d(1,1,1)) << endl;

  // Compute the volume of an elongated cube
  cout << "123 Cube Volume : " << CellVolume(
    SMLVec3d(0,0,0),SMLVec3d(0,0,3),SMLVec3d(0,2,0),SMLVec3d(0,2,3),
    SMLVec3d(1,0,0),SMLVec3d(1,0,3),SMLVec3d(1,2,0),SMLVec3d(1,2,3)) << endl;
}


/** Test differential geometry relationships */
int TestDifferentialGeometry(const char *fnMPDE)
{
  // Load a medial PDE for the test
  CartesianMPDE mp(2, 4, 25, 49);
  mp.LoadFromParameterFile(fnMPDE);

  // Pick a point to evaluate at
  double u = 0.6, v = 0.2, eps = 0.0001;

  // Compute the jet at this point
  SMLVec3d X00, X10, X01, X20, X11, X02;
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 0, 0, 0, 3, X00.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 1, 0, 0, 3, X10.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 0, 1, 0, 3, X01.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 2, 0, 0, 3, X20.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 1, 1, 0, 3, X11.data_block());
  mp.GetMedialSurface()->EvaluateDerivative(u, v, 0, 2, 0, 3, X02.data_block());

  // Compute the differential geometry
  GeometryDescriptor gd;
  gd.SetTwoJet(
    X00.data_block(), X10.data_block(), X01.data_block(),
    X20.data_block(), X11.data_block(), X02.data_block());

  // Check the symbols
  double E = dot_product(X10, X10);
  double F = dot_product(X10, X01);
  double G = dot_product(X01, X01);
  double Eu = 2.0 * dot_product(X10, X20);
  double Ev = 2.0 * dot_product(X10, X11);
  double Gu = 2.0 * dot_product(X01, X11);
  double Gv = 2.0 * dot_product(X01, X02);
  double Fu = dot_product(X20, X01) + dot_product(X11, X10);
  double Fv = dot_product(X11, X01) + dot_product(X02, X10);

  double g = E * G - F * F;

  double G111 = (G * Eu - 2 * F * Fu + F * Ev) / (2 * g);
  double G121 = (G * Ev - F * Gu) / (2 * g);
  double G221 = (2 * G * Fv - G * Gu - F * Gv) / (2 * g);
  double G112 = (2 * E * Fu - E * Ev - F * Eu) / (2 * g);
  double G122 = (E * Gu - F * Ev) / (2 * g);
  double G222 = (E * Gv - 2 * F * Fv + F * Gu) / (2 * g);

  cout << "E = " << E << " vs " << gd.xCovariantTensor[0][0] << endl;
  cout << "F = " << F << " vs " << gd.xCovariantTensor[0][1] << endl;
  cout << "G = " << G << " vs " << gd.xCovariantTensor[1][1] << endl;
  cout << "g = " << g << " vs " << gd.g << endl;

  cout << "G 11 1 = " << G111 << " vs " << gd.xChristoffelSecond[0][0][0] << endl;
  cout << "G 12 1 = " << G121 << " vs " << gd.xChristoffelSecond[0][1][0] << endl;
  cout << "G 22 1 = " << G221 << " vs " << gd.xChristoffelSecond[1][1][0] << endl;
  cout << "G 11 2 = " << G112 << " vs " << gd.xChristoffelSecond[0][0][1] << endl;
  cout << "G 12 2 = " << G122 << " vs " << gd.xChristoffelSecond[0][1][1] << endl;
  cout << "G 22 2 = " << G222 << " vs " << gd.xChristoffelSecond[1][1][1] << endl;

  gd.PrintSelf(cout);
  
  // Evaluate the nearby points
  /*
  SMLVec3d YPZ, YPP, YZP, YMP, YMZ, YMM, YZM, YPM;
  mp.GetSurface()->EvaluateDerivative(u + eps, v      , 0, 0, 0, 3, YPZ.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u + eps, v + eps, 0, 0, 0, 3, YPP.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u      , v + eps, 0, 0, 0, 3, YZP.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u - eps, v + eps, 0, 0, 0, 3, YMP.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u - eps, v      , 0, 0, 0, 3, YMZ.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u - eps, v - eps, 0, 0, 0, 3, YMM.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u      , v - eps, 0, 0, 0, 3, YZM.data_block()); 
  mp.GetSurface()->EvaluateDerivative(u + eps, v - eps, 0, 0, 0, 3, YPM.data_block()); 

  // Compute the partial derivatives
  SMLVec3d Y00, Y01, Y10, Y20, Y11, Y02;
  Y10 = (YPZ - YMZ) / (2.0 * eps);
  Y01 = (YZP - YZM) / (2.0 * eps);
  Y20 = (YPZ + YMZ - X00) / (eps * eps);
  Y02 = (YZP + YZM - X00) / (eps * eps);
  Y11 = (YPP + YMM - YPM - YMP) / (4 * eps * eps);

  // Compute the partial derivative tensor
  */
  return 0;
}

/**
 * This test routine checks the accuracy of basis function variation
 * computations, under various masks
 */
int TestBasisFunctionVariation(const char *fnMPDE)
{
  int iReturn = 0;

  // Load a medial PDE for the test
  CartesianMPDE mp(fnMPDE);

  // Get the model
  GenericMedialModel *model = mp.GetMedialModel();

  // Set up a selection mask
  vnl_vector<size_t> xMask(model->GetNumberOfCoefficients(), 0);
  vnl_random rnd;
  for(size_t i = 0; i < xMask.size(); i+=4)
    xMask.update(vnl_vector<size_t>(4, rnd.lrand32(0, 1)), i);

  // Set up an array of coefficient masks
  vector<CoefficientMapping *> vm;
  vm.push_back(new IdentityCoefficientMapping(model->GetNumberOfCoefficients()));
  vm.push_back(new AffineTransformCoefficientMapping(model));
  vm.push_back(new SubsetCoefficientMapping(xMask));
  
  char *nm[] = {
    "IdentityCoefficientMapping",
    "AffineTransformCoefficientMapping",
    "SubsetCoefficientMapping" };
  
  // Repeat for each mask
  for(size_t i=0; i<vm.size(); i++)
    {
    // Get the current value of the coefficients and parameters
    vnl_vector<double> C, C0 = model->GetCoefficientArray();
    vnl_vector<double> P(vm[i]->GetNumberOfParameters(), 0.0);

    // Create a random variation vector in P
    double eps = 0.0001;
    vnl_vector<double> dP(vm[i]->GetNumberOfParameters());
    for(size_t j=0; j < dP.size(); j++)
      dP[j] = rand() * 2.0 / RAND_MAX - 1.0;

    // Compute the surface derivative using each variation
    SMLVec4d f1, f2, dfNum, dfAn;

    // Evaluate the surface at finite differences
    model->SetCoefficientArray(vm[i]->Apply(C0, P + eps * dP));
    model->ComputeAtoms();
    mp.GetMedialSurface()->Evaluate(0.5, 0.5, f1.data_block());

    model->SetCoefficientArray(vm[i]->Apply(C0, P - eps * dP));
    model->ComputeAtoms();
    mp.GetMedialSurface()->Evaluate(0.5, 0.5, f2.data_block());

    dfNum = (f1 - f2) * 0.5 / eps;

    // Same evaluation using analytical derivative
    model->SetCoefficientArray(vm[i]->Apply(C0, P));
    model->ComputeAtoms();

    // Get the variation corresponding to dP
    IHyperSurface2D *xVariation = mp.GetMedialSurface()->GetVariationSurface(
      vm[i]->ApplyJacobianInParameters(C0, P, dP).data_block());
    xVariation->Evaluate(0.5, 0.5, dfAn.data_block());

    // Report
    double dMax = (dfAn - dfNum).inf_norm();
    if(dMax > eps)
      iReturn += -1;
    
    cout << "Testing mask " << nm[i] << endl;
    cout << "  right f.d.    : " << f1 << endl;
    cout << "  left f.d.     : " << f2 << endl;
    cout << "  central diff. : " << dfNum << endl;
    cout << "  numeric der.  : " << dfAn << endl;
    cout << "  maximum error : " << dMax << endl;
    }
 
  return iReturn;
}

/**
 * This test routine makes sure that the derivatives computed in the MedialPDE
 * solver and related classes are correct, bt comparing them to central
 * difference derivatives
 *
int TestDerivativesNoImage(const char *fnMPDE, const char *fnPCA)
{
  int iReturn = 0;

  // Create a cartesian medial model for testing
  CartesianMPDE mp(fnMPDE);

  // Extract the medial model itself
  GenericMedialModel *model = mp.GetMedialModel();

  // Test the most convoluted mask that we have
  cout << "**************************************************" << endl;
  cout << "** TESTIING Affine3DAndPCACoefficientMask       **" << endl;
  cout << "**************************************************" << endl;
 
  // We have to set the number of coefficients to match the PCA matrix
  mp.SetNumberOfCoefficients(8, 10);
  
  // Create a PCA based mask based on the test matrix
  double xAPCAPosn[] = 
    { 1.2, 0.1, 0.2,           // Affine Matrix
     -0.2, 0.9, 0.1, 
     -0.1, 0.1, 1.0, 
      5.0, 3.0,-2.0,           // Translation Vector
      0.5, 0.4, 0.3,-0.2, 0.3  // PCA offsets
    };

  // Read principal components from the file
  vnl_matrix<double> pcaMatrix;
  ReadMatrixFile(pcaMatrix, fnPCA);
  PrincipalComponents pca(pcaMatrix);
  cout << "PCA Matrix Size: " << pcaMatrix.rows() << " x " << pcaMatrix.columns() << endl;

  // Create the PCA/Affine optimizer with 5 modes
  PCAPlusAffineCoefficientMapping xPCAMask(model, &pca, 5);
  // PCACoefficientMapping xPCAMask(&pca, 5);

  // Create the starting point vector
  vnl_vector<double> pPCA(xAPCAPosn, 17);

  // Perform the test
  // iReturn += TestGradientComputation(model, &xPCAMask, pPCA, 3);
  iReturn += TestGradientComputation(model, &xPCAMask, 3);

  // We have to set the number of coefficients to match the PCA matrix
  mp.SetNumberOfCoefficients(2, 4);
  
  // Test straight-through gradient computation
  cout << "**************************************************" << endl;
  cout << "** TESTIING IdentityCoefficientMapping          **" << endl;
  cout << "**************************************************" << endl;
 
  IdentityCoefficientMapping xMapping(model);
  iReturn += TestGradientComputation(model, &xMapping, 3);

  // Test affine transform gradient computation
  cout << "**************************************************" << endl;
  cout << "** TESTIING AffineTransformCoefficientMapping   **" << endl;
  cout << "**************************************************" << endl;

  // Don't start at the default position, or we may get a false positive
  double xAffinePosn[] = 
    { 1.2, 0.2, 0.1, 
      0.1, 0.9,-0.2,
      0.3,-0.2, 1.4,
      7.0,-4.5, 2.3 };
  
  AffineTransformCoefficientMapping xAffineMask(model);
  vnl_vector<double> pAffine(xAffinePosn, 12);
  iReturn += TestGradientComputation(mp.GetMedialModel(), &xAffineMask, pAffine, 3);

  // Return the total of the test values
  return iReturn;
}*/

/**
 * This test routine makes sure that the derivatives computed in the MedialPDE
 * solver and related classes are correct, bt comparing them to central
 * difference derivatives
 */
int TestDerivativesNoImage(const char *fnMPDE)
{
  int iReturn = 0;

  // Create a cartesian medial model for testing
  MedialPDE mp(fnMPDE);

  // Extract the medial model itself
  GenericMedialModel *model = mp.GetMedialModel();

  // Test straight-through gradient computation
  cout << "**************************************************" << endl;
  cout << "** TESTIING IdentityCoefficientMapping          **" << endl;
  cout << "**************************************************" << endl;
 
  IdentityCoefficientMapping xMapping(model);
  iReturn += TestGradientComputation(model, &xMapping, 3);

  // Test affine transform gradient computation
  cout << "**************************************************" << endl;
  cout << "** TESTIING AffineTransformCoefficientMapping   **" << endl;
  cout << "**************************************************" << endl;

  // Don't start at the default position, or we may get a false positive
  double xAffinePosn[] = 
    { 1.2, 0.2, 0.1, 
      0.1, 0.9,-0.2,
      0.3,-0.2, 1.4,
      7.0,-4.5, 2.3 };
  
  AffineTransformCoefficientMapping xAffineMask(model);
  vnl_vector<double> pAffine(xAffinePosn, 12);
  iReturn += TestGradientComputation(mp.GetMedialModel(), &xAffineMask, pAffine, 3);

  // Test the most convoluted mask that we have
  cout << "**************************************************" << endl;
  cout << "** TESTIING Affine3DAndPCACoefficientMask       **" << endl;
  cout << "**************************************************" << endl;
 
  // Create a PCA based mask based on the test matrix
  double xAPCAPosn[] = 
    { 0.2, 0.1, 0.2,           // Affine Matrix
     -0.2,-0.1, 0.1, 
     -0.1, 0.1, 0.0, 
      5.0, 3.0,-2.0,           // Translation Vector
      0.5, 0.4, 0.3,-0.2, 0.3  // PCA offsets
    };

  // Create a PCA based on a random matrix
  vnl_matrix<double> pcaMatrix(15, model->GetNumberOfCoefficients());
  vnl_random randy;
  for(size_t q = 0; q < pcaMatrix.rows(); q++) 
    for(size_t k = 0; k < pcaMatrix.columns(); k++)
      pcaMatrix[q][k] = randy.drand32(-0.01, 0.01);
  PrincipalComponents pca(pcaMatrix);

  // Create the PCA/Affine optimizer with 5 modes
  PCAPlusAffineCoefficientMapping xPCAMask(model, &pca, 5);

  // Create the starting point vector
  vnl_vector<double> pPCA(xAPCAPosn, 17);

  // Perform the test
  iReturn += TestGradientComputation(model, &xPCAMask, pPCA, 3);

  // Return the total of the test values
  return iReturn;
}

int TestGradientTiming(const char *fnMPDE)
{
  // Create and read the MPDE
  CartesianMPDE mp(2, 4, 32, 80);
  mp.LoadFromParameterFile(fnMPDE);

  // Test straight-through gradient computation
  IdentityCoefficientMapping xMask(mp.GetMedialModel());
  return TestGradientComputation(mp.GetMedialModel(), &xMask);
}

void MakeFlatTemplate(FourierSurface *xSurface)
{
  size_t n = 10;
  size_t q = (n+1) * (n+1);
  double *u = new double[q], *v = new double[q];
  double *x = new double[q], *y = new double[q];
  double *z = new double[q], *rho = new double[q];

  size_t k = 0;
  for(size_t i = 0; i <= n; i++) for(size_t j = 0; j <= n; j++)
    {
    u[k] = i * 1.0 / n;
    v[k] = j * 1.0 / n;
    x[k] = 6 * u[k];
    y[k] = 12 * v[k];
    z[k] = 0.0;
    rho[k] = -0.45;
    k++;
    }

  xSurface->FitToData(q, 0, u, v, x);
  xSurface->FitToData(q, 1, u, v, y);
  xSurface->FitToData(q, 2, u, v, z);
  xSurface->FitToData(q, 3, u, v, rho);

  delete u; delete v; delete x; delete y; delete z; delete rho;
}

struct TermInfo
{
  string name;
  EnergyTerm *et;
  double step;
  TermInfo(const char *nm, EnergyTerm *e, double s)
    {
    name = nm;
    et = e;
    step = s;
    }
};

int TestDerivativesWithImage(const char *fnMPDE, FloatImage *img, const char *fnRefMPDE = NULL)
{
  // Return Code
  int iReturn = 0;

  // Load a cartesian medial PDE
  MedialPDE mp(fnMPDE);
  
  // Define a test image (this in not a real thing)
  GenericMedialModel *model = mp.GetMedialModel();
  model->ComputeAtoms();

  SMLVec3d C = model->GetCenterOfRotation();
  double rLogSum = 0;
  for(size_t ia = 0; ia < model->GetNumberOfAtoms(); ia++)
    rLogSum += log((C - model->GetAtomArray()[ia].X).magnitude());
  double rMean = exp(rLogSum / model->GetNumberOfAtoms());

  TestFloatImage testimg(C, rMean, rMean/10);

  if(img == NULL)
    img = &testimg;

  // Define a test point set and a test radius function for computation
  vnl_vector<double> px(model->GetNumberOfAtoms(), 0.0);
  vnl_vector<double> py(model->GetNumberOfAtoms(), 0.0);
  vnl_vector<double> pz(model->GetNumberOfAtoms(), 0.0);
  vnl_vector<double> rad(model->GetNumberOfAtoms(), 0.0);
  vnl_random randy;
  for(size_t i = 0; i < rad.size(); i++)
    {
    px[i] = model->GetAtomArray()[i].X[0] + randy.drand32(-0.1, 0.1);
    py[i] = model->GetAtomArray()[i].X[1] + randy.drand32(-0.1, 0.1);
    pz[i] = model->GetAtomArray()[i].X[2] + randy.drand32(-0.1, 0.1);
    rad[i] = model->GetAtomArray()[i].R * randy.drand32(0.9, 1.1);
    }

  // If the reference model is supplied, load it
  MedialPDE *ref = NULL;
  if(fnRefMPDE)
    ref = new MedialPDE(fnRefMPDE);

  // Create an array of image match terms
  vector<TermInfo> vt;
  if(ref)
    {
    vector<GenericMedialModel *> refvec;
    refvec.push_back(ref->GetMedialModel());
    vt.push_back(TermInfo(
        "", new LocalDistanceDifferenceEnergyTerm(model, refvec), 0.1));
    }
  vt.push_back(TermInfo(
      "", new DiffeomorphicEnergyTerm(model), 0.1));
  vt.push_back(TermInfo(
      "", new RadiusPenaltyTerm(0.01, 4, 100, 10), 0.1));
  vt.push_back(TermInfo(
      "", new BoundaryImageMatchTerm(model, img), 0.1));
  vt.push_back(TermInfo(
      "", new BoundaryJacobianEnergyTerm(), 1.0e-4));
  vt.push_back(TermInfo(
      "", new BoundaryCurvaturePenalty(model), 0.1));
  vt.push_back(TermInfo(
      "", new MedialCurvaturePenalty(), 0.1));
  vt.push_back(TermInfo(
      "", new VolumeOverlapEnergyTerm(model, img, 4), 0.1));
  vt.push_back(TermInfo(
      "", new MedialBendingEnergyTerm(model), 0.1));
  vt.push_back(TermInfo(
      "", new DistanceToPointSetEnergyTerm(
        model, px.data_block(), py.data_block(), pz.data_block()), 0.1));
  vt.push_back(TermInfo(
      "", new DistanceToRadiusFieldEnergyTerm(
        model, rad.data_block()), 0.1));
  vt.push_back(TermInfo(
      "", new MedialAnglesPenaltyTerm(model), 1.0e-4));
  vt.push_back(TermInfo(
      "", new MedialRegularityTerm(model), 0.1));
  vt.push_back(TermInfo(
      "", new BoundaryGradRPenaltyTerm(), 0.1));

  // Create an array of masks
  vector<CoefficientMapping *> vm;
  vm.push_back(new IdentityCoefficientMapping(model));
  vm.push_back(new SubsetCoefficientMapping(model->GetRadialCoefficientMask()));
  vm.push_back(new AffineTransformCoefficientMapping(model));

  char *nm[] = {
    "IDENTY", "RADIAL", "AFFINE" };

  printf("%12s  %12s :  %10s  %10s  %10s  %8s  %8s  P/F\n",
    "ENERGY-TERM ", "COEF-MAP",
    "ENERGY-VAL", 
    "ABS-ERROR", "REL-ERROR", "T(AGRAD)", "T(NGRAD)");

  // Loop over both options
  size_t i, j;
  for(i = 0; i < vt.size(); i++) for(j = 0; j < vm.size(); j++)
    {
    
    MedialPDE mp1(fnMPDE);
  
    // Set up the test
    // MedialOptimizationProblem mop(model, vm[j]);
    MedialOptimizationProblem mop(mp1.GetMedialModel(), vm[j]);
    mop.AddEnergyTerm(vt[i].et, 1.0);
    iReturn += TestOptimizerGradientComputation(
      mop, *vm[j], mp1.GetMedialModel(), vt[i].step, vt[i].et->GetShortName().c_str(), nm[j]);
    // iReturn += TestOptimizerGradientComputation(mop, *vm[j], model, stepsize[i]);
    }

  // Delete both pointers
  for(i = 0; i < vt.size(); i++) 
    delete vt[i].et;
  
  for(j = 0; j < vm.size(); j++)
    delete vm[j];

  return iReturn;
}

int TestAffineTransform(const char *fnMPDE)
{
  // Load the model from file
  MedialPDE mp(fnMPDE);
  GenericMedialModel *model = mp.GetMedialModel();
  vnl_vector<double> hint = model->GetHintArray();

  // Generate a random affine transform. Because the PDE solution may change
  // under shear and non-uniform scaling, we stick just to the similarity
  // transform (A is an orthogonal matrix)
  AffineTransformDescriptor::Vec b(3, 0.0), c(3, 0.0);
  AffineTransformDescriptor::Mat A(3, 3, 0.0);

  // Generate random vectors V1, V2 for making the matrix A, plus b and c
  SMLVec3d v1, v2;
  vnl_random rnd;
  for(size_t i = 0; i < 3; i++)
    {
    b[i] = rnd.drand32(-1.0, 1.0);
    c[i] = rnd.drand32(-1.0, 1.0);
    v1[i] = rnd.drand32(-1.0, 1.0);
    v2[i] = rnd.drand32(-1.0, 1.0);
    }

  // Create the orthogonal matrix A
  A.set_column(0, v1.normalize());
  A.set_column(1, vnl_cross_3d(A.get_column(0), v2).normalize());
  A.set_column(2, vnl_cross_3d(A.get_column(0), A.get_column(1)).normalize());
  double s = rnd.drand32(0.0, 1.0); A *= s;

  // Print the transform
  cout << "Affine Matrix: " << endl << A << endl;
  cout << "A.t() * A: " << endl << A.transpose() * A << endl;
  cout << "Translation: " << b << endl;
  cout << "Center of Rotation: " << c << endl;

  // Get a list of interior points
  vector<SMLVec3d> X1;

  MedialInternalPointIterator it = model->GetInternalPointIterator(6);
  for( ; !it.IsAtEnd(); ++it)
    {
    SMLVec3d x = GetInternalPoint(it, model->GetAtomArray());
    SMLVec3d Ax = (A * (x - c)) + b + c;
    X1.push_back(Ax); 
    }

  cout << "Obtained initial point data" << endl;

  // Get the coefficients of the model
  vnl_vector<double> xCoeff = model->GetCoefficientArray();

  // Scale the hint by the square of s (hint = phi)
  hint *= s * s;

  // Apply the transform to the model
  const AffineTransformDescriptor *adt = model->GetAffineTransformDescriptor();
  model->SetCoefficientArray(adt->ApplyAffineTransform(xCoeff, A, b, c));
  model->ComputeAtoms(hint.data_block());

  // Max difference accumulator
  double maxdiff = 0.0;

  size_t q = 0;
  for(it = model->GetInternalPointIterator(6); !it.IsAtEnd(); ++it, ++q)
    {
    SMLVec3d x = GetInternalPoint(it, model->GetAtomArray());
    double diff = (x - X1[q]).two_norm();
    maxdiff = std::max(maxdiff, diff);
    }

  // Report
  cout << "Maximum difference in Affine Transform computations: " << maxdiff << endl;
  return (maxdiff > 1.0e-6);
}

int tquan(const char *name, double val, double expval)
{
  cout << name << " : " << val << " vs. " << expval << endl;
  if(fabs(val - expval) < 0.00001)
    return 0;
  else return 1;
}

/** 
 * This is a simple canned atom math test. It just checks that the methods in 
 * the medial atom computation are correct. The data supplied below is from a
 * Mathematica notebook
 */
int TestAtomMath()
{
  int rc = 0;

  MedialAtom a;
  MedialAtom da;

  a.X = SMLVec3d(134.685, 136.622, 107.994);
  a.Xu = SMLVec3d(0.999853, -0.000184419, 8.73641);
  a.Xv = SMLVec3d(0.0000424167, 1.00006, -3.10348);
  a.Xuu = SMLVec3d(0.000608536, 0.00105602, -16.3668);
  a.Xuv = SMLVec3d(-0.000217664, -0.000382739, 5.38932);
  a.Xvv = SMLVec3d(0.0000943396, 0.000155847, -2.58367);
  a.F = 0.722479; a.Fu = 0.77349; a.Fv = -0.229935; 
  a.Fuu = -1.47484; a.Fuv = 0.466966; 
  a.Fvv = -0.220057; 

  da.X = SMLVec3d(-0.347434, -0.189541, -0.111304);
  da.Xu = SMLVec3d(0.59202, -0.159739, 0.730395);
  da.Xv = SMLVec3d(-0.260392, -0.137631, -0.421155);
  da.Xuu = SMLVec3d(-0.722212, -5.02168, -2.0399);
  da.Xuv = SMLVec3d(0.253515, 2.4806, 0.721954);
  da.Xvv = SMLVec3d(-0.121972, -1.01809, -0.303003);
  da.F = 0.111147; da.Fu = -2.59863; da.Fv = 0.716595; 
  da.Fuu = 6.58777; da.Fuv = -2.70117; da.Fvv = 1.20756;

  a.ComputeDifferentialGeometry();
  a.ComputeNormalVector();
  a.ComputeBoundaryAtoms(false);
  a.ComputeBoundaryCurvature();

  MedialAtom::DerivativeTerms dt;
  a.ComputeCommonDerivativeTerms(dt);
  a.ComputeMetricTensorDerivatives(da);
  a.ComputeChristoffelDerivatives(da);
  a.ComputeBoundaryAtomDerivatives(da, dt);
  a.ComputeBoundaryCurvatureDerivative(da);

  // Check that the values match what we expect
  rc += tquan("Bnd. Mean. Curv[0]", a.xBnd[0].curv_mean, 0.574173);
  rc += tquan("Bnd. Gauss. Crv[0]", a.xBnd[0].curv_gauss, 0.00155289);
  rc += tquan("Bnd. Mean. Curv[0] Deriv", da.xBnd[0].curv_mean, 1.99231);
  rc += tquan("Bnd. Gauss. Crv[0] Deriv", da.xBnd[0].curv_gauss, 0.00449958);

  return rc;
}

int TestSmoothedImageSampler(const char *fn_image)
{
  // Read the image of interest
  FloatImage fimg;
  fimg.LoadFromFile(fn_image);

  // Get the underlying itk::Image
  typedef FloatImage::WrapperType::ImageType ImageType;
  ImageType::Pointer im = fimg.GetInternalImage()->GetInternalImage();

  // Compute the bounding box (in pixel coordinates)
  double bb_start[] = {0.0, 0.0, 0.0};
  double bb_end[3];
  for(size_t d = 0; d < 3; d++)
    bb_end[d] = im->GetBufferedRegion().GetSize()[d];

  // Create the smooth image sampler
  SmoothedImageSampler sissy(0.8, 3.5, im->GetBufferPointer(), bb_start, bb_end);

  // Sample the image at some random location
  SMLVec3d X1, X2, dX;
  for(size_t d = 0; d < 3; d++)
    {
    X1[d] = bb_start[d] + rand() * (bb_end[d] - bb_start[d]) / RAND_MAX;
    X2[d] = bb_start[d] + rand() * (bb_end[d] - bb_start[d]) / RAND_MAX;
    }

  X1[0] = 210.417;
  X1[1] = 51.7757;
  X1[2] = 110.159;
  X2[0] = 202.093;
  X2[1] = 51.8476;  
  X2[2] = 103.568;


  for(double t = 0; t < 1.0; t+=0.001)
    {
    SMLVec3d X, Xp, Xm, G, Gd, Ga;
    double f, fp[3], fm[3];
    
    // Set X
    X = X1 * (1-t) + X2 * t;

    // Do finite difference computation
    double eps = 0.0001;
    f = sissy.Sample(X.data_block(), G.data_block());
    for(size_t d = 0; d < 3; d++)
      {
      Xp = X; Xp[d] += eps;
      Xm = X; Xm[d] -= eps;
      fp[d] = sissy.Sample(Xp.data_block(), Gd.data_block());
      fm[d] = sissy.Sample(Xm.data_block(), Gd.data_block());
      Ga[d] = (fp[d] - fm[d]) / (2 * eps);
      }

    // Compute the difference
    SMLVec3d dG = G - Ga;
    double maxdiff = dG.inf_norm();

    printf("%9g %9g %9g : %9.4g / %9.4g %9.4g %9.4g / %9.4g %9.4g %9.4g / %9.4g\n",
      X[0], X[1], X[2], f, G[0], G[1], G[2], Ga[0], Ga[1], Ga[2], maxdiff);

    // cout << X << " " << f << " " << maxdiff << endl;

    if(maxdiff > eps)
      {
      cout << "X : " << X << endl;
      cout << "f : " << f << endl;
      cout << "G : " << G << endl;
      cout << "Ga : " << Ga << endl;
      cout << "fp : " << fp[0] << " " << fp[1] << " " << fp[2] << endl;
      cout << "fm : " << fm[0] << " " << fm[1] << " " << fm[2] << endl;
      cout << "MaxDiff : " << maxdiff << endl;
      // return -1;
      }
    }
  return 0;
}


int usage()
{
  cout << "testpde: MedialPDE Test Module" << endl;
  cout << "  usage: testpde TEST_ID [parameters] " << endl;
  cout << "  tests: " << endl;
  cout << "    DERIV1 XX.mpde XX.mat      Check analytic derivative PDEs." << endl;
  cout << "    DERIV2 XX.mpde             Check gradient computation in image match terms." << endl;
  cout << "    DERIV3 XX.mpde             Check variations on basis functions." << endl;
  cout << "    DERIV4 XX.mpde             Test diff. geom. operators." << endl;
  cout << "    AFFINE XX.mpde             Test affine transform computation." << endl;
  cout << "    ATOMMATH                   Test local medial geometry." << endl;
  cout << "    SAMPLE XX.img              Test image sampler code (SmoothedImageSampler)" << endl;
  cout << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Different tests that can be executed
  if(argc == 1) return usage();

  // Choose a test depending on the parameters
  if(0 == strcmp(argv[1], "DERIV1") && argc > 2)
    return TestDerivativesNoImage(argv[2]);
  else if(0 == strcmp(argv[1], "DERIV2") && argc > 2)
    {
    FloatImage *image = NULL;
    if(argc > 4)
      {
      try
        {
        image = new FloatImage();
        image->LoadFromPath(argv[3], argv[4]);
        }
      catch(...)
        {
        image = NULL;
        cerr << "Couldn't load image; proceeding without" << endl;
        }
      }
    const char *refmod = NULL;
    if(argc > 5)
      refmod = argv[5];

    return TestDerivativesWithImage(argv[2],image,refmod);
    }
  else if(0 == strcmp(argv[1], "DERIV3") && argc > 2)
    return TestBasisFunctionVariation(argv[2]);
  else if(0 == strcmp(argv[1], "DERIV4") && argc > 2)
    return TestDifferentialGeometry(argv[2]);
  else if(0 == strcmp(argv[1], "DERIV5") && argc > 2)
    return TestGradientTiming(argv[2]);
  else if(0 == strcmp(argv[1], "WEDGE"))
    return TestWedgeVolume();
  else if(0 == strcmp(argv[1], "VOLUME1") && argc > 2)
    return TestVolumeComputation(argv[2]);
  else if(0 == strcmp(argv[1], "AFFINE") && argc > 2)
    return TestAffineTransform(argv[2]);
  else if(0 == strcmp(argv[1], "ATOMMATH"))
    return TestAtomMath();
  else if(0 == strcmp(argv[1], "SAMPLE"))
    {
    if(argc < 3)
      { 
      cerr << "Usage: " << argv[0] << " SAMPLE image" << endl;
      return -1;
      }
    return TestSmoothedImageSampler(argv[2]);
    }
  else if(0 == strcmp(argv[1], "TETGEN"))
    {
    MedialPDE mp(argv[2]);
    TestTetGen(mp.GetMedialModel());
    }

  else 
    return usage();
}
