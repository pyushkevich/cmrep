#include "MedialAtomGrid.h"
#include "MedialPDESolver.h"
#include "FourierSurface.h"
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "mspline.h"
#include "glengine.h"
#include "vnl/vnl_cross.h"

// #include "imaging.txx"
#include "ImageCubeITK.h"
#include "MedialPDERenderer.h"

using namespace std;

// Define an arbitrary surface to test the differential geometry
// computation, and later to test the medial PDEs
class MyProblem1 : virtual public IMedialPDEProblem 
{
public:
  void Initialize(MedialAtom &xAtom)
    {
    float u = xAtom.u;
    float v = xAtom.v;

    // Compute the partial derivatives
    xAtom.X[0] = 0.3*u;
    xAtom.X[1] = v;
    xAtom.X[2] = 0.03*(-3 - 3*u + u*u + 4*v + 2*u*v - v*v);
    xAtom.Xu[0] = 0.3;
    xAtom.Xu[1] = 0;
    xAtom.Xu[2] = 0.03*(-3 + 2*u + 2*v);
    xAtom.Xv[0] = 0;
    xAtom.Xv[1] = 1;
    xAtom.Xv[2] = 0.03*(4 + 2*u - 2*v);
    xAtom.Xuu[0] = 0;
    xAtom.Xuu[1] = 0;
    xAtom.Xuu[2] = 0.06;
    xAtom.Xuv[0] = 0;
    xAtom.Xuv[1] = 0;
    xAtom.Xuv[2] = 0.06;
    xAtom.Xvv[0] = 0;
    xAtom.Xvv[1] = 0;
    xAtom.Xvv[2] = -0.06;

    // Set the riemannian laplacian
    xAtom.xLapR = -1.0;

    // Compute the normal vector
    xAtom.ComputeNormalVector();
    }
};

// Define an arbitrary surface to test the differential geometry
// computation, and later to test the medial PDEs
class MyProblem2 : virtual public IMedialPDEProblem {
public:
  void ComputeJet2(MedialPoint *mp)
    {
    float u = mp->u;
    float v = mp->v;
    mp->F[0] = u; mp->F[1] = v; mp->F[2] = 0.0;
    mp->Fu[0] = 1; mp->Fu[1] = 0; mp->Fu[2] = 0;
    mp->Fv[0] = 0; mp->Fv[1] = 1; mp->Fv[2] = 0;
    mp->Fuu[0] = 0; mp->Fuu[1] = 0; mp->Fuu[2] = 0;
    mp->Fuv[0] = 0; mp->Fuv[1] = 0; mp->Fuv[2] = 0;
    mp->Fvv[0] = 0; mp->Fvv[1] = 0; mp->Fvv[2] = 0;
    mp->F3[0] = mp->F3raw[0] = 0.0;
    mp->F3[1] = mp->F3raw[1] = 0.0;
    mp->F3[2] = mp->F3raw[2] = 1.0;
    }
  double ComputeLaplacian(double u, double v)
    { return -1.0; }
};

class PDESplineWrapper : virtual public IMedialPDEProblem
{
public:
  /** A wrapper around a medial spline */
  PDESplineWrapper(DynamicBSpline2D *spline, unsigned int level)
  { 
    this->spline = spline;
    this->level = level;
    sgd = SplineGridDefinition::newSGDForCache(spline, level); 
  }
    
  /* Function to get the 2-jet at u, v pair */
  void ComputeJet2(MedialPoint *mp)
  {
    // Find the u, v index for the input coordinates
    unsigned int iu = sgd->getIndexForUV(0, mp->u);
    unsigned int iv = sgd->getIndexForUV(1, mp->v);

    cout << mp->u << " : " << iu << " ;" << mp->v << " : " << iv << endl;

    // Interpolate the grid point for this medial point
    spline->interpolateGridPoint(*sgd, iu, iv, 0, 0, 0, 2, mp->F.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 1, 0, 0, 2, mp->Fu.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 0, 1, 0, 2, mp->Fv.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 2, 0, 0, 2, mp->Fuu.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 1, 1, 0, 2, mp->Fuv.data_block());
    spline->interpolateGridPoint(*sgd, iu, iv, 0, 2, 0, 2, mp->Fvv.data_block());

    // Compute the normal
    SMLVec3f Nraw = cross_3d(mp->Xu(), mp->Xv());
    SMLVec3f N = Nraw / Nraw.magnitude();
    memcpy(mp->F3raw.data_block(), Nraw.data_block(), 3 * sizeof(float));
    memcpy(mp->F3.data_block(), N.data_block(), 3 * sizeof(float));
  }
  
  double ComputeLaplacian(double u, double v)
    { return -0.7; }  
    
private:
  // Reference to the medial spline
  DynamicBSpline2D *spline;
  SplineGridDefinition *sgd;
  unsigned int level;
};

class PDEFourierWrapper : virtual public IMedialPDEProblem
{
public:
  PDEFourierWrapper(FourierSurface *xSurface)
    {
    this->xSurface = xSurface;
    }

  void ComputeJet2(MedialPoint *mp)
    {
    float u = mp->u, v = mp->v;

    for(unsigned int d = 0; d < 3; d++)
      {
      mp->F[d]   = xSurface->Evaluate(u, v, 0, 0, d);
      mp->Fu[d]  = xSurface->Evaluate(u, v, 1, 0, d);
      mp->Fv[d]  = xSurface->Evaluate(u, v, 0, 1, d);
      mp->Fuu[d] = xSurface->Evaluate(u, v, 2, 0, d);
      mp->Fuv[d] = xSurface->Evaluate(u, v, 1, 1, d);
      mp->Fvv[d] = xSurface->Evaluate(u, v, 0, 2, d);
      }

    // Compute the normal
    SMLVec3f Nraw = cross_3d(mp->Xu(), mp->Xv());
    SMLVec3f N = Nraw / Nraw.magnitude();
    memcpy(mp->F3raw.data_block(), Nraw.data_block(), 3 * sizeof(float));
    memcpy(mp->F3.data_block(), N.data_block(), 3 * sizeof(float));
    }

  double ComputeLaplacian(double u, double v)
    {
    double ramp = (0.5 - abs(0.5 - u)) * (0.5 - abs(0.5 - v));
    double test = rand() * 0.5 / RAND_MAX;
    return -1.0; 
    }
  
private:
  FourierSurface *xSurface;
};

void testFourierFit(FourierSurface *xFourier)
{
  // Decide how many points to interpolate
  unsigned int nSide = 11, nPoints = nSide * nSide;
  double uStep = 1.0 / (nSide - 1);

  // Allocate arrays of points and coordinates
  float xPoints[nPoints], yPoints[nPoints], zPoints[nPoints];
  float uPoints[nPoints], vPoints[nPoints];

  // Create an array of points
  unsigned int i = 0;
  for(unsigned int u = 0; u < nSide; u++)
    for(unsigned int v = 0; v < nSide; v++)
      {
      float uu = u * uStep, vv = v * uStep;

      uPoints[i] = uu;
      vPoints[i] = vv;
      xPoints[i] = 0.5 * uu + 0.25;
      yPoints[i] = vv;
      zPoints[i] = ((uu - 0.5) * (uu - 0.5) + (vv - 0.5) * (vv - 0.5)) * 0.25;
      ++i;
      }

  // Peform the fit
  xFourier->FitData(0, nPoints, uPoints, 1, vPoints, 1, xPoints, 1, 4, 4);
  xFourier->FitData(1, nPoints, uPoints, 1, vPoints, 1, yPoints, 1, 4, 4);
  xFourier->FitData(2, nPoints, uPoints, 1, vPoints, 1, zPoints, 1, 4, 4);
}

struct DiscreteAtom 
{ 
  double u, v, x, y, z; 
  unsigned int iu, iv; 
};

void FitDiscreteMRep(const char *file, FourierSurface *surface)
{
  ifstream fin(file, ios_base::in);
  int iu, iv, iend;
  double x, y, z, d;

  // Vector to store atoms
  typedef vector<DiscreteAtom> AList;
  AList xAtoms;

  // Get the max u and v values
  unsigned int uMax = 0, vMax = 0;

  // Read atoms
  bool done = false;
  while(!done)
    {
    DiscreteAtom atom;
    
    // Read the atom
    fin >> atom.iu; fin >> atom.iv; fin >> iend; 
    fin >> atom.x; fin >> atom.y; fin >> atom.z; 
    fin >> d; fin >> d; fin >> d; fin >> d; 
    fin >> d; fin >> d; fin >> d;

    if(uMax < atom.iu) uMax = atom.iu;
    if(vMax < atom.iv) vMax = atom.iv;

    if(!fin.good())
      { done = true; }
    else
      {
      xAtoms.push_back(atom);
      }
    }

  // Scale by u and v to unit square
  for(AList::iterator it = xAtoms.begin(); it!=xAtoms.end(); ++it)
    { it->u = it->iu * 1.0 / uMax; it->v = it->iv * 1.0 / vMax; }

  // Create double arrays
  float *xx = new float[xAtoms.size()];
  float *yy = new float[xAtoms.size()];
  float *zz = new float[xAtoms.size()];
  float *uu = new float[xAtoms.size()];
  float *vv = new float[xAtoms.size()];

  for(unsigned int i = 0; i < xAtoms.size(); i++)
    {
    xx[i] = xAtoms[i].x;
    yy[i] = xAtoms[i].y;
    zz[i] = xAtoms[i].z;
    uu[i] = xAtoms[i].u;
    vv[i] = xAtoms[i].v;
    }

  // Perform the fitting on x, y and z
  unsigned int stride = sizeof(DiscreteAtom);
  surface->FitData(0, xAtoms.size(), uu, 1, vv, 1, xx, 1, 5, 5);
  surface->FitData(1, xAtoms.size(), uu, 1, vv, 1, yy, 1, 5, 5);
  surface->FitData(2, xAtoms.size(), uu, 1, vv, 1, zz, 1, 5, 5);

  // Clean up
  delete xx; delete yy; delete zz; delete uu; delete vv;
}

/** Write a match class for matching our solutions to images */
class MRepImageMatcher : virtual public BoundaryMeasure
{
public:
  // Compute the match between an m-rep and image
  virtual float ComputeMatch(IMedialSurfacePatch *xPatch)
    {
    float xMatch, xArea;
    xMatch = xPatch->IntegrateBoundaryMeasure(this, xArea);
    return xMatch / xArea;
    }
};

class MRepEdgeMatcher : public MRepImageMatcher
{
public:
  // Initialize with an image
  void SetEdgeImage(AbstractImage3D *xImage)
    { this->xImage = xImage; }

  // Compute the match
  float computeBoundaryMeasure(const MedialPoint &mp, int d)
    {
    const float *x = mp.bp[d].X.data_block();
    return xImage->interpolateVoxel(x[0], x[1], x[2]);
    }
  
  float computeCrestBoundaryMeasure(const MedialPoint &mp)
    { return computeBoundaryMeasure(mp, 0); }
  
private:
  AbstractImage3D *xImage;
};

#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"

ImageCubeITK<float> *TestImageIO(const char *fname)
{
  cout << "Reading image " << fname << endl;

  // Load an image from ITK
  typedef itk::Image<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(fname);
  fltReader->Update();
  ImageType::Pointer imgInput = fltReader->GetOutput();

  cout << "Creating image cube" << fname << endl;

  // Create an image cube from this image
  ImageCubeITK<float> *cubeInput = new ImageCubeITK<float>;
  cubeInput->SetImage(imgInput, 1.0);

  // Get the size object
  ImageType::SizeType szImage = imgInput->GetBufferedRegion().GetSize();

  cout << "Generating random indices" << fname << endl;

  // Generate indices around a sphere of radius 0.6 in the image (some
  // samples fall outside, some inside)
  unsigned int nSamples = 300;
  unsigned int nTests = nSamples * nSamples * nSamples, iTest, q;
  vector<SMLVec3f> xIndex;
  for(unsigned int iTheta = 0; iTheta < nSamples; ++iTheta)
    for(unsigned int iPhi = 0; iPhi < nSamples; ++iPhi)
      for(unsigned int iRad = 0; iRad < nSamples; ++iRad)
        {
        double theta = iTheta * M_PI / nSamples; 
        double phi = iPhi * 2.0 * M_PI / nSamples; 
        double r = iRad * 0.6 / nSamples;

        double x = r * cos(phi) * sin(theta);
        double y = r * sin(phi) * sin(theta);
        double z = r * cos(theta);

        xIndex.push_back(SMLVec3f(x,y,z));
        }
  /*
  //  Generate an array of random indices into the image
  unsigned int nTests = 2000000, iTest, q;
  vector<SMLVec3f> xIndex(nTests);
  for(iTest = 0; iTest < nTests; iTest++)
    {
    xIndex[iTest][0] = rand() * 1.0 / RAND_MAX * szImage[0];
    xIndex[iTest][1] = rand() * 1.0 / RAND_MAX * szImage[1];
    xIndex[iTest][2] = rand() * 1.0 / RAND_MAX * szImage[2];
    }
  */

  // Interpolate the image at these locations using the cube
  double tCube, tITK;
  float xTestCube = 0, xTestITK = 0;

  cout << "Running ITK interpolation" << fname << endl;
  
  // Interpolate using the ITK code
  typedef itk::LinearInterpolateImageFunction<ImageType, float> FuncType;
  FuncType::Pointer funInterp = FuncType::New();
  funInterp->SetInputImage(imgInput);

  cout << imgInput->GetOrigin() << endl;
  cout << imgInput->GetSpacing() << endl;
  cout << imgInput->GetBufferedRegion() << endl;

  tITK = clock();
  for(q = 0; q < 10; q++)
    for(iTest = 0; iTest < nTests; iTest++)
      {
      itk::ContinuousIndex<float, 3> idx(xIndex[iTest].data_block());
      if(funInterp->IsInsideBuffer(idx))
        {
      	float xVal = funInterp->EvaluateAtContinuousIndex(idx);
        xTestITK += xVal;
        }
      }
  tITK = (clock() - tITK) * 1.0 / CLOCKS_PER_SEC;

  cout << "Running image cube interpolation" << fname << endl;
  tCube = clock();
  for(q = 0; q < 10; q++)
    for(iTest = 0; iTest < nTests; iTest++)
      {
      float *p = xIndex[iTest].data_block();
      float xVal = cubeInput->interpolateVoxel(p[0], p[1], p[2]);
      xTestCube += xVal;
      }
  tCube = (clock() - tCube) * 1.0 / CLOCKS_PER_SEC;
  
  cout << "ITK returned " << xTestITK << " in " << (1000 * tITK) << " ms." << endl;
  cout << "CUB returned " << xTestCube << " in " << (1000 * tCube) << " ms." << endl;

  return cubeInput;
}

int main(int argc, char *argv[])
{
  TestCartesianGrid(); // Remove! 
  
  // Command line parameter variables
  unsigned int xResolution = 8;
  char *fnImage = NULL, *fnModel = NULL;

  // Read command line parameters 
  for(unsigned int iArg = 1; iArg < argc; iArg++)
    {
    if(strcmp(argv[iArg],"-r") == 0 && iArg < argc - 1)
      xResolution = atoi(argv[++iArg]);
    if(strcmp(argv[iArg],"-i") == 0 && iArg < argc - 1)
      fnImage = argv[++iArg];
    if(strcmp(argv[iArg],"-m") == 0 && iArg < argc - 1)
      fnModel = argv[++iArg];
    }
  
  if(fnImage)
    {
    ImageCubeITK<float> *cube = TestImageIO(fnImage);return 0;
    }

  /* Step 1. Test the Differential Geometry code */
  double u = 0.3, v = 0.4;
  
  // Build and display a spline
  unsigned int xPointsU = 6, xPointsV = 6;
  DynamicBSpline2D *spline = new DynamicBSpline2D(xPointsU, xPointsV, 3, 0);
  for(unsigned int i = 0; i <= xPointsU; i++)
    for(unsigned int j = 0; j <= xPointsV; j++)
      spline->setControl(i, j, 2, 0.0);
    
  PDESplineWrapper psw(spline, 3);
    
  // Compute the Jet
  MyProblem1 p1;
  MedialPoint mp; mp.u = u; mp.v = v;
  p1.ComputeJet2(&mp);

  // Another problem and jet
  MyProblem2 p2;
  p2.ComputeJet2(&mp);

  // Compute the differential geometry
  GeometryDescriptor<float> gd;
  gd.SetJet(&mp);
  gd.PrintSelf(cout);

  // Try a fourier problem
  FourierSurface xFourier(5,5,3);
  srand(clock());
  for(unsigned int k = 0; k < 4; k++)
    for(unsigned int l = 0; l < 4; l++)
      {
      double v = (rand() * 1.0 / RAND_MAX) - 0.5;
      v *= pow(0.5, (double)(k+l));
      // v *= exp((double)(-k - l)) / exp(2.0);
      xFourier.SetCoefficient(k, l, 2, v);
      cout << "u: " << k << ";  v: " << l << "  v: " << v << endl;
      }

  // Test data fitting
  if(fnModel)
    FitDiscreteMRep(fnModel, &xFourier);
  else
    testFourierFit(&xFourier);

  // Create a wrapper around the Fourier problem
  PDEFourierWrapper pfw(&xFourier);
  
  /* Step 2. Evaluate the equation at each site */
  // Initialize the solver
  unsigned int m = xResolution * xPointsU;
  unsigned int n = xResolution * xPointsV;
  MedialPDESolver mps(m+1, n+1);

  time_t t = clock();
  mps.Solve(&pfw);
  cout << "*** Elapsed Time: " << ((clock() - t) * 1.0 / CLOCKS_PER_SEC) << " ms." << " ***" << endl;
  
  // Test the image matching code
  if(fnImage)
    {
    ImageCubeITK<float> *cube = TestImageIO(fnImage);
    MRepEdgeMatcher matcher;
    matcher.SetEdgeImage(cube);
    cout << "Image Match: " << matcher.ComputeMatch(&mps) << endl;
    return -1;
    }

  // Create a renderer
  PDESplineRenderer *rndPDESpline = new PDESplineRenderer(&mps);

  // Initialize the GL environment
  GLDisplayDriver::init(argc, argv);

  // Add the important renderers
  GLDisplayDriver::addRenderer(new DefaultLightRenderer());
  GLDisplayDriver::addRenderer(new StarfieldRenderer());

  // Add our spline renderer
  GLDisplayDriver::addRenderer(rndPDESpline);

  // Start GLUT
  glutMainLoop();
}




