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

void ExportMedialMeshToVTK(
    GenericMedialModel *xModel, ITKImageWrapper<float> *xImage, const char *file);

void ExportBoundaryMeshToVTK(
    GenericMedialModel *xModel, ITKImageWrapper<float> *xImage, const char *file);

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

using namespace Ipopt;
using namespace gnlp;


typedef itk::Image<float, 3> ImageType;






std::vector<SMLVec3d> find_closest(GenericMedialModel *model, vtkPolyData *target)
{
  // Create locator for finding closest points
  vtkCellLocator *loc = vtkCellLocator::New();
  loc->SetDataSet(target);
  loc->CacheCellBoundsOn();
  loc->BuildLocator();

  // Output vector
  std::vector<SMLVec3d> cp(model->GetNumberOfBoundaryPoints());

  // Compute thickness values
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    SMLVec3d x = GetBoundaryPoint(it, model->GetAtomArray()).X;

    double xs[3], d2, d;
    int subid;
    vtkIdType cellid;

    loc->FindClosestPoint(x.data_block(), xs, cellid, subid, d2);

    cp[it.GetIndex()] = SMLVec3d(xs);
    }

  return cp;
}

#include <vtkPointLocator.h>

typedef std::pair<int, SMLVec3d> PointMatch;

std::vector<PointMatch>
find_closest_to_mesh(vtkPolyData *target, GenericMedialModel *model, int nSamples)
{
  // Create a VTK points object
  vtkPoints *out_pts = vtkPoints::New();
  out_pts->Allocate(model->GetNumberOfBoundaryPoints());

  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = it.GetIndex();
    out_pts->InsertNextPoint(GetBoundaryPoint(it,model->GetAtomArray()).X.data_block());
    }

  vtkPolyData *poly = vtkPolyData::New();
  poly->SetPoints(out_pts);

  // Create locator for finding closest points
  vtkPointLocator *loc = vtkPointLocator::New();
  loc->SetDataSet(poly);
  loc->BuildLocator();

  // Sample points from the target mesh
  std::vector<PointMatch> result;
  for(int i = 0; i < nSamples; i++)
    {
    int q = rand() % target->GetNumberOfPoints();
    SMLVec3d xTarget(target->GetPoint(q));

    vtkIdType id = loc->FindClosestPoint(xTarget.data_block());

    result.push_back(std::make_pair((int) id, xTarget));
    }

  return result;
}

// Callback functions for TOMS
IPOptProblemInterface *globopt;
void toms_calcf(int &n, double *x, int &nf, double &f, int *dummy1, double *dummy2, void *info)
{
  globopt->eval_f(n, x, true, f);
}

void toms_calcg(int &n, double *x, int &nf, double *g, int *, double *, void *info)
{
  globopt->eval_grad_f(n, x, true, g);
}


void run_toms(IPOptProblemInterface *ip, ConstrainedNonLinearProblem *p)
{
  // SOLVE USING TOMS
  globopt = ip;
  int nCoeff = p->GetNumberOfVariables();
  double *scaling = (double *) malloc(nCoeff * sizeof(double));
  std::fill_n(scaling, nCoeff, 1.0);

  double *x = (double *) malloc(nCoeff * sizeof(double));
  for(int i = 0; i < nCoeff; i++)
    x[i] = p->GetVariableValue(i);

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
  iv[mxiter_ - 1] = 100;
  iv[mxfcal_ - 1] = 10 * 100;
  // iv[19 - 1] = 0;
  // iv[22 - 1] = 0;
  // iv[24 - 1] = 0;

  sumsl_(
        nCoeff, scaling, x,
        &toms_calcf, &toms_calcg,
        iv, liv, lv, v,
        NULL, NULL, NULL);

}

Expression *TetraHedronVolume(
    Problem *p, GenericMedialModel *model,
    std::vector<Expression *> a,
    std::vector<Expression *> b,
    std::vector<Expression *> c,
    std::vector<Expression *> d)
{
  // Compute the three edge vectors
  std::vector<Expression *> Q(3, NULL), R(3, NULL), S(3, NULL);
  for(int j = 0; j < 3; j++)
    {
    Q[j] = new BinaryDifference(p, b[j], a[j]);
    R[j] = new BinaryDifference(p, c[j], a[j]);
    S[j] = new BinaryDifference(p, d[j], a[j]);
    }

  // Compute the vector triple product
  Expression *vol =
      new TernarySum(
        p,
        new BinaryProduct(
          p, Q[0], new BinaryDifference(
            p,
            new BinaryProduct(p, R[1], S[2]),
            new BinaryProduct(p, R[2], S[1]))),
        new BinaryProduct(
          p, Q[1], new BinaryDifference(
            p,
            new BinaryProduct(p, R[2], S[0]),
            new BinaryProduct(p, R[0], S[2]))),
        new BinaryProduct(
          p, Q[2], new BinaryDifference(
            p,
            new BinaryProduct(p, R[0], S[1]),
            new BinaryProduct(p, R[1], S[0]))));

  return vol;
}


double GetCentralDifference(Expression *ex, Variable *v, double delta=1e-5)
{
  // Do the perturbation
  double val = v->Evaluate();
  v->SetValue(val + delta);
  double f2 = ex->Evaluate();
  v->SetValue(val - delta);
  double f1 = ex->Evaluate();
  v->SetValue(val);

  return (f2 - f1) / (2 * delta);
}

void TestExpressionRandomDerivative(
    Problem *p,
    Expression *ex,
    const char *nickname,
    int order)
{
  // Check order
  if(order == 0)
    return;

  // Get the list of dependent variables
  const Problem::Dependency &depvar = p->GetDependentVariables(ex);
  if(depvar.size() == 0)
    return;

  // Pick a random dependent variable to test
  int q = rand() % depvar.size();
  Problem::Dependency::const_iterator it = depvar.begin();
  for(int p = 0; p < q; p++) ++it;
  Variable *v = *it;

  // Test with respect to this variable
  Expression *pd = p->GetPartialDerivative(ex, v);
  if(!pd)
    return;

  double dAnalytic = pd ? pd->Evaluate() : 0.0;
  double dCentralDiff = GetCentralDifference(ex, v);

  printf("D[%10s,%10s,%d]: %12.8f  %12.8f  %12.8f\n",
         nickname, v->GetName().c_str(), order,
         dAnalytic, dCentralDiff, std::fabs(dAnalytic-dCentralDiff));

  // Test higher order derivatives
  TestExpressionRandomDerivative(p, pd, nickname, order-1);
}

void DerivativeTest(ConstrainedNonLinearProblem *p, int nTests)
{
  printf("TEST [%12s]: %12s  %12s  %12s\n",
         "Variable", "Analytic", "CentralDiff", "Delta");

  // The derivatives of the objective function
  for(int i = 0; i < nTests; i++)
    {
    TestExpressionRandomDerivative(p, p->GetObjective(), "obj", 2);
    }

  // Test th derivatives of the constraints
  for(int i = 0; i < nTests; i++)
    {
    // Pick a random constraint
    int iCon = rand() % p->GetNumberOfConstraints();
    Expression *con = p->GetConstraint(iCon);

    char buffer[16];
    sprintf(buffer, "Con_%d", iCon);

    TestExpressionRandomDerivative(p, con, buffer, 2);
    }
}




#include "vtkFloatArray.h"
#include "vtkPointData.h"

void SaveSamples(std::vector<std::vector<Expression *> > sampleX,
                  std::vector<Expression *> sampleF,
                  const char *filename)
{
  vtkPoints *pts = vtkPoints::New();
  pts->Allocate(sampleX.size());

  vtkFloatArray *arr = vtkFloatArray::New();
  arr->SetNumberOfComponents(1);
  arr->Allocate(sampleX.size());

  for(int i = 0; i < sampleX.size(); i++)
    {
    pts->InsertNextPoint(sampleX[i][0]->Evaluate(),
                         sampleX[i][1]->Evaluate(),
                         sampleX[i][2]->Evaluate());
    arr->InsertNextTuple1(sampleF[i]->Evaluate());
    }

  vtkPolyData *poly = vtkPolyData::New();
  poly->SetPoints(pts);
  poly->GetPointData()->SetScalars(arr);

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(poly);
  writer->SetFileName(filename);
  writer->Update();

}

void SaveGradient(
    ConstrainedNonLinearProblem *p,
    SubdivisionMedialModel *model,
    std::vector<std::vector<Expression *> > X,
    Expression *f,
    const char *filename)
{
  vtkPoints *pts = vtkPoints::New();
  pts->Allocate(model->GetNumberOfBoundaryPoints());

  vtkFloatArray *arr = vtkFloatArray::New();
  arr->SetNumberOfComponents(3);
  arr->Allocate(model->GetNumberOfBoundaryPoints());
  arr->SetName("Gradient");

  for(int i = 0; i < X.size(); i++)
    {
    pts->InsertNextPoint(X[i][0]->Evaluate(),
                         X[i][1]->Evaluate(),
                         X[i][2]->Evaluate());
    Expression *dx = p->GetPartialDerivative(f, (Variable *)X[i][0]);
    Expression *dy = p->GetPartialDerivative(f, (Variable *)X[i][1]);
    Expression *dz = p->GetPartialDerivative(f, (Variable *)X[i][2]);
    arr->InsertNextTuple3(dx->Evaluate(), dy->Evaluate(), dz->Evaluate());
    }

  vtkPolyData *poly = vtkPolyData::New();
  poly->SetPoints(pts);
  poly->GetPointData()->SetScalars(arr);

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(poly);
  writer->SetFileName(filename);
  writer->Update();
}

#include "itkOrientedRASImage.h"
#include "itk_to_nifti_xform.h"
#include "vtkCubeSource.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkCell.h"
#include "vtkSmartPointer.h"
#include "vtkTriangleFilter.h"
#include "vtkLinearSubdivisionFilter.h"

/**
  Generate a tetrahedral mesh for the model exerior using TetGen.
  */
void CreateTetgenMesh(GenericMedialModel *model,
                      FloatImage *image,
                      VarVecArray &X,
                      ConstrainedNonLinearProblem *problem)
{
  // Create a mesh object to populate
  tetgenio in;
  in.initialize();

  // Create the cube
  vtkSmartPointer<vtkCubeSource> cube = vtkCubeSource::New();
  cube->SetBounds(-0.5, 0.5 + image->GetInternalImage()->GetImageSize(0),
                  -0.5, 0.5 + image->GetInternalImage()->GetImageSize(1),
                  -0.5, 0.5 + image->GetInternalImage()->GetImageSize(2));

  // Subdivide the cube
  vtkSmartPointer<vtkTriangleFilter> fltTri = vtkTriangleFilter::New();
  fltTri->SetInputConnection(cube->GetOutputPort());

  vtkSmartPointer<vtkLinearSubdivisionFilter> fltSub = vtkLinearSubdivisionFilter::New();
  fltSub->SetInputConnection(fltTri->GetOutputPort());
  fltSub->SetNumberOfSubdivisions(2);

  // Get the transform matrix
  vnl_matrix_fixed<double, 4, 4> TS = ConstructNiftiSform(
        image->GetInternalImage()->GetInternalImage()->GetDirection().GetVnlMatrix(),
        image->GetInternalImage()->GetInternalImage()->GetOrigin().GetVnlVector(),
        image->GetInternalImage()->GetInternalImage()->GetSpacing().GetVnlVector());

  vtkSmartPointer<vtkCleanPolyData> clean = vtkCleanPolyData::New();
  clean->SetInputConnection(fltSub->GetOutputPort());

  vtkSmartPointer<vtkTransform> tran = vtkTransform::New();
  tran->SetMatrix(TS.data_block());

  // Create the transform filter
  vtkSmartPointer<vtkTransformFilter> tf = vtkTransformFilter::New();
  tf->SetInputConnection(clean  ->GetOutputPort());
  tf->SetTransform(tran);
  tf->Update();

  // Get the transformed cube
  vtkSmartPointer<vtkPolyData> tcube = dynamic_cast<vtkPolyData *>(tf->GetOutput());

  // Initialize all the points (number of points, plus six cube vertices)
  in.numberofpoints = model->GetNumberOfBoundaryPoints() + tcube->GetNumberOfPoints();
  in.pointlist = new REAL[in.numberofpoints * 3];

  in.pointmarkerlist = new int[in.numberofpoints];

  std::fill(in.pointmarkerlist,
            in.pointmarkerlist + model->GetNumberOfBoundaryPoints(), 1);
  std::fill(in.pointmarkerlist + model->GetNumberOfBoundaryPoints(),
            in.pointmarkerlist + in.numberofpoints, 2);

  /*
  in.numberofpointattributes = 1;
  in.pointattributelist = new REAL[in.numberofpoints];
  std::fill(in.pointattributelist,
            in.pointattributelist + model->GetNumberOfBoundaryPoints(), 1.0);
  std::fill(in.pointattributelist + model->GetNumberOfBoundaryPoints(),
            in.pointattributelist + in.numberofpoints, 2.0);
            */

  // Fill out the point array
  for(MedialBoundaryPointIterator bit = model->GetBoundaryPointIterator();
      !bit.IsAtEnd(); ++bit)
    {
    REAL *p = in.pointlist + 3 * bit.GetIndex();
    SMLVec3d X = GetBoundaryPoint(bit, model->GetAtomArray()).X;
    p[0] = X[0]; p[1] = X[1]; p[2] = X[2];
    }

  // Pass in the cube vertices
  REAL *p = in.pointlist + model->GetNumberOfBoundaryPoints() * 3;
  for(int i = 0; i < tcube->GetNumberOfPoints(); i++)
    {
    for(int j = 0; j < 3; j++)
      *(p++) = tcube->GetPoint(i)[j];
    }

  // Create the faces in the mesh
  in.numberoffacets = model->GetNumberOfBoundaryTriangles() + tcube->GetNumberOfCells();
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

  const int TRIMARK = 100000;
  const int CUBMARK = 200000;

  // Initialize all of the facets
  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Set up the facet
    tetgenio::facet &f = in.facetlist[trit.GetIndex()];
    f.numberofpolygons = 1;
    f.numberofholes = 0;
    f.holelist = NULL;
    f.polygonlist = new tetgenio::polygon[1];
    f.polygonlist[0].numberofvertices = 3;
    f.polygonlist[0].vertexlist = new int[3];
    f.polygonlist[0].vertexlist[0] = trit.GetBoundaryIndex(0);
    f.polygonlist[0].vertexlist[1] = trit.GetBoundaryIndex(1);
    f.polygonlist[0].vertexlist[2] = trit.GetBoundaryIndex(2);

    // Pass the facet into list
    in.facetmarkerlist[trit.GetIndex()] = TRIMARK + trit.GetIndex();
    }

  // Initialize the cube's facets
  tetgenio::facet *fp = in.facetlist + model->GetNumberOfBoundaryTriangles();
  int *fmp = in.facetmarkerlist + model->GetNumberOfBoundaryTriangles();
  for(int i = 0; i < tcube->GetNumberOfCells(); i++)
    {
    vtkCell *cell = tcube->GetCell(i);
    fp->numberofpolygons = 1;
    fp->numberofholes = 0;
    fp->holelist = NULL;
    fp->polygonlist = new tetgenio::polygon[1];
    fp->polygonlist[0].numberofvertices = cell->GetNumberOfPoints();
    fp->polygonlist[0].vertexlist = new int[cell->GetNumberOfPoints()];
    for(int j = 0; j < cell->GetNumberOfPoints(); j++)
      {
      fp->polygonlist[0].vertexlist[j] =
          cell->GetPointId(j) + model->GetNumberOfBoundaryPoints();
      }
    *fmp++ = CUBMARK + i;
    ++fp;
    }

  // Add a hole for the volume
  /*
  in.numberofholes = 1;
  in.holelist = new double[3];
  in.holelist[0] = model->GetAtomArray()[0].X[0];
  in.holelist[1] = model->GetAtomArray()[0].X[1];
  in.holelist[2] = model->GetAtomArray()[0].X[2];
  */

  // Save the mesh
  in.save_nodes(const_cast<char *>("mytest"));
  in.save_poly(const_cast<char *>("mytest"));

  // Create an output mesh object
  tetgenio out;
  out.initialize();

  // Create the options
  tetgenbehavior tb;
  tb.parse_commandline(const_cast<char *>("-p -q3.2 -YY"));

  // Perform tetrahedralization
  tetrahedralize(&tb, &in, &out);

  // Now let's see what we got out
  printf("TETGEN result: %d tets, %d points\n",
         out.numberoftetrahedra,
         out.numberofpoints);

  // We need to convert the tetrahedralization into a set of variables and
  // constraints. There are three types of nodes, with attributes as follows:
  //   0  -  nodes inserted by tetgen on the interior.
  //   1  -  nodes in the model
  //   2  -  nodes in the cube

  // We want to create new variables for the 0-nodes. We create constants for
  // the coordinates of the 1-nodes.
  VarVecArray Y(out.numberofpoints, VarVec(3, NULL));
  for(int i = 0; i < out.numberofpoints; i++)
    {
    int type = out.pointmarkerlist[i];
    for(int j = 0; j < 3; j++)
      {
      if(type == 1)
        {
        // Medial mesh point - copy from X
        Y[i][j] = X[i][j];
        }
      else if(type == 2)
        {
        // Cube mesh point - create a new constant
        Y[i][j] = new Constant(problem, out.pointlist[i*3+j]);
        }
      else if(type == 0)
        {
        // Point created by tetgen
        Y[i][j] = problem->AddVariable("Yij", out.pointlist[i*3+j]);
        }
      }
    }

  // Finally, we create volume constraints for all the tetrahedra
  for(int i = 0; i < out.numberoftetrahedra; i++)
    {
    int *tet = out.tetrahedronlist + i * 4;

    Expression *vol =
        TetraHedronVolume(problem, model, Y[tet[0]], Y[tet[1]], Y[tet[2]], Y[tet[3]]);

    double tv = vol->Evaluate();

    problem->AddConstraint(vol, 0.1 * tv, 100 * tv);
    }
}



int main(int argc, char *argv[])
{
  // The first parameter is the cm-rep to start from
  char *cmrepfile = argv[1];

  // The second is the VTK mesh for ICP
  char *targetmesh = argv[2];

  // The image to fit to
  char *targetimage = argv[3];

  // Load and process the mrep
  MedialPDE mrep(cmrepfile);
  SubdivisionMedialModel *model =
      dynamic_cast<SubdivisionMedialModel *>(mrep.GetMedialModel());

  // Load the target mesh
  vtkPolyData *target = ReadVTKMesh(targetmesh);

  // Load the target image
  BinaryImage imgBinary;
  imgBinary.LoadFromFile(targetimage);

  // Smooth and turn into floating point
  FloatImage *image = new FloatImage();
  image->SetToBlurredBinary(&imgBinary, 0.2);
  image->SetOutsideValue(-1.0);

  // The new jet interpolator
  ImageJetInterpolator interp_new;
  interp_new.SetInputImage(image->GetInternalImage()->GetInternalImage(), -1.0);

  // Test the interpolator
  SMLVec3d xProbe(-69.6316, -118.364, 56.4);

  ImageJetInterpolator::Vec3 gNew;
  ImageJetInterpolator::Mat3 hNew;
  double newval;
  interp_new.EvaluateAtPhysicalPoint(xProbe, newval, gNew, hNew);
  printf("New value %f, gradient [%f %f %f]\n", newval, gNew[0], gNew[1], gNew[2]);

  for(int i = 0; i < 3; i++)
    {
    SMLVec3d eps(0, 0, 0); eps[i] = 0.001;
    ImageJetInterpolator::Vec3 g1, g2;
    ImageJetInterpolator::Mat3 hDummy;
    double f1, f2;
    interp_new.EvaluateAtPhysicalPoint(xProbe+eps, f2, g2, hDummy);
    interp_new.EvaluateAtPhysicalPoint(xProbe-eps, f1, g1, hDummy);
    printf("Dx[%d]:   A = %f   N = %f\n", i, gNew[i], (f2-f1) / 0.002);

    for(int j = 0; j < 3; j++)
      printf("DDx[%d][%d]:   A = %f   N = %f\n", i, j, hNew[i][j], (g2[j]-g1[j]) / 0.002);
    }

  // Create the optimization problem
  ConstrainedNonLinearProblem *p = new ConstrainedNonLinearProblem();

  // For each point on the model, find the closest target point
  std::vector<SMLVec3d> targetPoint = find_closest(model, target);

  // The boundary positions
  VarVecArray X(model->GetNumberOfBoundaryPoints(), VarVec(3, NULL));

  // The medial positions
  VarVecArray M(model->GetNumberOfAtoms(), VarVec(3, NULL));

  // The radius values
  VarVec R(model->GetNumberOfAtoms(), NULL);

  // The boundary normal vectors
  VarVecArray N(model->GetNumberOfBoundaryPoints(), VarVec(3, NULL));

  // The vectors U from the medial axis to the boundary. They are the same
  // as -R*N. They appear in multiple places though, so we should store them
  VarVecArray U(model->GetNumberOfBoundaryPoints(), VarVec(3, NULL));

  // The triangle areas on the boundary surface (to avoid square roots)
  VarVec taX(model->GetNumberOfBoundaryTriangles(), NULL);

  // The triangle areas on the medial surface
  VarVec taM(model->GetNumberOfTriangles(), NULL);

  // A buffer for making variable names
  char buffer[64];


  // Configure the medial point variables
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = (int) it.GetIndex();
    SMLVec3d x = GetBoundaryPoint(it, model->GetAtomArray()).X;
    SMLVec3d n = GetBoundaryPoint(it, model->GetAtomArray()).N;

    // Set up the medial coordinate and normal
    for(int j = 0; j < 3; j++)
      {
      sprintf(buffer, "X[%d,%d]", i, j);
      X[i][j] = p->AddVariable(buffer, x[j]);

      sprintf(buffer, "N[%d,%d]", i, j);
      N[i][j] = p->AddVariable(buffer, n[j]);
      }
    }

  // ------------------------------------------------------------------------
  // Configure the medial axis variables
  // ------------------------------------------------------------------------
  for(MedialAtomIterator it = model->GetAtomIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = it.GetIndex();
    MedialAtom &a = model->GetAtomArray()[i];
    for(int j = 0; j < 3; j++)
      {
      // Create the medial point
      sprintf(buffer, "M[%d,%d]", i, j);
      M[i][j] = p->AddVariable(buffer, a.X[j]);
      }

    sprintf(buffer, "R[%d]", i);
    // R[i] = p->AddVariable(buffer, a.R,
    //                       it.IsEdgeAtom() ? a.R : 0,
    //                       it.IsEdgeAtom() ? a.R : 32);
    R[i] = p->AddVariable(buffer, a.R, 0.1);
    }

  // Compute the U vectors
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = (int) it.GetIndex();
    SMLVec3d x = GetBoundaryPoint(it, model->GetAtomArray()).X;
    SMLVec3d n = GetBoundaryPoint(it, model->GetAtomArray()).N;

    // Set up the medial coordinate and normal
    for(int j = 0; j < 3; j++)
      U[i][j] = new BinaryDifference(p, X[i][j], M[it.GetAtomIndex()][j]);
    }

  // ------------------------------------------------------------------------
  // Configure the boundary and medial triangle area variables
  // ------------------------------------------------------------------------

  // Create and initialize expression arrays for boundary area element
  VarVec AeltX(model->GetNumberOfBoundaryPoints(), NULL);
  for(int i = 0; i < AeltX.size(); i++) AeltX[i] = new BigSum(p);

  // Create and initialize expression arrays for medial area element
  VarVec AeltM(model->GetNumberOfAtoms(), NULL);
  for(int i = 0; i < AeltM.size(); i++) AeltM[i] = new BigSum(p);

  // Process boundary triangles
  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Get the triangle index
    int k = trit.GetIndex();

    // Get the expression for (2A)^2
    Expression *areaSq = TriangleTwiceAreaSqr(p,
                                              X[trit.GetBoundaryIndex(0)],
                                              X[trit.GetBoundaryIndex(1)],
                                              X[trit.GetBoundaryIndex(2)]);

    // Get the actual area
    double a = sqrt(areaSq->Evaluate()) / 2.0;

    // Create the area variable
    sprintf(buffer, "taX[%d]", k);
    taX[k] = p->AddVariable(buffer, a, 0.0);

    // Accumulate the area element
    for(int j = 0; j < 3; j++)
      static_cast<BigSum *>(AeltX[trit.GetBoundaryIndex(j)])->AddSummand(taX[k]);

    // Create the constraint on the triangle area (2*A)^2 = areaSq
    p->AddConstraint(
          new BinaryDifference(
            p,
            new ScalarProduct(p, new Square(p, taX[k]), 2.0),
            areaSq), 0.0, 0.0);
    }

  // Process medial triangles
  for(MedialTriangleIterator trit = model->GetMedialTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Get the triangle index
    int k = trit.GetIndex();

    // Get the expression for (2A)^2
    Expression *areaSq = TriangleTwiceAreaSqr(p,
                                              M[trit.GetAtomIndex(0)],
                                              M[trit.GetAtomIndex(1)],
                                              M[trit.GetAtomIndex(2)]);

    // Get the actual area
    double a = sqrt(areaSq->Evaluate()) / 2.0;

    // Create the area variable
    sprintf(buffer, "taM[%d]", k);
    taM[k] = p->AddVariable(buffer, a, 0.0);

    // Accumulate the area element
    for(int j = 0; j < 3; j++)
      static_cast<BigSum *>(AeltM[trit.GetAtomIndex(j)])->AddSummand(taM[k]);

    // Create the constraint on the triangle area (2*A)^2 = areaSq
    p->AddConstraint(
          new BinaryDifference(
            p,
            new ScalarProduct(p, new Square(p, taM[k]), 2.0),
            areaSq), 0.0, 0.0);
    }

  // ------------------------------------------------------------------------
  // Configure the constraints on boundary normal and also compute curvatures
  // ------------------------------------------------------------------------

  // Create a LoopScheme for specifying normal vector constraints
  TriangleMesh *bmesh = model->GetIterationContext()->GetBoundaryMesh();
  LoopTangentScheme lts;
  lts.SetMesh(bmesh);

  // Expressions describing Xu, Xv, Nu, Nv at each point on the boundary
  VarVecArray Xd[] = {
    VarVecArray(model->GetNumberOfBoundaryPoints(), VarVec(3, NULL)),
    VarVecArray(model->GetNumberOfBoundaryPoints(), VarVec(3, NULL)) };

  VarVecArray Nd[] = {
    VarVecArray(model->GetNumberOfBoundaryPoints(), VarVec(3, NULL)),
    VarVecArray(model->GetNumberOfBoundaryPoints(), VarVec(3, NULL)) };

  VarVec curv_mean(model->GetNumberOfBoundaryPoints(), NULL);
  BigSum *objBending = new BigSum(p);
  BigSum *objSimpleBending = new BigSum(p);

  // Add the constraints relating each normal to the neighboring vertices
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = it.GetIndex();

    // Repeat for u and v directions
    for(int d = 0; d < 2; d++)
      {
      BigSum *Xdi[] = { new BigSum(p), new BigSum(p), new BigSum(p) };
      BigSum *Ndi[] = { new BigSum(p), new BigSum(p), new BigSum(p) };

      double wi = lts.GetOwnWeight(d, it.GetIndex());
      for(int j = 0; j < 3; j++)
        {
        Xdi[j]->AddSummand(new ScalarProduct(p, X[i][j], wi));
        Ndi[j]->AddSummand(new ScalarProduct(p, N[i][j], wi));
        }

      for(EdgeWalkAroundVertex walk(bmesh, it.GetIndex()); !walk.IsAtEnd(); ++walk)
        {
        double wij = lts.GetNeighborWeight(d, walk);
        for(int j = 0; j < 3; j++)
          {
          Xdi[j]->AddSummand(new ScalarProduct(p, X[walk.MovingVertexId()][j], wij));
          Ndi[j]->AddSummand(new ScalarProduct(p, N[walk.MovingVertexId()][j], wij));
          }
        }

      // Store these expressions
      for(int j = 0; j < 3; j++)
        {
        Xd[d][i][j] = Xdi[j];
        Nd[d][i][j] = Ndi[j];
        }

      // Create the constrait Xu . N = 0
      Expression *constrNormXu = DotProduct(p, Xd[d][i], N[i]);

      // Add the constraint to the problem
      p->AddConstraint(constrNormXu, 0, 0);
      }

    // While we are here, add the constraint on the normal being norm 1
    Expression *constrNormMag =
        new TernaryGradientMagnitudeSqr(p,
                                        N[it.GetIndex()][0],
                                        N[it.GetIndex()][1],
                                        N[it.GetIndex()][2]);
    p->AddConstraint(constrNormMag, 1.0, 1.0);

    // Compute the first/second fundamental form at the point
    Expression *FF1[2][2], *FF2[2][2];
    for(int q = 0; q < 2; q++)
      {
      for(int r = 0; r < 2; r++)
        {
        FF1[q][r] = DotProduct(p, Xd[q][i], Xd[r][i]);
        FF2[q][r] = new Negation(p, DotProduct(p, Xd[q][i], Nd[r][i]));
        }
      }

    // Compute the mean curvature
    curv_mean[i] =
        new BinaryFraction(
          p,
          new BinaryDifference(
            p,
            new BinarySum(
              p,
              new BinaryProduct(p, FF2[0][0], FF1[1][1]),
              new BinaryProduct(p, FF2[1][1], FF1[0][0])),
            new BinarySum(
              p,
              new BinaryProduct(p, FF2[0][1], FF1[1][0]),
              new BinaryProduct(p, FF2[1][0], FF1[0][1]))),
          new ScalarProduct(
            p,
            new BinaryDifference(
              p,
              new BinaryProduct(p, FF1[0][0], FF1[1][1]),
              new Square(p, FF1[0][1])),
            2.0)) ;

    objBending->AddSummand(new Square(p, new BinaryDifference(
                                        p,
                                        new BinarySum(
                                          p,
                                          new BinaryProduct(p, FF2[0][0], FF1[1][1]),
                                          new BinaryProduct(p, FF2[1][1], FF1[0][0])),
                                        new BinarySum(
                                          p,
                                          new BinaryProduct(p, FF2[0][1], FF1[1][0]),
                                          new BinaryProduct(p, FF2[1][0], FF1[0][1])))));

    // Compute the simple bending objective
    BigSum *Xj[] = { new BigSum(p), new BigSum(p), new BigSum(p) };
    EdgeWalkAroundVertex walk(bmesh, it.GetIndex());
    walk.Valence();

    for(; !walk.IsAtEnd(); ++walk)
      {
      int k = walk.MovingVertexId();
      for(int j = 0; j < 3; j++)
        {
        Xj[j]->AddSummand(X[k][j]);
        }
      }

    double scale = 1.0 / walk.Valence();
    objSimpleBending->AddSummand(new TernarySum(p,
                                   new Square(
                                     p, new BinaryDifference(
                                       p,
                                       new ScalarProduct(p, Xj[0], scale),
                                       X[i][0])),
                                   new Square(
                                     p, new BinaryDifference(
                                       p,
                                       new ScalarProduct(p, Xj[1], scale),
                                       X[i][1])),
                                   new Square(
                                     p, new BinaryDifference(
                                       p,
                                       new ScalarProduct(p, Xj[2], scale),
                                       X[i][2]))));
    }

  // ------------------------------------------------------------------------
  // Add the actual medial constraints
  // ------------------------------------------------------------------------

  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int iBnd = it.GetIndex();
    int iAtom = it.GetAtomIndex();

    // Code up X - r * N - M = 0
    for(int j = 0; j < 3; j++)
      {
      Expression *constMedial =
          new BinaryDifference(p,
                               X[iBnd][j],
                               new BinarySum(p,
                                             M[iAtom][j],
                                             new BinaryProduct(p,
                                                               R[iAtom],
                                                               N[iBnd][j])));
      p->AddConstraint(constMedial, 0, 0);
      }
    }


  // ------------------------------------------------------------------------
  // Define the objective on the basis
  // ------------------------------------------------------------------------

  // Define a basis for the surface
  int nBasis = 20;
  MeshBasisCoefficientMapping basismap(
        model->GetIterationContext()->GetBoundaryMesh(), nBasis, 3);

  // Create the coefficient variables
  VarVecArray XC(nBasis, VarVec(3, NULL));
  for(int i = 0; i < nBasis; i++)
    {
    for(int j = 0; j < 3; j++)
      {
      XC[i][j] = p->AddVariable("XC", 0.0);
      }
    }

  // Define the objective on the basis
  BigSum *objBasisResidual = new BigSum(p);
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int iBnd = it.GetIndex();
    SMLVec3d Xfixed = GetBoundaryPoint(it, model->GetAtomArray()).X;

    for(int j = 0; j < 3; j++)
      {
      BigSum *xfit = new BigSum(p);

      xfit->AddSummand(new Constant(p, Xfixed[j]));
      for(int i = 0; i < nBasis; i++)
        {
        xfit->AddSummand(new ScalarProduct(p, XC[i][j],
                                           basismap.GetBasisComponent(i, iBnd)));
        }

      // Xfit is the approximation of X using the basis
      objBasisResidual->AddSummand(
            new Square(p, new BinaryDifference(p, xfit, X[iBnd][j])));
      }
    }


  // ------------------------------------------------------------------------
  // Create a total volume objective -
  // ------------------------------------------------------------------------
  BigSum *objVolume = new BigSum(p);

  // Measure the wedge volume for each boundary triangle
  VarVec wedgeVol(model->GetNumberOfBoundaryTriangles(), NULL);

  // Add the Jacobian constraints - all tetrahedra must have positive volume
  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Get the boundary vertices
    VarVec x[] = { X[trit.GetBoundaryIndex(0)],
                   X[trit.GetBoundaryIndex(1)],
                   X[trit.GetBoundaryIndex(2)] };

    // And the medial vertices
    VarVec m[] = { M[trit.GetAtomIndex(0)],
                   M[trit.GetAtomIndex(1)],
                   M[trit.GetAtomIndex(2)] };

    // There are several ways to cut a wedge into tetras, we choose one
    Expression *c1 = TetraHedronVolume(p, model, m[2], x[0], x[1], x[2]);
    Expression *c2 = TetraHedronVolume(p, model, m[1], x[0], x[1], m[2]);
    Expression *c3 = TetraHedronVolume(p, model, m[2], x[0], m[0], m[1]);

    // printf("Tetra Vols: %f, %f, %f\n", c1->Evaluate(), c2->Evaluate(), c3->Evaluate());

    // Save the total wedge volume for use in integration
    wedgeVol[trit.GetIndex()] = new TernarySum(p, c1, c2, c3);

    // Each tetra should have positive volume
    /*
    p->AddConstraint(c1, 0.1, 40);
    p->AddConstraint(c2, 0.1, 40);
    p->AddConstraint(c3, 0.1, 40);
    */

    // Total volume integral
    objVolume->AddSummand(wedgeVol[trit.GetIndex()]);
    }

  // ------------------------------------------------------------------------
  // Create the medial/boundary Jacobian constraint (normals point in same direction)
  // ------------------------------------------------------------------------

  // Lower bound for Nm . Nx / (Nm . Nm)
  double constJacFact = 0.1;

  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Get the boundary vertices
    VarVec x[] = { X[trit.GetBoundaryIndex(0)],
                   X[trit.GetBoundaryIndex(1)],
                   X[trit.GetBoundaryIndex(2)] };

    // And the medial vertices
    VarVec m[] = { M[trit.GetAtomIndex(0)],
                   M[trit.GetAtomIndex(1)],
                   M[trit.GetAtomIndex(2)] };

    // Get the expression for the triangle normal
    VarVec xu = VectorApplyPairwise<BinaryDifference>(p, x[1], x[0]);
    VarVec xv = VectorApplyPairwise<BinaryDifference>(p, x[2], x[0]);
    VarVec Nx = CrossProduct(p, xu, xv);

    VarVec mu = VectorApplyPairwise<BinaryDifference>(p, m[1], m[0]);
    VarVec mv = VectorApplyPairwise<BinaryDifference>(p, m[2], m[0]);
    VarVec Nm = CrossProduct(p, mu, mv);

    // Get the dot products
    Expression *NmNx = DotProduct(p, Nx, Nm);
    Expression *NmNm = DotProduct(p, Nm, Nm);

    // Get the Jacobian constraint expression
    Expression *con = new BinaryDifference(p, NmNx,
                                           new ScalarProduct(p, NmNm, constJacFact));

    // Set the constraint
    p->AddConstraint(con, 0.0, ConstrainedNonLinearProblem::UBINF);

    // We can also constrain the size of the medial triangle from below
    p->AddConstraint(NmNm, 0.1, ConstrainedNonLinearProblem::UBINF);

    // Print the ratio
    // printf("JACCON: %f %f \n", NmNm->Evaluate(), NmNx->Evaluate() / NmNm->Evaluate());

    }

  // ------------------------------------------------------------------------
  // Create the MIB constraint
  // ------------------------------------------------------------------------
  TriangleMesh *mesh = model->GetIterationContext()->GetBoundaryMesh();
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int iBnd = it.GetIndex();
    int iAtom = it.GetAtomIndex();

    for(EdgeWalkAroundVertex walk(mesh, iBnd); !walk.IsAtEnd(); ++walk)
      {
      int k = walk.MovingVertexId();
      Expression *distsq = new TernarySum(
            p,
            new Square(p, new BinaryDifference(p, M[iAtom][0], X[k][0])),
            new Square(p, new BinaryDifference(p, M[iAtom][1], X[k][1])),
            new Square(p, new BinaryDifference(p, M[iAtom][2], X[k][2])));
      p->AddConstraint(
            new BinaryDifference(p, distsq, new Square(p, R[iAtom])),
            0, ConstrainedNonLinearProblem::UBINF);
      }
    }

  // ------------------------------------------------------------------------
  // Construct the first part of the objective function
  // ------------------------------------------------------------------------
  BigSum *objSqDist = new BigSum(p);
  for(int i = 0; i < X.size(); i++)
    {
    for(int j = 0; j < 3; j++)
      {
      objSqDist->AddSummand(
            new Square(p,
                       new BinaryDifference(p, X[i][j],
                                            new Constant(p, targetPoint[i][j]))));
      }
    }

  // ------------------------------------------------------------------------
  // Construct an opposite objective: distance from a selected set of points
  // to the closest point on the medial mesh
  // ------------------------------------------------------------------------
  std::vector<PointMatch> meshToModel = find_closest_to_mesh(target, model, 1000);
  BigSum *objRecipSqDist = new BigSum(p);
  for(int i = 0; i < meshToModel.size(); i++)
    {
    SMLVec3d xMesh = meshToModel[i].second;
    int iModel = meshToModel[i].first;
    for(int j = 0; j < 3; j++)
      {
      objRecipSqDist->AddSummand(
            new Square(p,
                       new BinaryDifference(p, X[iModel][j],
                                            new Constant(p, xMesh[j]))));
      }
    }

  // ------------------------------------------------------------------------
  // Construct the total surface area objective
  // ------------------------------------------------------------------------
  BigSum *objSurfArea = new BigSum(p);
  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Add the area to the objective
    objSurfArea->AddSummand(taX[trit.GetIndex()]);
    }

  // ------------------------------------------------------------------------
  // Add a constraint on triangle aspect ratios
  // ------------------------------------------------------------------------

  // To create this constraint, we need to first create a variable for the
  // length of every edge in the mesh.
  typedef std::pair<int, int> EdgeType;
  typedef std::map<EdgeType, Expression *> EdgeLenMap;
  EdgeLenMap edgeLenMap, edgeLenMapSqr;

  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Look across from each vertex in the triangle
    VarVec edge(3, NULL);
    VarVec edgesq(3, NULL);

    for(int j = 0; j < 3; j++)
      {
      int k1 = trit.GetBoundaryIndex((j + 1) % 3);
      int k2 = trit.GetBoundaryIndex((j + 2) % 3);
      EdgeType e = std::make_pair(std::min(k1, k2), std::max(k1, k2));

      EdgeLenMap::const_iterator eit = edgeLenMap.find(e);
      if(eit == edgeLenMap.end())
        {
        // Compute the squared length of the egde
        Expression *len2 = new TernaryGradientMagnitudeSqr(
              p,
              new BinaryDifference(p, X[k1][0],X[k2][0]),
              new BinaryDifference(p, X[k1][1],X[k2][1]),
              new BinaryDifference(p, X[k1][2],X[k2][2]));

        // Create the variable for this edge
        edge[j] = p->AddVariable("len", sqrt(len2->Evaluate()), 0);
        edgesq[j] = new Square(p, edge[j]);

        // Constrain the variable to the actual length
        p->AddConstraint(new BinaryDifference(p, len2, edgesq[j]), 0, 0);

        // Store the edge length variable
        edgeLenMap[e] = edge[j];
        edgeLenMapSqr[e] = edgesq[j];
        }
      else
        {
        edge[j] = eit->second;
        edgesq[j] = edgeLenMapSqr[e];
        }
      }

    // Create constraints based on the cosine rule
    double min_angle = vnl_math::pi * 20 / 180;
    double max_cos = cos(min_angle);

    for(int j = 0; j < 3; j++)
      {
      Expression *c = edge[j], *a = edge[(j+1) % 3], *b = edge[(j+2) % 3];
      Expression *c2 = edgesq[j], *a2 = edgesq[(j+1) % 3], *b2 = edgesq[(j+2) % 3];

      // (c^2 - 2 * a * b * cosTheta) - (a^2 + b^2)
      Expression *constr = new BinaryDifference(
            p,
            new BinarySum(
              p, c2, new ScalarProduct(p, new BinaryProduct(p, a, b), 2 * max_cos)),
            new BinarySum(p, a2, b2));
      p->AddConstraint(constr, 0.0, ConstrainedNonLinearProblem::UBINF);
      }
    }

  /*
  double constEdgeRatio = 0.1;
  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Compute the squared edge length opposite each vertex
    Expression *edge[3];
    for(int j = 0; j < 3; j++)
      {
      int k1 = trit.GetBoundaryIndex((j + 1) % 3);
      int k2 = trit.GetBoundaryIndex((j + 2) % 3);
      edge[j] = new TernaryGradientMagnitudeSqr(p,
                                                new BinaryDifference(p, X[k1][0],X[k2][0]),
                                                new BinaryDifference(p, X[k1][1],X[k2][1]),
                                                new BinaryDifference(p, X[k1][2],X[k2][2]));
      }

    // Create the constraint for each edge
    for(int j = 0; j < 3; j++)
      {
      p->AddConstraint(
            new BinaryDifference(p, edge[j],
                                 new ScalarProduct(p, edge[(j+1) % 3], constEdgeRatio)),
            0.0, ConstrainedNonLinearProblem::UBINF);
      }
    }
  */

  // ------------------------------------------------------------------------
  // Derive an image match objective
  // ------------------------------------------------------------------------
  VarVecArray sampleX;
  VarVec sampleF;
  Expression *objObjectIntegral, *objVolumeIntegral;

  CreateOverlapObjective(p, model, &interp_new, X, M, U,
                         &objObjectIntegral, &objVolumeIntegral,
                         sampleX, sampleF);

  /*
  CreateAreaElementScaledOverlapObjective(p, model, &interp_new, X, M, R,
                                          AeltX, AeltM,
                                          &objObjectIntegral, &objVolumeIntegral,
                                          sampleX, sampleF);
  */


  /*
  VolumeOverlapEnergyTerm compterm(model, image, 3);
  SolutionData sd(model->GetIterationContext(), model->GetAtomArray());
  sd.ComputeIntegrationWeights();
  compterm.ComputeEnergy(&sd);

  printf("*** Overlap term report ***\n");
  printf("ObjInt = %f, VolInt = %f, Samples = %d\n",
         objObjectIntegral->Evaluate(),
         objVolumeIntegral->Evaluate(),
         sampleF.size());
  compterm.PrintReport(std::cout);
  */




  /*

  // Integrate match over wedges. The factor of 18 here accounts for the
  // fact that wedgeVol is scaled by six
  Expression *sampleScaleFactor = new Constant(
        p, model->GetNumberOfBoundaryPoints() * 1.0 / (18 * sampleF.size()));

  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    Expression *meanWedgeMatch = new BinaryProduct(
          p,
          sampleScaleFactor,
          new TernarySum(
            p,
            spokeMatch[trit.GetBoundaryIndex(0)],
            spokeMatch[trit.GetBoundaryIndex(1)],
            spokeMatch[trit.GetBoundaryIndex(2)]));

    objImageMatch->AddSummand(
          new BinaryProduct(p, meanWedgeMatch, wedgeVol[trit.GetIndex()]));
    }
    */

  // Compute a Dice-like objective value
  double volTarget = image->ComputeObjectVolume();
  Expression *objDice = new BinaryDifference(
        p, new Constant(p, 1.0),
        new ScalarProduct(p, objObjectIntegral, 1.0 / volTarget));

  // As an initial objective, we just want to resolve all the constraints
  BigSum *objDisplacement = new BigSum(p);
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = it.GetIndex();
    objDisplacement->AddSummand(
          new TernaryGradientMagnitudeSqr(
            p,
            new BinaryDifference(p, X[i][0], new Constant(p, X[i][0]->Evaluate())),
            new BinaryDifference(p, X[i][1], new Constant(p, X[i][1]->Evaluate())),
            new BinaryDifference(p, X[i][2], new Constant(p, X[i][2]->Evaluate()))));
    }

  // ------------------------------------------------------------------------
  // Derive the final objective
  // ------------------------------------------------------------------------
  double scaleRecip = model->GetNumberOfBoundaryPoints() * 1.0 / 1000;
  BigSum *obj = new BigSum(p);
  // obj->AddSummand(objSqDist);
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, scaleRecip), objRecipSqDist));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 20), objImageMatch));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 0.01), objSurfArea));


  obj->AddSummand(new ScalarProduct(p, objDice, 2000));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 0.5), objSimpleBending));

  // Set of objectives for fitting the model to itself
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 1), objDisplacement));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 0.1), objSimpleBending));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 0.02), objSurfArea));

  obj->AddSummand(new ScalarProduct(p, objBasisResidual, 1.0));

  // ------------------------------------------------------------------------
  // Create diffeomorphic constraints for the mesh
  // ------------------------------------------------------------------------
  // CreateTetgenMesh(model, image, X, p);

  // ------------------------------------------------------------------------
  // Create a kernel for smoothing the gradient
  // ------------------------------------------------------------------------
  /*
  typedef ImmutableSparseMatrix<double> SparseMatrix;
  SparseMatrix::STLSourceType kSource(p->GetNumberOfVariables());

  double alpha = 0.2;
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int iBnd = it.GetIndex();
    int iAtom = it.GetAtomIndex();

    for(int j = 0; j < 3; j++)
      {
      int idx_i = static_cast<Variable *>(X[iBnd][j])->GetIndex();
      for(EdgeWalkAroundVertex walk(mesh, iBnd); !walk.IsAtEnd(); ++walk)
        {
        int k = walk.MovingVertexId();
        int idx_j = static_cast<Variable *>(X[k][j])->GetIndex();
        SparseMatrix::STLEntryType entry;
        entry.first = idx_j;
        entry.second = alpha / walk.Valence();
        kSource[idx_i].push_back(entry);
        }
      kSource[idx_i].push_back(make_pair(idx_i, 1.0 - alpha));
      }
    }

  // Initialize the rest of the matrix to identity
  for(int i = 0; i < kSource.size(); i++)
    {
    if(kSource[i].size() == 0)
      kSource[i].push_back(make_pair(i, 1.0));
    }

  // Create a sparse matrix and multiply by itself a couple times
  SparseMatrix K;
  K.SetFromSTL(kSource, kSource.size());
  for(int q = 0; q < 2; q++)
    SparseMatrix::Multiply(K, K, K);
  p->SetGradientSmoothingKernel(K);
  */


  // Save the sample and image values at samples
  SaveSamples(sampleX, sampleF, "samples_before.vtk");

  // Plot the gradient of the problem
  SaveGradient(p, model, X, obj, "grad_obj_before.vtk");
  SaveGradient(p, model, X, objDice, "grad_obj_dice_before.vtk");
  SaveGradient(p, model, X, objSimpleBending, "grad_obj_bend_before.vtk");


  // Evaluate the objective
  std::cout << "MSD to target: " << objSqDist->Evaluate() << std::endl;
  std::cout << "MSD to model: " << scaleRecip * objRecipSqDist->Evaluate() << std::endl;
  std::cout << "Surface area: " << objSurfArea->Evaluate() << std::endl;
  std::cout << "Model volume: " << objVolume->Evaluate() / 6 << std::endl;
  std::cout << "Target volume: " << volTarget << std::endl;
  std::cout << "Image match: " << objObjectIntegral->Evaluate() << std::endl;
  std::cout << "Volume integral: " << objVolumeIntegral->Evaluate() << std::endl;
  std::cout << "Bending energy: " << objBending->Evaluate() << std::endl;
  std::cout << "Simple bending energy: " << objSimpleBending->Evaluate() << std::endl;
  std::cout << "Dice objective: " << 1- objDice->Evaluate() << std::endl;
  std::cout << "Displacement objective: " << objDisplacement->Evaluate() << std::endl;
  std::cout << "Total objective: " << obj->Evaluate() << std::endl;

  // Configure the problem
  p->SetObjective(obj);
  p->SetupProblem(true);
  p->ClearDerivativeCaches();

  // Test some derivatives;
  // DerivativeTest(p, 100);

  // Stats on how many expressions we have
  const Problem::ExpressionSet &expr = p->GetChildExpressions();
  int nConst = 0;
  for(Problem::ExpressionSet::const_iterator it = expr.begin(); it != expr.end(); ++it)
    {
    Expression *e = (*it);
    if(dynamic_cast<Constant *>(e))
      nConst++;
    }

  printf("Expressions: %d, constants: %d\n", expr.size(), nConst);


  // Solve the problem
  SmartPtr<IPOptProblemInterface> ip = new IPOptProblemInterface(p);

  // Set up the IPopt problem
  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-8);
  // app->Options()->SetNumericValue("mu_init", 1e-3);
  // app->Options()->SetNumericValue("mu_target", 1e-5);
  // app->Options()->SetStringValue("mu_strategy", "adaptive");
  // app->Options()->SetStringValue("output_file", "ipopt.out");

  app->Options()->SetIntegerValue("max_iter", 200);
  // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  // app->Options()->SetStringValue("derivative_test", "second-order");
  // app->Options()->SetStringValue("derivative_test_print_all", "yes");

  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
    }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(GetRawPtr(ip));

  std::cout << "MSD to target: " << objSqDist->Evaluate() << std::endl;
  std::cout << "MSD to model: " << scaleRecip * objRecipSqDist->Evaluate() << std::endl;
  std::cout << "Surface area: " << objSurfArea->Evaluate() << std::endl;
  std::cout << "Model volume: " << objVolume->Evaluate() / 6 << std::endl;
  std::cout << "Target volume: " << volTarget << std::endl;
  std::cout << "Image match: " << objObjectIntegral->Evaluate() << std::endl;
  std::cout << "Volume integral: " << objVolumeIntegral->Evaluate() << std::endl;
  std::cout << "Bending energy: " << objBending->Evaluate() << std::endl;
  std::cout << "Simple bending energy: " << objSimpleBending->Evaluate() << std::endl;
  std::cout << "Dice objective: " << 1- objDice->Evaluate() << std::endl;
  std::cout << "Displacement objective: " << objDisplacement->Evaluate() << std::endl;
  std::cout << "Residual objective: " << objBasisResidual->Evaluate() << std::endl;
  std::cout << "Total objective: " << obj->Evaluate() << std::endl;


  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
    }
  else {
    printf("\n\n*** The problem FAILED!\n");
    }


  // Save the mesh for output
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    for(int j = 0; j < 3; j++)
      {
      int i = it.GetIndex();
      BoundaryAtom &B = GetBoundaryPoint(it, model->GetAtomArray());
      B.X[j] = X[i][j]->Evaluate();
      B.N[j] = N[i][j]->Evaluate();
      B.curv_mean = curv_mean[i]->Evaluate();
      }
    }

  // Configure the medial axis variables
  for(MedialAtomIterator it = model->GetAtomIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = it.GetIndex();
    MedialAtom &a = model->GetAtomArray()[i];
    for(int j = 0; j < 3; j++)
      {
      a.X[j] = M[i][j]->Evaluate();
      }
    a.R = R[i]->Evaluate();
    }

  // Save the resulting meshes
  ExportMedialMeshToVTK(model, NULL, "test.med.vtk");
  ExportBoundaryMeshToVTK(model, NULL, "test.bnd.vtk");

  // Save the sample and image values at samples
  SaveSamples(sampleX, sampleF, "samples_after.vtk");

  // Plot the gradient of the problem
  SaveGradient(p, model, X, obj, "grad_obj_after.vtk");
  SaveGradient(p, model, X, objDice, "grad_obj_dice_after.vtk");
  SaveGradient(p, model, X, objBasisResidual, "grad_obj_residual_after.vtk");

  vtkPoints *out_pts = vtkPoints::New();
  out_pts->Allocate(model->GetNumberOfBoundaryPoints());
  vtkPolyData *out_poly = vtkPolyData::New();
  out_poly->Allocate(model->GetNumberOfBoundaryTriangles());
  out_poly->SetPoints(out_pts);
  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = it.GetIndex();
    out_pts->InsertNextPoint(
          X[i][0]->Evaluate(), X[i][1]->Evaluate(), X[i][2]->Evaluate());
    }

  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    vtkIdType ids[] = {trit.GetBoundaryIndex(0),
                       trit.GetBoundaryIndex(1),
                       trit.GetBoundaryIndex(2)};

    out_poly->InsertNextCell(VTK_TRIANGLE, 3, ids);
    }

  out_poly->BuildCells();
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(out_poly);
  writer->SetFileName("testout.vtk");
  writer->Update();


  // Save the mesh for output
  vtkPoints *out_mpts = vtkPoints::New();
  out_mpts->Allocate(model->GetNumberOfAtoms());
  vtkPolyData *out_mpoly = vtkPolyData::New();
  out_mpoly->Allocate(model->GetNumberOfTriangles());
  out_mpoly->SetPoints(out_mpts);
  for(MedialAtomIterator it = model->GetAtomIterator();
      !it.IsAtEnd(); ++it)
    {
    int i = it.GetIndex();
    out_mpts->InsertNextPoint(
          M[i][0]->Evaluate(), M[i][1]->Evaluate(), M[i][2]->Evaluate());
    }

  for(MedialTriangleIterator trit = model->GetMedialTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    vtkIdType ids[] = {trit.GetAtomIndex(0),
                       trit.GetAtomIndex(1),
                       trit.GetAtomIndex(2)};

    out_mpoly->InsertNextCell(VTK_TRIANGLE, 3, ids);
    }

  out_mpoly->BuildCells();
  vtkPolyDataWriter *mwriter = vtkPolyDataWriter::New();
  mwriter->SetInput(out_mpoly);
  mwriter->SetFileName("testmout.vtk");
  mwriter->Update();


  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.
  delete p;

  return (int) status;
}
