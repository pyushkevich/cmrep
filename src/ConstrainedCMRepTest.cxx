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
#include "vtkPolyDataReader.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"


#include "tetgen.h"

#include "vtkPolyData.h"
#include "vtkCellLocator.h"
#include "vtkSmartPointer.h"
#include <vtkPointLocator.h>

#include <vector>
#include <map>
#include <utility>
#include "itk_to_nifti_xform.h"

#include "ConstrainedCMRepObjectives.h"
#include "itksys/SystemTools.hxx"

/**
  ---------
  TODO LIST
  ---------

  * Add a constraint on the dihedral angle of medial triangles

*/
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

class ClosestPointMatcher
{
public:

  struct MatchLocation
  {
    // The point on the target mesh
    SMLVec3d xTarget;

    // The triangle that is closest to the point on the target mesh
    int iTriangle;

    // The barycentric coordinates of the closest point on the triangle
    SMLVec3d xBary;
  };

  // A pair consisting of a
  typedef std::pair<MatchLocation, SMLVec3d> PointMatch;

  ClosestPointMatcher(vtkPolyData *target, int nClusterDivisions = 32);

  std::vector<SMLVec3d> FindClosestToTarget(VarVecArray &X);

  std::vector<MatchLocation> FindClosestToSource(TriangleMesh *mesh, VarVecArray &X);

  int GetNumberOfTargetPointsUsed()
  { return m_ReducedTarget->GetNumberOfPoints(); }

protected:
  vtkSmartPointer<vtkPolyData> m_Target;
  vtkSmartPointer<vtkCellLocator> m_TargetLocator;
  vtkSmartPointer<vtkPoints> m_ReducedTarget;
};

#include <vtkQuadricClustering.h>

ClosestPointMatcher
::ClosestPointMatcher(vtkPolyData *target, int nClusterDivisions)
{
  // Store the target
  m_Target = target;

  // Create a locator
  m_TargetLocator = vtkSmartPointer<vtkCellLocator>::New();
  m_TargetLocator->SetDataSet(m_Target);
  m_TargetLocator->CacheCellBoundsOn();
  m_TargetLocator->BuildLocator();

  // Need the bounds on the target
  target->ComputeBounds();

  // Set up the clustering
  vtkSmartPointer<vtkQuadricClustering> clu =
      vtkSmartPointer<vtkQuadricClustering>::New();
  clu->SetInputData(target);
  clu->SetDivisionOrigin(target->GetCenter());
  double spacing = target->GetLength() / nClusterDivisions;
  clu->SetDivisionSpacing(spacing, spacing, spacing);
  clu->Update();

  // Get the reduced target
  m_ReducedTarget = clu->GetOutput()->GetPoints();

  // Save the samples
  vtkSmartPointer<vtkPolyDataWriter> wr =
      vtkSmartPointer<vtkPolyDataWriter>::New();
  wr->SetInputConnection(clu->GetOutputPort());
  wr->SetFileName("clusty.vtk");
  wr->Update();
}

std::vector<SMLVec3d>
ClosestPointMatcher::FindClosestToTarget(VarVecArray &X)
{
  // Output vector
  std::vector<SMLVec3d> cp(X.size());

  // Compute thickness values
  for(int i = 0; i < X.size(); i++)
    {
    SMLVec3d xi = gnlp::VectorEvaluate(X[i]);

    double xs[3], d2, d;
    int subid;
    vtkIdType cellid;

    m_TargetLocator->FindClosestPoint(xi.data_block(), xs, cellid, subid, d2);

    cp[i] = SMLVec3d(xs);
    }

  return cp;
}

std::vector<ClosestPointMatcher::MatchLocation>
ClosestPointMatcher::FindClosestToSource(TriangleMesh *mesh, VarVecArray &X)
{
  // Create a VTK points object
  vtkSmartPointer<vtkPoints> out_pts = vtkSmartPointer<vtkPoints>::New();
  out_pts->Allocate(X.size());

  for(int i = 0; i < X.size(); i++)
    {
    out_pts->InsertNextPoint(
          X[i][0]->Evaluate(), X[i][1]->Evaluate(), X[i][2]->Evaluate());
    }

  // Create the polydata
  vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
  poly->SetPoints(out_pts);
  poly->Allocate(mesh->triangles.size());

  for(int i = 0; i < mesh->triangles.size(); i++)
    {
    Triangle &t = mesh->triangles[i];
    vtkIdType v[] = {(vtkIdType) t.vertices[0], (vtkIdType) t.vertices[1], (vtkIdType) t.vertices[2]};
    poly->InsertNextCell(VTK_TRIANGLE, 3, v);
    }

  // Build everything
  poly->BuildCells();

  // Create locator for finding closest points
  vtkSmartPointer<vtkCellLocator> locator = vtkSmartPointer<vtkCellLocator>::New();
  locator->SetDataSet(poly);
  locator->BuildLocator();

  // Create a point set for debugging
  vtkSmartPointer<vtkFloatArray> vMatch = vtkSmartPointer<vtkFloatArray>::New();
  vMatch->SetNumberOfComponents(3);
  vMatch->Allocate(3 * m_ReducedTarget->GetNumberOfPoints());
  vMatch->SetName("ClosestPoint");

  // Sample points from the target mesh
  std::vector<MatchLocation> result;
  for(int i = 0; i < m_ReducedTarget->GetNumberOfPoints(); i++)
    {
    MatchLocation loc;
    loc.xTarget.set(m_ReducedTarget->GetPoint(i));

    double xs[3], d2, d;
    int subid;
    vtkIdType cellid;

    // Find the closest point
    locator->FindClosestPoint(loc.xTarget.data_block(), xs, cellid, subid, d2);

    // Solve a system for the barycentric coordinates
    Triangle &T = mesh->triangles[cellid];
    SMLVec3d A = VectorEvaluate(X[T.vertices[0]]);
    SMLVec3d B = VectorEvaluate(X[T.vertices[1]]);
    SMLVec3d C = VectorEvaluate(X[T.vertices[2]]);

    vnl_matrix<double> W(3, 2);
    W.set_column(0, B-A);
    W.set_column(1, C-A);
    SMLVec3d v = SMLVec3d(xs) - A;

    vnl_matrix<double> WtW = W.transpose() * W;
    vnl_vector<double> Wtv = W.transpose() * v;
    vnl_vector<double> q = vnl_inverse(WtW) * Wtv;

    // The barycentric coordinates are (1-q1-q2, q1, q2)
    loc.iTriangle = cellid;
    loc.xBary[0] = 1 - (q[0] + q[1]);
    loc.xBary[1] = q[0];
    loc.xBary[2] = q[1];

    result.push_back(loc);

    SMLVec3d xss = loc.xBary[0] * A + loc.xBary[1] * B + loc.xBary[2] * C - loc.xTarget;
    vMatch->InsertNextTuple3(xss[0], xss[1], xss[2]);
    }


  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  pd->SetPoints(m_ReducedTarget);
  pd->GetPointData()->AddArray(vMatch);

  // Save the samples
  vtkSmartPointer<vtkPolyDataWriter> wr =
      vtkSmartPointer<vtkPolyDataWriter>::New();
  wr->SetInputData(pd);
  wr->SetFileName("pointmatch.vtk");
  wr->Update();

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
    Problem *p,
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


double GetCentralDifference(Problem *p, Expression *ex, Variable *v, double delta=1e-5)
{
  double val = v->Evaluate();

  // Do the perturbation
  ex->MakeTreeDirty();
  v->SetValue(val + delta);
    double f2 = ex->Evaluate();

  ex->MakeTreeDirty();
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
  double dCentralDiff = GetCentralDifference(p, ex, v);

  printf("D[%10s,%10s,%d]: %12.8f  %12.8f  %12.8f\n",
         nickname, v->GetName().c_str(), order,
         dAnalytic, dCentralDiff, std::fabs(dAnalytic-dCentralDiff));

  // Test higher order derivatives
  TestExpressionRandomDerivative(p, pd, nickname, order-1);
}

void DerivativeTest(ConstrainedNonLinearProblem *p, int nTests)
{
  p->MakeChildrenDirty();
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

  p->MakeChildrenDirty();
}




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
  writer->SetInputData(poly);
  writer->SetFileName(filename);
  writer->Update();

}

void SaveGradient(
    ConstrainedNonLinearProblem *p,
    std::vector<std::vector<Expression *> > X,
    Expression *f,
    const char *filename)
{
  vtkPoints *pts = vtkPoints::New();
  pts->Allocate(X.size());

  vtkFloatArray *arr = vtkFloatArray::New();
  arr->SetNumberOfComponents(3);
  arr->Allocate(X.size());
  arr->SetName("Gradient");

  for(int i = 0; i < X.size(); i++)
    {
    pts->InsertNextPoint(X[i][0]->Evaluate(),
                         X[i][1]->Evaluate(),
                         X[i][2]->Evaluate());
    Expression *dx = p->GetPartialDerivative(f, (Variable *)X[i][0]);
    Expression *dy = p->GetPartialDerivative(f, (Variable *)X[i][1]);
    Expression *dz = p->GetPartialDerivative(f, (Variable *)X[i][2]);
    arr->InsertNextTuple3(dx ? dx->Evaluate() : 0,
                          dy ? dy->Evaluate() : 0,
                          dz ? dz->Evaluate() : 0);
    }

  vtkPolyData *poly = vtkPolyData::New();
  poly->SetPoints(pts);
  poly->GetPointData()->SetScalars(arr);

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(poly);
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
        TetraHedronVolume(problem, Y[tet[0]], Y[tet[1]], Y[tet[2]], Y[tet[3]]);

    double tv = vol->Evaluate();

    problem->AddConstraint(vol, "TETVOL", 0.1 * tv, 100 * tv);
    }
}







void SaveCircumcenterMesh(VarVecArray &CC, VarVec &CR, VarVecArray &CCBC)
{

  vtkPoints *out_pts = vtkPoints::New();
  out_pts->Allocate(CC.size());

  vtkPolyData *out_poly = vtkPolyData::New();
  out_poly->SetPoints(out_pts);

  vtkFloatArray *arrRad = vtkFloatArray::New();
  arrRad->SetNumberOfComponents(1);
  arrRad->Allocate(CC.size());
  arrRad->SetName("Radius");

  vtkFloatArray *arrBC = vtkFloatArray::New();
  arrBC->SetNumberOfComponents(3);
  arrBC->Allocate(CC.size());
  arrBC->SetName("BC");

  for(int i = 0; i < CC.size(); i++)
    {
    out_pts->InsertNextPoint(
          CC[i][0]->Evaluate(), CC[i][1]->Evaluate(), CC[i][2]->Evaluate());

    arrRad->InsertNextTuple1(CR[i]->Evaluate());
    arrBC->InsertNextTuple3(
          CCBC[i][0]->Evaluate(),
          CCBC[i][1]->Evaluate(),
          CCBC[i][2]->Evaluate());
    }

  out_poly->GetPointData()->SetScalars(arrRad);
  out_poly->GetPointData()->AddArray(arrBC);

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(out_poly);
  writer->SetFileName("circumcenter.vtk");
  writer->Update();
}




/**
 * Compute edge and triangle properties (areas, normals, lengths) using
 * only quadratic expressions and constraints
 */
VarVecArray dummyArray;

void ComputeTriangleAndEdgeProperties(
    ConstrainedNonLinearProblem *p,
    TriangleMesh *mesh,
    const VarVecArray &X,   // Vertices
    VarVecArray &NT,        // Triangle normals
    VarVec &AT,             // Triangle areas
    double minArea,         // Smallest area allowed
    bool doEdges = false,   // Triangle edge lengths (opposite each vertex)
    VarVecArray &TEL = dummyArray)
{
  // Initialize to arrays of NULLS
  std::fill(AT.begin(), AT.end(), (Expression *) NULL);
  std::fill(NT.begin(), NT.end(), VarVec(3, NULL));
  AT.resize(mesh->triangles.size(), NULL);
  NT.resize(mesh->triangles.size(), VarVec(3, NULL));

  if(doEdges)
    {
    std::fill(TEL.begin(), TEL.end(), VarVec(3, NULL));
    TEL.resize(mesh->triangles.size(), VarVec(3, NULL));
    }

  // Iterate over all the triangles in this mesh
  for(int it = 0; it < mesh->triangles.size(); it++)
    {
    // Here is a triangle
    Triangle &t = mesh->triangles[it];
    size_t *v = t.vertices;

    // The un-normalized normal of the triangle
    VarVec Xu = VectorApplyPairwise<BinaryDifference>(p, X[v[1]], X[v[0]]);
    VarVec Xv = VectorApplyPairwise<BinaryDifference>(p, X[v[2]], X[v[0]]);
    VarVec Xu_cross_Xv = CrossProduct(p, Xu, Xv);

    // The area is half the norm of this cross product
    vnl_vector_fixed<double, 3> v_Xu_cross_Xv = VectorEvaluate(Xu_cross_Xv);
    double v_area = 0.5 * v_Xu_cross_Xv.magnitude();
    vnl_vector_fixed<double, 3> v_normal = 0.5 * v_Xu_cross_Xv / v_area;

    // Create the variables for the triangle area and normal
    AT[it] = p->AddVariable("AT", v_area, minArea);
    for(int d = 0; d < 3; d++)
      NT[it][d] = p->AddVariable("NT", v_normal[d]);

    // Create the constraint relating the area and the normal
    for(int d = 0; d < 3; d++)
      {
      Expression *con = new BinaryDifference(
            p,
            new ScalarProduct(p, new BinaryProduct(p, AT[it], NT[it][d]), 2.0),
            Xu_cross_Xv[d]);

      if(fabs(con->Evaluate()) > 1e-6)
        std::cout << "Con_TA-TN: " << con->Evaluate() << std::endl;

      p->AddConstraint(con, "TA-TN", 0, 0);
      }

    // Normal is length one
    Expression *normlen = DotProduct(p, NT[it], NT[it]);
    p->AddConstraint(normlen, "TN.TN", 1.0, 1.0);

    // Compute the edges
    if(doEdges)
      {
      for(int d = 0; d < 3; d++)
        {
        // The edge may have been set up by the opposite triangle
        if(TEL[it][d] != NULL)
          continue;

        // The opposite vertices
        size_t v1 = v[(d + 1) % 3], v2 = v[(d + 2) % 3];

        // Set up the egde length expression
        Expression *edgeLenSq = DistanceSqr(p, X[v1], X[v2]);

        // Create a variable for the edge length
        Expression *edgeLen = p->AddVariable("EL", sqrt(edgeLenSq->Evaluate()), 0);

        // Create the constraint linking the two
        Expression *con = new BinaryDifference(
              p, new Square(p, edgeLen), edgeLenSq);
        p->AddConstraint(con, "EDGELEN", 0, 0);

        // Assign the edge length variable to the current and opposite triangle
        TEL[it][d] = edgeLen;
        if(t.neighbors[d] != NOID)
          {
          TEL[t.neighbors[d]][t.nedges[d]] = edgeLen;
          }
        }
      }
    }
}

#include "BruteForceSubdivisionMedialModel.h"

Expression *ComputeDistanceToCorrespondingMeshObjective(
  ConstrainedNonLinearProblem *p,
  VarVecArray &X,
  std::vector<SMLVec3d> &targetPoint)
{
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

  return objSqDist;
}

Expression *ComputeDistanceToMeshObjective(
    ConstrainedNonLinearProblem *p,
    ClosestPointMatcher *cpm,
    VarVecArray &X)
{
  // For each point on the model, find the closest target point
  std::vector<SMLVec3d> targetPoint = cpm->FindClosestToTarget(X);

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

  return objSqDist;
}

Expression *ComputeDistanceToModelObjective(
    ConstrainedNonLinearProblem *p,
    ClosestPointMatcher *cpm,
    TriangleMesh *mesh,
    VarVecArray &X)
{
  std::vector<ClosestPointMatcher::MatchLocation> meshToModel
      = cpm->FindClosestToSource(mesh, X);

  BigSum *objRecipSqDist = new BigSum(p);
  for(int i = 0; i < meshToModel.size(); i++)
    {
    ClosestPointMatcher::MatchLocation loc = meshToModel[i];
    size_t *v = mesh->triangles[loc.iTriangle].vertices;

    // The closest point based on the current barycentric coordinates
    for(int d = 0; d < 3; d++)
      {
      WeightedSumGenerator wsg(p);
      for(int j = 0; j < 3; j++)
        wsg.AddTerm(X[v[j]][d], loc.xBary[j]);

      wsg.AddConstant(-loc.xTarget[d]);
      objRecipSqDist->AddSummand(new Square(p, wsg.GenerateSum()));
      }

    }

  return objRecipSqDist;
}

void SaveBoundaryMesh(const char *file,
                      ConstrainedNonLinearProblem *p,
                      TriangleMesh *bmesh,
                      std::vector<int> &mIndex,
                      std::vector<std::vector<int> > &mtbIndex,
                      VarVecArray &X,
                      VarVecArray &N,
                      const VarVec &R)
{
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  pts->Allocate(X.size());

  vtkSmartPointer<vtkFloatArray> rad = vtkSmartPointer<vtkFloatArray>::New();
  rad->SetNumberOfComponents(1);
  rad->Allocate(X.size());
  rad->SetName("Radius");

  vtkSmartPointer<vtkIntArray> mix = vtkSmartPointer<vtkIntArray>::New();
  mix->SetNumberOfComponents(1);
  mix->Allocate(X.size());
  mix->SetName("MedialIndex");

  vtkSmartPointer<vtkIntArray> mult = vtkSmartPointer<vtkIntArray>::New();
  mult->SetNumberOfComponents(1);
  mult->Allocate(X.size());
  mult->SetName("Tangency");

  vtkSmartPointer<vtkFloatArray> norm = vtkSmartPointer<vtkFloatArray>::New();
  norm->SetNumberOfComponents(3);
  norm->Allocate(X.size());

  for(int i = 0; i < X.size(); i++)
    {
    int j = mIndex[i];
    pts->InsertNextPoint(X[i][0]->Evaluate(), X[i][1]->Evaluate(), X[i][2]->Evaluate());
    norm->InsertNextTuple3(N[i][0]->Evaluate(), N[i][1]->Evaluate(), N[i][2]->Evaluate());
    rad->InsertNextTuple1(R[j]->Evaluate());
    mix->InsertNextTuple1(j);
    mult->InsertNextTuple1(mtbIndex[j].size());
    }

  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  pd->Allocate(bmesh->triangles.size());
  pd->SetPoints(pts);
  pd->GetPointData()->SetNormals(norm);
  pd->GetPointData()->AddArray(mix);
  pd->GetPointData()->AddArray(mult);
  pd->GetPointData()->AddArray(rad);

  for(int i = 0; i < bmesh->triangles.size(); i++)
    {
    vtkIdType vtx[3];
    for(int j = 0; j < 3; j++)
      vtx[j] = bmesh->triangles[i].vertices[j];
    pd->InsertNextCell(VTK_TRIANGLE, 3, vtx);
    }

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetInputData(pd);
  writer->SetFileName(file);
  writer->Update();
}



void SaveMedialMesh(const char *file,
                    ConstrainedNonLinearProblem *p,
                    TriangleMesh *bmesh,
                    std::vector<int> &mIndex,
                    std::vector<std::vector<int> > &mtbIndex,
                    const VarVecArray &M,
                    const VarVec &R,
                    const VarVecArray &X,
                    const VarVec &TA)
{
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  pts->Allocate(M.size());

  vtkSmartPointer<vtkFloatArray> rad = vtkSmartPointer<vtkFloatArray>::New();
  rad->SetNumberOfComponents(1);
  rad->Allocate(M.size());
  rad->SetName("Radius Function");

  vtkSmartPointer<vtkFloatArray> spoke[3];
 
  spoke[0] = vtkSmartPointer<vtkFloatArray>::New();
  spoke[0]->SetNumberOfComponents(3);
  spoke[0]->Allocate(M.size());
  spoke[0]->SetName("Spoke1");

  spoke[1] = vtkSmartPointer<vtkFloatArray>::New();
  spoke[1]->SetNumberOfComponents(3);
  spoke[1]->Allocate(M.size());
  spoke[1]->SetName("Spoke2");

  spoke[2] = vtkSmartPointer<vtkFloatArray>::New();
  spoke[2]->SetNumberOfComponents(3);
  spoke[2]->Allocate(M.size());
  spoke[2]->SetName("Spoke3");

  for(int i = 0; i < M.size(); i++)
    {
    vnl_vector_fixed<double, 3> xm;
    for(int j = 0; j < 3; j++)
      xm[j] = M[i][j]->Evaluate();

    pts->InsertNextPoint(xm[0], xm[1], xm[2]);
    rad->InsertNextTuple1(R[i]->Evaluate());
    
    std::vector<int> &ix = mtbIndex[i];

    vnl_vector_fixed<double, 3> xs;
    for(int k = 0; k < std::min((int) ix.size(), 3); k++)
      {
      for(int j = 0; j < 3; j++)
        xs[j] = X[ix[k]][j]->Evaluate() - xm[j];
      spoke[k]->InsertNextTuple3(xs[0], xs[1], xs[2]);
      }
    for(int q = ix.size(); q < 3; q++)
      {
      spoke[q]->InsertNextTuple3(xs[0], xs[1], xs[2]);
      }
    }

  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  pd->Allocate(bmesh->triangles.size());
  pd->SetPoints(pts);
  pd->GetPointData()->SetScalars(rad);
  pd->GetPointData()->AddArray(spoke[0]);
  pd->GetPointData()->AddArray(spoke[1]);
  pd->GetPointData()->AddArray(spoke[2]);

  for(int i = 0; i < bmesh->triangles.size(); i++)
    {
    vtkIdType vtx[3];
    for(int j = 0; j < 3; j++)
      vtx[j] = mIndex[bmesh->triangles[i].vertices[j]];
    pd->InsertNextCell(VTK_TRIANGLE, 3, vtx);
    }

  if(TA.size())
    {
    vtkSmartPointer<vtkFloatArray> t_area = vtkSmartPointer<vtkFloatArray>::New();
    t_area->SetNumberOfComponents(1);
    t_area->Allocate(bmesh->triangles.size());
    t_area->SetName("Triangle area");

    for(int i = 0; i < bmesh->triangles.size(); i++)
      t_area->InsertNextTuple1(TA[i]->Evaluate());

    pd->GetCellData()->AddArray(t_area);
    }

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetInputData(pd);
  writer->SetFileName(file);
  writer->Update();
}


/**
 * Data extracted from the input VTK file
 */
class BCMTemplate
{
public:

  SubdivisionSurface::MeshLevel bmesh;
  std::vector<SMLVec3d> x, Nx;
  std::vector<double> R;
  std::vector<int> mIndex;

  void Load(const char *file);
  void Save(const char *file);

  void Subdivide();
  void ImportFromCMRep(const char *file);

};

void BCMTemplate::Load(const char *file)
{
  // Load the input mesh
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(file);
  reader->Update();
  vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();

  // Generate the mesh object
  TriangleMeshGenerator gen(&bmesh, pd->GetNumberOfPoints());
  for(int i = 0; i < pd->GetNumberOfCells(); i++)
    {
    vtkCell *c = pd->GetCell(i);
    if(c->GetNumberOfPoints() != 3)
      throw MedialModelException("Bad cell in input");
    gen.AddTriangle(c->GetPointId(0), c->GetPointId(1), c->GetPointId(2));
    }
  gen.GenerateMesh();

  vtkSmartPointer<vtkDataArray> mix = pd->GetPointData()->GetArray("MedialIndex");
  if(!mix)
    throw MedialModelException("Missing a medial index in input mesh");

  // Get the coordinates and the medial index
  x.resize(pd->GetNumberOfPoints());
  mIndex.resize(pd->GetNumberOfPoints());
  for(int i = 0; i < pd->GetNumberOfPoints(); i++)
    {
    x[i].set(pd->GetPoint(i));
    mIndex[i] = floor(0.5 + mix->GetTuple1(i));
    }
}

void BCMTemplate::Save(const char *file)
{
  // Save as a VTK mesh
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  pts->Allocate(bmesh.nVertices);

  vtkSmartPointer<vtkIntArray> mix = vtkSmartPointer<vtkIntArray>::New();
  mix->SetNumberOfComponents(1);
  mix->Allocate(bmesh.nVertices);
  mix->SetName("MedialIndex");

  vtkSmartPointer<vtkFloatArray> nrm = vtkSmartPointer<vtkFloatArray>::New();
  nrm->SetNumberOfComponents(3);
  nrm->Allocate(bmesh.nVertices);
  
  vtkSmartPointer<vtkFloatArray> rad = vtkSmartPointer<vtkFloatArray>::New();
  rad->SetNumberOfComponents(1);
  rad->Allocate(bmesh.nVertices);
  rad->SetName("Radius");

  for(int i = 0; i < bmesh.nVertices; i++)
    {
    pts->InsertNextPoint(x[i].data_block());
    nrm->InsertNextTuple(Nx[i].data_block());
    rad->InsertNextTuple1(R[i]);
    mix->InsertNextTuple1(mIndex[i]);
    }

  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  pd->Allocate(bmesh.triangles.size());
  pd->SetPoints(pts);
  pd->GetPointData()->AddArray(mix);
  pd->GetPointData()->AddArray(rad);
  pd->GetPointData()->SetNormals(nrm);

  for(int i = 0; i < bmesh.triangles.size(); i++)
    {
    vtkIdType vtx[3];
    for(int j = 0; j < 3; j++)
      vtx[j] = bmesh.triangles[i].vertices[j];
    pd->InsertNextCell(VTK_TRIANGLE, 3, vtx);
    }

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetInputData(pd);
  writer->SetFileName(file);
  writer->Update();
}

struct Edge : public std::pair<int, int>
{
  Edge(int a, int b) :
    std::pair<int, int>()
  {
    this->first = std::min(a, b);
    this->second = std::max(a, b);
  }
};

void BCMTemplate::Subdivide()
{
  // Make a copy of the boundary mesh
  bmesh.SetAsRoot();
  SubdivisionSurface::MeshLevel src = bmesh;

  // Subdivide into the current mesh
  SubdivisionSurface::Subdivide(&src, &bmesh);

  // Get the subdivision matrix
  typedef ImmutableSparseMatrix<double> SMat;
  SMat &W = bmesh.weights;

  // Find all unique edges
  int iNext = *std::max_element(mIndex.begin(), mIndex.end()) + 1;
  typedef std::map<Edge, int> EdgeIndex;
  EdgeIndex edgeIndex;

  // Assign a new mIndex to every unique edge in every triangle
  for(int t = 0; t < src.triangles.size(); t++)
    {
    size_t *v = src.triangles[t].vertices;
    for(int j = 0; j < 3; j++)
      {
      Edge e(mIndex[v[(j+1) % 3]], mIndex[v[(j+2) % 3]]);
      std::pair<EdgeIndex::iterator, bool> ret = edgeIndex.insert(std::make_pair(e, 0));
      if(ret.second)
        {
        ret.first->second = iNext++;
        std::cout << "Edge " << e.first << "," << e.second << " to " << iNext << std::endl;
        }
      else
        std::cout << "Already found " << e.first << "," << e.second << std::endl;
      }
    }

  // Apply the subdivision to the mesh coordinates
  std::vector<SMLVec3d> xsrc = x;
  std::vector<SMLVec3d> Nsrc = Nx;
  std::vector<double> Rsrc = R;
  x.clear(); Nx.clear(); Rsrc.clear();

  for(int i = 0; i < bmesh.nVertices; i++)
    {
    SMLVec3d xi(0.0, 0.0, 0.0), Ni(0.0, 0.0, 0.0);
    double Ri = 0.0;
    for(SMat::RowIterator it = W.Row(i); !it.IsAtEnd(); ++it)
      {
      xi += xsrc[it.Column()] * it.Value();
      Ni += Nsrc[it.Column()] * it.Value();
      Ri += R[it.Column()] * it.Value();
      }
    x.push_back(xi); Nx.push_back(Ni); R.push_back(Ri);
    }

  // Compute the m-indices
  mIndex.resize(bmesh.nVertices, -1);
  for(int t = 0; t < src.triangles.size(); t++)
    {
    // The new vertices inserted into this triangle
    size_t *vnew = bmesh.triangles[4 * t + 3].vertices;
    size_t *vold = src.triangles[t].vertices;

    // Find the edge for each triangle
    for(int j = 0; j < 3; j++)
      {
      if(mIndex[vnew[j]] == -1)
        {
        // Find the index assigned to this edge
        Edge e(mIndex[vold[(j+1) % 3]], mIndex[vold[(j+2) % 3]]);
        mIndex[vnew[j]] = edgeIndex[e];
        std::cout << "Vertex " << vnew[j] << " MIndex " << mIndex[vnew[j]] << std::endl;
        }
      }
    }

  // Set the current level as root
  bmesh.SetAsRoot();
}

int MatchToOppositeTriangle(int index, int level)
{
  // Build a list of index modulo 4 at each level
  std::vector<int> offset;
  int q = index;
  for(int i = 0; i < level; i++)
    {
    offset.push_back(q % 4);
    q = q / 4;
    }

  // Swap offsets 0 and 2
  for(int i = 0; i < level; i++)
    {
    if(offset[i] == 0)
      offset[i] = 2;
    else if(offset[i] == 2)
      offset[i] = 0;
    }

  // Rebuild the index
  for(int i = level-1; i >= 0; i--)
    {
    q = q * 4 + offset[i];
    }

  return q;
}

void BCMTemplate::ImportFromCMRep(const char *file)
{
  // Load and process the mrep (TODO: remove)
  MedialPDE mrep(file);
  SubdivisionMedialModel *tmpmodel =
      dynamic_cast<SubdivisionMedialModel *>(mrep.GetMedialModel());

  // Get the coefficient medial mesh - we will subdivide this mesh
  const SubdivisionSurface::MeshLevel *mc = tmpmodel->GetCoefficientMesh();
  const SubdivisionSurface::MeshLevel *ma = tmpmodel->GetAtomMesh();

  // Subdivision level is log_4(ma/mc)
  int tratio = ma->triangles.size() / mc->triangles.size();
  int level = 0;
  while(tratio > (1 << (level * 2)))
    level++;

  // Create a coarse-level boundary mesh
  SubdivisionSurface::MeshLevel *bc = new SubdivisionSurface::MeshLevel();

  // Create a mapping from atom vertices to boundary vertices that is monotonic
  typedef vnl_vector_fixed<int, 2> IntPair;
  std::vector<IntPair> mtobVertexMap;
  int boundaryVertexIndex = 0;
  for(int i = 0; i < mc->nVertices; i++)
    {
    EdgeWalkAroundVertex walk(mc, i);
    IntPair b;
    b[0] = boundaryVertexIndex++;
    b[1] = walk.IsOpen() ? b[0] : boundaryVertexIndex++;
    mtobVertexMap.push_back(b);
    }

  // We will use a generator to create the boundary mesh
  TriangleMeshGenerator gen(bc, boundaryVertexIndex);

  // Add boundary triangles in two passes
  for(int d = 0; d < 2; d++)
    {
    for(int t = 0; t < mc->triangles.size(); t++)
      {
      const size_t *v = mc->triangles[t].vertices;
      int b1 = mtobVertexMap[v[0]][d];
      int b2 = mtobVertexMap[v[1]][d];
      int b3 = mtobVertexMap[v[2]][d];
      if(d == 0)
        gen.AddTriangle(b3, b2, b1);
      else
        gen.AddTriangle(b1, b2, b3);
      }
    }

  // Generate the coarse boundary mesh
  gen.GenerateMesh();

  // Subdivide the boundary mesh
  SubdivisionSurface::MeshLevel *ba = new SubdivisionSurface::MeshLevel();
  SubdivisionSurface::RecursiveSubdivide(bc, ba, level);

  // At this point, we have the boundary mesh of proper topology, but we don't
  // yet know which medial atoms the boundary atoms link to! We need to find a
  // medial atom and side for each boundary atom
  x.resize(ba->nVertices);
  Nx.resize(ba->nVertices);
  R.resize(ba->nVertices);
  mIndex.resize(ba->nVertices, -1);

  // Visit all vertices in all triangles
  for(int t = 0; t < ba->triangles.size(); t++)
    {
    // What side of the medial axis we're on
    int side = t / (ba->triangles.size() / 2);

    // Get the corresponding medial triangle
    int mt = (side == 0)
        ? MatchToOppositeTriangle(t, level)
        : t - (ba->triangles.size() / 2);

    // Loop over all three vertices
    for(int j = 0; j < 3; j++)
      {
      // Get the corresponding medial vertex
      int mj = (side == 1) ? j : 2 - j;

      // Find the information from the medial data
      int vb = ba->triangles[t].vertices[j];
      if(mIndex[vb] < 0)
        {
        mIndex[vb] = ma->triangles[mt].vertices[mj];
        x[vb] = tmpmodel->GetAtomArray()[mIndex[vb]].xBnd[side].X;

        // Compute normal and radius
        Nx[vb] = x[vb] - tmpmodel->GetAtomArray()[mIndex[vb]].X;
        R[vb] = Nx[vb].magnitude();
        Nx[vb] /= R[vb];
        }
      }
    }

  // Store as the mesh
  this->bmesh = *ba;

  // Delete intermediates
  delete bc; delete ba;
}


int ConvertCMRepToBoundaryRepresentation(std::string fnInput, std::string fnOutput)
{
  // We want to load a triangular mesh in which each point is associated with
  // a tag (corresponding medial atom).
  BCMTemplate tmpl;
  tmpl.ImportFromCMRep(fnInput.c_str());
  tmpl.Save(fnOutput.c_str());
  return 0;
}

void FixCmrepMedialMesh(std::string fnInput, std::string fnOutput)
{
  // Load the m-rep
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(fnInput.c_str());
  reader->Update();

  // Convert it into a triangle mesh
  vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();
  SubdivisionSurface::MeshLevel mesh;
  SubdivisionSurface::ImportLevelFromVTK(pd, mesh);

  vtkSmartPointer<vtkDataArray> ar = pd->GetPointData()->GetArray("Radius");

  // Get the points and the radii
  std::vector<SMLVec3d> x;
  std::vector<double> r;
  for(int i = 0; i < mesh.nVertices; i++)
    {
    x.push_back(SMLVec3d(pd->GetPoint(i)));
    r.push_back(ar ? ar->GetTuple1(i) : 0.5);
    }

  // Search for triangles to split
  std::set<int> splits;
  int iLast = mesh.nVertices;

  SubdivisionSurface::MeshLevel mfix;

  // Create a mesh generator
  TriangleMeshGenerator gen(&mfix, 0);

  for(int i = 0; i < mesh.triangles.size(); i++)
    {
    // Triangles already split should be ignored
    if(splits.find(i) != splits.end())
      continue;

    bool isBad = false;
    size_t *v = mesh.triangles[i].vertices;
    for(int j = 0; j < 3; j++)
      {
      size_t v1 = v[(j+1) % 3];
      size_t v2 = v[(j+2) % 3];
      EdgeWalkAroundVertex w1(&mesh, v1);
      EdgeWalkAroundVertex w2(&mesh, v2);
      size_t nj = mesh.triangles[i].neighbors[j];
      if(w1.IsOpen() && w2.IsOpen() && nj != NOID)
        {
        // This triangle is bad - we need to split it
        std::cout << "Found bad edge " << v1 << ", " << v2 << std::endl;
        isBad = true;

        // Add the four new triangles
        int vopp = mesh.triangles[nj].vertices[mesh.triangles[i].nedges[j]];
        gen.AddTriangle(v[j], v1, iLast);
        gen.AddTriangle(v[j], iLast, v2);
        gen.AddTriangle(v1, vopp, iLast);
        gen.AddTriangle(iLast, vopp, v2);

        // Create a new point
        x.push_back((x[v1] + x[v2]) * 0.5);
        r.push_back((r[v[j]] + r[vopp]) * 0.5);
        ++iLast;

        // Mark the opposite triangle to be ignored
        splits.insert(nj);
        }
      }

    if(!isBad)
      {
      gen.AddTriangle(v[0], v[1], v[2]);
      }
    }

  mfix.nVertices = iLast;
  gen.GenerateMesh();



  // Save the new mesh
  vtkSmartPointer<vtkFloatArray> ar_new = vtkSmartPointer<vtkFloatArray>::New();
  ar_new->SetNumberOfComponents(1);
  ar_new->Allocate(mfix.nVertices);
  ar_new->SetName("Radius");

  vtkSmartPointer<vtkPoints> p_new = vtkSmartPointer<vtkPoints>::New();
  p_new->Allocate(mfix.nVertices);


  for(int i = 0; i < mfix.nVertices; i++)
    {
    p_new->InsertNextPoint(x[i].data_block());
    ar_new->InsertNextTuple1(r[i]);
    }

  vtkSmartPointer<vtkPolyData> pd_new = vtkSmartPointer<vtkPolyData>::New();
  pd_new->SetPoints(p_new);
  pd_new->GetPointData()->AddArray(ar_new);

  SubdivisionSurface::ExportLevelToVTK(mfix, pd_new);

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(pd_new);
  writer->SetFileName(fnOutput.c_str());
  writer->Update();

}

int SubdivideBoundaryRepresentation(
    std::string fnTemplate,
    int subdivisionLevel,
    std::string fnOutput)
{
  // Load the template!
  BCMTemplate tmpl;
  tmpl.Load(fnTemplate.c_str());

  // Create the output template
  for(int i = 0; i < subdivisionLevel; i++)
    tmpl.Subdivide();

  // Save the mesh
  tmpl.Save(fnOutput.c_str());

  return 0;
}

#include <vnl/algo/vnl_svd.h>

/**
 * Un-subdivide a mesh by one level. This only computes the triangles and
 * nVertices, nothing else.
 */
void UnSubdivide(
    SubdivisionSurface::MeshLevel &mesh,
    SubdivisionSurface::MeshLevel &parent)
{
  // Initialize the parent mesh
  parent.triangles.clear();
  parent.nVertices = 0;

  // Create each of the parent triangles
  for(int i = 0; i < mesh.triangles.size(); i+=4)
    {
    // Create the coarse triangle
    Triangle tc;

    // Repeat for each vertex
    for(int j = 0; j < 3; j++)
      {
      // Configure the vertices using the subdivided triangles
      tc.vertices[j] = mesh.triangles[i + j].vertices[j];

      // Update the max vertex
      parent.nVertices = std::max(tc.vertices[j]+1, parent.nVertices);

      // Compute the neighbors
      tc.neighbors[j] = mesh.triangles[i + ((j+1)%3)].neighbors[j] / 4;

      // Compute the neighbor edge indexes
      tc.nedges[j] = mesh.triangles[i + ((j+1)%3)].nedges[j];
      }

    // Add the triangle
    parent.triangles.push_back(tc);
    }
}

/**
 * Reverse engineer a boundary mesh. For meshes that were generated by
 * subdivision, this method will create a coarse-level mesh and fit its
 * vertex values to the data in the current mesh
 */
bool ReverseEngineerSubdivisionMesh(
    BCMTemplate &tChild, BCMTemplate &tParent, int level)
{
  // The number of triangles must be divisible by four^level
  if(tChild.bmesh.triangles.size() % (1 << (2 * level)))
    return false;

  // For each level, produce a parent mesh
  SubdivisionSurface::MeshLevel mlCurrent = tChild.bmesh;
  SubdivisionSurface::MeshLevel mlParent;

  for(int i = 0; i < level; i++)
    {
    UnSubdivide(mlCurrent, mlParent);
    if(!SubdivisionSurface::CheckMeshLevel(&mlParent))
      return false;
    mlCurrent = mlParent;
    }

  // Finish initializing the parent
  tParent.bmesh = mlParent;
  tParent.bmesh.SetAsRoot();
  tParent.bmesh.ComputeWalks();

  // Check that the parent subdivides back into the child
  SubdivisionSurface::MeshLevel cmp;
  SubdivisionSurface::RecursiveSubdivide(&tParent.bmesh, &cmp, level);

  if(tChild.bmesh.nVertices != cmp.nVertices)
    {
    printf("Vertex count mismatch %d vs %d\n", (int) tChild.bmesh.nVertices, (int) cmp.nVertices);
    return false;
    }

  for(int i = 0; i < cmp.triangles.size(); i++)
    {
    Triangle t1 = cmp.triangles[i];
    Triangle t2 = tChild.bmesh.triangles[i];
    for(int j = 0; j < 3; j++)
      {
      if(t1.vertices[j] != t2.vertices[j])
        {
        printf("vertex mismatch %d:%d\n", i, j);
        return false;
        }
      if(t1.neighbors[j] != t2.neighbors[j])
        {
        printf("neighbors mismatch %d:%d\n", i, j);
        return false;
        }
      if(t1.nedges[j] != t2.nedges[j])
        {
        printf("nedges mismatch %d:%d\n", i, j);
        return false;
        }
      }
    }

  // At this point, we just want to take the weight matrix from cmp and
  // put it into mesh
  tChild.bmesh.parent = cmp.parent;
  tChild.bmesh.weights = cmp.weights;

  // Find the coordinates for the parent mesh that best fit the coordinates
  // in the child mesh. This involves solving the problem
  vnl_matrix<double> W = tChild.bmesh.weights.GetDenseMatrix();

  // TODO: lame!
  vnl_svd<double> svd(W.transpose() * W);

  // Compute the x, y, z coordinates
  tParent.x.resize(tParent.bmesh.nVertices);
  for(int j = 0; j < 3; j++)
    {
    vnl_vector<double> xi(tChild.x.size());
    for(int k = 0; k < xi.size(); k++)
      xi[k] = tChild.x[k][j];

    vnl_vector<double> bi = tChild.bmesh.weights.MultiplyTransposeByVector(xi);

    vnl_vector<double> yi = svd.solve(bi);
    for(int k = 0; k < yi.size(); k++)
      tParent.x[k][j] = yi[k];
    }

  // Just copy the mindex
  tParent.mIndex.resize(tParent.bmesh.nVertices);
  for(int i = 0; i < tParent.bmesh.nVertices; i++)
    tParent.mIndex[i] = tChild.mIndex[i];

  return true;
}

/**
 * A function to save the current model as well as different gradients,
 * and how a step in the gradient direction affects certain constraints
 */
void SaveConstraint()
{


}

enum ProgramAction {
  ACTION_FIT_TARGET,
  ACTION_FIT_SELF,
  ACTION_SUBDIVIDE,
  ACTION_CONVERT_CMREP,
  ACTION_FIX_CMREP,
  ACTION_NONE
};

int usage()
{
  const char *usage =
      "bcmrep - Boundary-constrained medial representation utility\n"
      "usage:\n"
      "  bcmrep [options]\n"
      "options that select an action (use one):\n"
      "  -icp template.vtk target.vtk   : fit template mesh to target mesh\n"
      "  -lsq template.vtk target.vtk   : least squares fit to target mesh\n"
      "  -sub template.vtk factor       : subdivide template by a factor\n"
      "  -cmr input.cmrep               : import a cm-rep template\n"
      "  -cfx input.cmrep               : fix a cm-rep template with bad triangles\n"
      "other required options:\n"
      "  -o output.vtk                  : mesh to save the result\n"
      "optional parameters:\n"
      "  -reg-lb NxNxN                  : use Laplace basis regularization \n"
      "                                           with N basis functions at each \n"
      "                                           iteration. (default, with N=20)\n"
      "  -reg-ss L                      : use subdivision surface-based \n"
      "                                           regularlization with L levels of \n"
      "                                           subdivision. \n"
      "  -reg-weight NxNxN              : weight of the regularization term\n"
      "                                           at each iteration\n"
      "  -reg-el W                      : penalize the length of the crest\n"
      "                                           curve with weight specified\n"
      "triangle shape constraints:\n"
      "  -bc alpha,beta,min_area        : Set boundary triangle constraints \n"
      "                                     alpha:    minimum triangle angle (def: 12)\n"
      "                                     beta:     minimum dihedral angle (def: 120)\n"
      "                                     min_area: minimum triangle area  (def: 0.1)\n"
      "  -mc alpha,beta,min_area        : Set medial triangle constraints (see above) \n"
      "                                     by default all three are set to 0 (no constraint)\n"
      "IPOpt options:\n"
      "  -solver NAME                   : select which solver to use (see IPOpt options. Def: ma86)\n"
      "  -no-hess                       : Turn off analytical hessian (does not work well)\n";
  std::cout << usage;
  return -1;
}

struct RegularizationOptions
{
  enum RegularizationMode { SUBDIVISION, SPECTRAL };

  int SubdivisionLevel, BasisSize;
  double Weight;
  RegularizationMode Mode;
  std::string Solver;
  bool UseHessian;

  double EdgeLengthWeight;

  RegularizationOptions()
    : SubdivisionLevel(0), BasisSize(20),
      Weight(4.0), Mode(SPECTRAL), Solver("ma86"),
      EdgeLengthWeight(0.0), UseHessian(true) {}
};

class ObjectiveCombiner
{
public:
  struct Entry {
    Expression *exp;
    double weight;
    double initValue;
  };

  void AddObjectiveTerm(Expression *obj, std::string name, double weight)
  {
    Entry ent;
    ent.exp = obj;
    ent.weight = weight;
    ent.initValue = obj->Evaluate();
    m_Generator.AddTerm(obj, weight);
    m_ObjectiveInfo[name] = ent;
  }

  Expression *GetTotalObjective()
  {
    if(m_TotalObjective == NULL)
      {
      m_TotalObjective = m_Generator.GenerateSum();
      m_InitTotalValue = m_TotalObjective->Evaluate();
      }
    return m_TotalObjective;
  }

  void PrintReport()
  {
    printf("%30s   %10s %10s   %10s   %10s\n", "TERM", "RawValue", "Weight", "WgtValue", "Delta");
    for(InfoMap::iterator it = m_ObjectiveInfo.begin(); it != m_ObjectiveInfo.end(); ++it)
      {
      double val = it->second.exp->Evaluate();
      printf("%30s   %10.4f %10.4f   %10.4f   %10.4f\n",
             it->first.c_str(),
             val, it->second.weight,
             val * it->second.weight,
             val * it->second.weight - it->second.initValue * it->second.weight);
      }
    double totVal = m_TotalObjective->Evaluate();
    printf("%30s   %10.4s %10.4s   %10.4f   %10.4f\n",
           "TOTAL OBJECTIVE", "", "", totVal, totVal - m_InitTotalValue);
  }

  ObjectiveCombiner(ConstrainedNonLinearProblem *p)
    : m_Generator(p)
  {
    m_Problem = p;
    m_TotalObjective = NULL;
  }

protected:
  WeightedSumGenerator m_Generator;
  ConstrainedNonLinearProblem *m_Problem;
  Expression *m_TotalObjective;
  double m_InitTotalValue;
  typedef std::map<std::string, Entry> InfoMap;
  InfoMap m_ObjectiveInfo;
};

struct TriangleConstraint
{
  // Minimum area of triangles
  double min_tri_area;

  // Minimum dihedral angle (120 is a good value)
  double min_dih_angle;

  // Minimum triangle angle (12 is a good value)
  double min_tri_angle;

  TriangleConstraint(double area, double tri_ang, double dih_ang)
    : min_tri_area(area), min_tri_angle(tri_ang), min_dih_angle(dih_ang) {}
};

int main(int argc, char *argv[])
{
  // Usage help
  if(argc < 2) return usage();

  // The action specified for the program
  ProgramAction action = ACTION_NONE;
  std::string fnTemplate, fnTarget, fnOutput, fnImportSource;
  int subdivisionLevel = 0;

  RegularizationOptions regOpts;

  // Triangle constraints for medial and boundary triangles
  TriangleConstraint tc_med(0,0,0);
  TriangleConstraint tc_bnd(0.1, 12, 120);

  for(int p = 1; p < argc; p++)
    {
    std::string cmd = argv[p];
    if(cmd == "-icp")
      {
      action = ACTION_FIT_TARGET;
      fnTemplate = argv[++p];
      fnTarget = argv[++p];
      }
    else if(cmd == "-lsq")
      {
      action = ACTION_FIT_SELF;
      fnTemplate = argv[++p];
      fnTarget = argv[++p];
      }
    else if(cmd == "-sub")
      {
      action = ACTION_SUBDIVIDE;
      fnTemplate = argv[++p];
      subdivisionLevel = atoi(argv[++p]);
      }
    else if(cmd == "-cmr")
      {
      action = ACTION_CONVERT_CMREP;
      fnImportSource = argv[++p];
      }
    else if(cmd == "-cfx")
      {
      action = ACTION_FIX_CMREP;
      fnImportSource = argv[++p];
      }
    else if(cmd == "-o")
      {
      fnOutput = argv[++p];
      }
    else if(cmd == "-o")
      {
      return usage();
      }
    else if(cmd == "-reg-lb")
      {
      regOpts.Mode = RegularizationOptions::SPECTRAL;
      regOpts.BasisSize = atoi(argv[++p]);
      }
    else if(cmd == "-reg-ss")
      {
      regOpts.Mode = RegularizationOptions::SUBDIVISION;
      regOpts.SubdivisionLevel = atoi(argv[++p]);
      }
    else if(cmd == "-reg-weight")
      {
      regOpts.Weight = atof(argv[++p]);
      }
    else if(cmd == "-reg-el")
      {
      regOpts.EdgeLengthWeight = atof(argv[++p]);
      }
    else if(cmd == "-solver")
      {
      regOpts.Solver = argv[++p];
      }
    else if(cmd == "-no-hess")
      {
      regOpts.UseHessian = false;
      }
    else if(cmd == "-bc")
      {
      tc_bnd.min_tri_angle = atof(argv[++p]);
      tc_bnd.min_dih_angle = atof(argv[++p]);
      tc_bnd.min_tri_area = atof(argv[++p]);
      }
    else if(cmd == "-mc")
      {
      tc_med.min_tri_angle = atof(argv[++p]);
      tc_med.min_dih_angle = atof(argv[++p]);
      tc_med.min_tri_area = atof(argv[++p]);
      }
    else
      {
      std::cerr << "Unknown command " << cmd << std::endl;
      return -1;
      }
    }

  // Decide what to do based on the action!
  if(action == ACTION_CONVERT_CMREP)
    {
    return ConvertCMRepToBoundaryRepresentation(fnImportSource, fnOutput);
    }

  else if(action == ACTION_FIX_CMREP)
    {
    FixCmrepMedialMesh(fnImportSource, fnOutput);
    return 0;
    }

  else if(action == ACTION_SUBDIVIDE)
    {
    return SubdivideBoundaryRepresentation(fnTemplate, subdivisionLevel, fnOutput);
    }

  else if(action == ACTION_NONE)
    {
    std::cerr << "No action specified" << std::endl;
    return -1;
    }

  // Load the template
  BCMTemplate tmpl;
  tmpl.Load(fnTemplate.c_str());

  // We will restrict our operations to a boundary triangle mesh and a set of medial
  // atom indices for each point on the triangle mesh
  TriangleMesh *bmesh = &tmpl.bmesh;

  // Get the number of boundary points
  int nb = bmesh->nVertices;

  // Initialize the data we are extracting from the boundary mesh
  std::vector<int> &mIndex = tmpl.mIndex;
  std::vector<SMLVec3d> &xInput = tmpl.x;

  // Load the target vertex coordintates if requesting LSQ fit
  std::vector<SMLVec3d> xLSQTarget = xInput;
  if(action == ACTION_FIT_SELF)
    {
    BCMTemplate lsq_target;
    lsq_target.Load(fnTarget.c_str());
    xLSQTarget = lsq_target.x;
    }

  // At this point, we are only working with bMesh and mIndex and xInput

  // Get the number of medial points. TODO: we need to make sure that every
  // index between 0 and nm is represented.
  int nm = 1 + *std::max_element(mIndex.begin(), mIndex.end());

  // Create a list of boundary atoms for each medial atom
  typedef std::vector<std::vector<int> > MedialToBoundaryIndex;
  MedialToBoundaryIndex mtbIndex(nm);
  for(int i = 0; i < nb; i++)
    {
    mtbIndex[mIndex[i]].push_back(i);
    }


  // Load the target mesh
  vtkPolyData *target = ReadVTKMesh(fnTarget.c_str());

  // Load the target image
  /*
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
    */

  // Create the optimization problem
  ConstrainedNonLinearProblem *p = new ConstrainedNonLinearProblem();

  // The boundary positions
  VarVecArray X(nb, VarVec(3, NULL));

  // The medial positions
  VarVecArray M(nm, VarVec(3, NULL));

  // The radius values
  VarVec R(nm, NULL);

  // The boundary normal vectors
  VarVecArray N(nb, VarVec(3, NULL));

  // The vectors U from the medial axis to the boundary. They are the same
  // as -R*N. They appear in multiple places though, so we should store them
  VarVecArray U(nb, VarVec(3, NULL));

  // The triangle areas, normals and edge lengths on the boundary surface
  // and on the medial surface (to avoid square roots)
  VarVec taX, taM;
  VarVecArray NT_X, TEL_X, NT_M, TEL_M;

  // A buffer for making variable names
  char buffer[64];

  // ------------------------------------------------------------------------
  // Configure the boundary point variables
  // ------------------------------------------------------------------------

  // Configure the medial point variables
  for(int i = 0; i < nb; i++)
    {
    SMLVec3d x = xInput[i];

    // Set up the medial coordinate and normal
    for(int j = 0; j < 3; j++)
      {
      sprintf(buffer, "X[%d,%d]", i, j);
      X[i][j] = p->AddVariable(buffer, x[j]);
      }
    }

  // ------------------------------------------------------------------------
  // Configure the constraints on the boundary normal
  // ------------------------------------------------------------------------

  // Create a LoopScheme for specifying normal vector constraints
  LoopTangentScheme lts;
  lts.SetMesh(bmesh);

  // Expressions describing Xu, Xv, Nu, Nv at each point on the boundary
  VarVecArray Xd[] = { VarVecArray(nb, VarVec(3, NULL)), VarVecArray(nb, VarVec(3, NULL)) };
  VarVecArray Nd[] = { VarVecArray(nb, VarVec(3, NULL)), VarVecArray(nb, VarVec(3, NULL)) };

  // Expression for the first fundamental and the shape operator
  VarVec curv_mean(nb, NULL);
  VarVec curv_gauss(nb, NULL);
  VarVec curv_k1(nb, NULL);
  VarVec curv_k2(nb, NULL);

  // Add the constraints relating each normal to the neighboring vertices
  for(int i = 0; i < nb; i++)
    {
    // Repeat for u and v directions
    for(int d = 0; d < 2; d++)
      {
      BigSum *Xdi[] = { new BigSum(p), new BigSum(p), new BigSum(p) };

      double wi = lts.GetOwnWeight(d, i);
      for(int j = 0; j < 3; j++)
        {
        Xdi[j]->AddSummand(new ScalarProduct(p, X[i][j], wi));
        }

      for(EdgeWalkAroundVertex walk(bmesh, i); !walk.IsAtEnd(); ++walk)
        {
        double wij = lts.GetNeighborWeight(d, walk);
        for(int j = 0; j < 3; j++)
          {
          Xdi[j]->AddSummand(new ScalarProduct(p, X[walk.MovingVertexId()][j], wij));
          }
        }

      // Store these expressions
      for(int j = 0; j < 3; j++)
        {
        Xd[d][i][j] = Xdi[j];
        }
      }

    // Compute the initial value of the normal
    SMLVec3d v_xu = VectorEvaluate(Xd[0][i]);
    SMLVec3d v_xv = VectorEvaluate(Xd[1][i]);
    SMLVec3d v_n = vnl_cross_3d(v_xu, v_xv).normalize();

    // Create the variables for the normal
    for(int j = 0; j < 3; j++)
      {
      sprintf(buffer, "N[%d,%d]", i, j);
      N[i][j] = p->AddVariable(buffer, v_n[j]);
      }

    // Create the constraints on the normal
    for(int d = 0; d < 2; d++)
      {
      // Create the constrait Xu . N = 0
      Expression *constrNormXu = DotProduct(p, Xd[d][i], N[i]);

      // Add the constraint to the problem
      p->AddConstraint(constrNormXu, "N.Xu", 0, 0);
      }

    // Add the constraint on the normal being norm 1
    Expression *constrNormMag = MagnitudeSqr(p, N[i]);
    p->AddConstraint(constrNormMag, "N.N", 1.0, 1.0);
    }


  // ------------------------------------------------------------------------
  // Configure the constraints relating boundary curvature to the radius fn
  // ------------------------------------------------------------------------

  for(int i = 0; i < nb; i++)
    {

    // For non-edge atoms, we are done, no need to compute the rest
    if(mtbIndex[mIndex[i]].size() > 1)
      {
      curv_mean[i] = new Constant(p, 0);
      curv_gauss[i] = new Constant(p, 0);
      curv_k1[i] = new Constant(p, 0);
      curv_k2[i] = new Constant(p, 0);
      continue;
      }

    // Compute the expressions Nu, Nv
    // Repeat for u and v directions
    for(int d = 0; d < 2; d++)
      {
      BigSum *Ndi[] = { new BigSum(p), new BigSum(p), new BigSum(p) };

      double wi = lts.GetOwnWeight(d, i);
      for(int j = 0; j < 3; j++)
        {
        Ndi[j]->AddSummand(new ScalarProduct(p, N[i][j], wi));
        }

      for(EdgeWalkAroundVertex walk(bmesh, i); !walk.IsAtEnd(); ++walk)
        {
        double wij = lts.GetNeighborWeight(d, walk);
        for(int j = 0; j < 3; j++)
          {
          Ndi[j]->AddSummand(new ScalarProduct(p, N[walk.MovingVertexId()][j], wij));
          }
        }

      // Store these expressions
      for(int j = 0; j < 3; j++)
        {
        Nd[d][i][j] = Ndi[j];
        }
      }

    // Compute the first fundamental form at the point.
    vnl_matrix_fixed<double, 2, 2> mFF1, mFF2, mSO;
    Expression *FF1[2][2], *FF2_neg[2][2], *SO[2][2];
    for(int q = 0; q < 2; q++)
      {
      for(int r = 0; r < 2; r++)
        {
        // Add the expression for the first f. f. as a constrained variable
        FF1[q][r] = p->AddExpressionAsConstrainedVariable(DotProduct(p, Xd[q][i], Xd[r][i]), "FF1");
        mFF1[q][r] = FF1[q][r]->Evaluate();

        // Minus the second f.f.
        FF2_neg[q][r] = DotProduct(p, Xd[q][i], Nd[r][i]);
        mFF2[q][r] = FF2_neg[q][r]->Evaluate();
        }
      }

    // Numerically solve for the shape operator
    mSO = - vnl_inverse(mFF1) * mFF2;

    // Add the expression for the shape operator
    for(int q = 0; q < 2; q++)
      {
      for(int r = 0; r < 2; r++)
        {
        SO[q][r] = p->AddVariable("SO", mSO[q][r]);
        }
      }

    // Solve for the shape operator
    for(int q = 0; q < 2; q++)
      {
      for(int r = 0; r < 2; r++)
        {
        // Minus the second f.f.
        Expression *sff = DotProduct(p, Xd[q][i], Nd[r][i]);

        // Expression for FFF * S - SFF
        Expression *con = new TernarySum(
              p,
              new BinaryProduct(p, FF1[q][0], SO[0][r]),
              new BinaryProduct(p, FF1[q][1], SO[1][r]),
              sff);

        if(fabs(con->Evaluate()) > 1e-6)
          std::cout << "Con_SO = " << con->Evaluate() << std::endl;

        // Set as a hard constraint
        p->AddConstraint(con, "SO", 0, 0);
        }
      }

    // Numerically solve for k1
    double mH = vnl_trace(mSO) / 2;
    double mK = vnl_det(mSO);
    double mk1 = mH - sqrt(mH*mH - mK);

    // Solve the characteristic polynomial for kappa1
    Expression *k1 = p->AddVariable("k1", mk1);

    // The constraint on kappa1
    Expression * con = new BinaryDifference(p,
                                            new BinaryProduct(p,
                                                              new BinaryDifference(p, SO[0][0], k1),
                                                              new BinaryDifference(p, SO[1][1], k1)),
                                            new BinaryProduct(p, SO[0][1], SO[1][0]));

    // Evaluate the constraint
    p->AddConstraint(con, "Kappa-eq", 0, 0);

    if(fabs(con->Evaluate()) > 1e-6)
      std::cout << "Con_K1 = " << con->Evaluate() << std::endl;

    // We want kappa1 to be the larger in magnitude of the two curvatures. This is equivalent
    // to having kappa1 > trace(SO)/2. K1 is negative.
    Expression *H = new ScalarProduct(p, new BinarySum(p, SO[0][0], SO[1][1]), 0.5);
    Expression *con2 = new BinaryDifference(p, k1, H);

    p->AddConstraint(con2, "Kappa-ineq", ConstrainedNonLinearProblem::LBINF, 0);

    if(con2->Evaluate() > -1e-6)
      std::cout << "Con_K1_sign = " << con2->Evaluate() << std::endl;

    // Store these curvatures for future reference
    curv_mean[i] = H;
    curv_k1[i] = k1;
    curv_k2[i] = new BinaryDifference(p, new ScalarProduct(p, H, 2), k1);
    curv_gauss[i] = new BinaryProduct(p, k1, curv_k2[i]);

    // Create the radius variable for this node
    int iAtom = mIndex[i];
    sprintf(buffer, "R[%d]", iAtom);
    R[iAtom] = p->AddVariable(buffer, -1.0 / k1->Evaluate(), 0.1);

    // Now enforce the link between R and K1
    Expression *conR = new BinaryProduct(p, R[iAtom], k1);
    p->AddConstraint(conR, "R*kappa", -1.0, -1.0);

    // Initialize the medial axis position
    for(int j = 0; j < 3; j++)
      {
      sprintf(buffer, "M[%d,%d]", iAtom, j);
      double mval = X[i][j]->Evaluate() - R[iAtom]->Evaluate() * N[i][j]->Evaluate();
      M[iAtom][j] = p->AddVariable(buffer, mval);
      }
    }

  // ------------------------------------------------------------------------
  // Compute the initial medial atoms and radii for non-edge atoms
  // ------------------------------------------------------------------------

  for(int i = 0; i < nm; i++)
    {
    int k = mtbIndex[i].size();
    if(k > 1)
      {
      // For atoms that are bi-tangent or greater, we find R and M that
      // minimize the expression Sum_j ||Xj - r Nj - M)||^2, i.e, the radius
      // such that the medial atoms computed for each boundary atom idependently
      // are as close together as possible. The r is found as
      // r = (Sum_{ij} (Xi-Xj)^t N_i) / (k^2 - Sum_{ij} Nj^t N_i)
      double numerator = 0.0, denominator = k*k;
      SMLVec3d sumX(0.0, 0.0, 0.0), sumN(0.0, 0.0, 0.0);
      for(int q = 0; q < k; q++)
        {
        int iq = mtbIndex[i][q];
        SMLVec3d Xq = VectorEvaluate(X[iq]), Nq = VectorEvaluate(N[iq]);
        for(int p = 0; p < k; p++)
          {
          int ip = mtbIndex[i][p];
          SMLVec3d Xp = VectorEvaluate(X[ip]), Np = VectorEvaluate(N[ip]);

          numerator += dot_product(Xq - Xp, Nq);
          denominator -= dot_product(Np, Nq);
          }

        sumX += Xq; sumN += Nq;
        }

      // Compute the best fit r and m
      double v_r = numerator / denominator;
      SMLVec3d v_m = (sumX - sumN * v_r) / ((double) k);

      // Store those values
      sprintf(buffer, "R[%d]", i);
      R[i] = p->AddVariable(buffer, v_r, 0.1);

      // Initialize the medial axis position
      for(int j = 0; j < 3; j++)
        {
        sprintf(buffer, "M[%d,%d]", i, j);
        M[i][j] = p->AddVariable(buffer, v_m[j]);
        }
      }
    }

  // ------------------------------------------------------------------------
  // Add the actual medial constraints
  // ------------------------------------------------------------------------

  for(int iBnd = 0; iBnd < nb; iBnd++)
    {
    int iAtom = mIndex[iBnd];

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
      p->AddConstraint(constMedial, "X-rNM", 0, 0);
      }
    }



  // ------------------------------------------------------------------------
  // Configure the boundary and medial triangle area variables
  // ------------------------------------------------------------------------

  // Create a medial mesh that duplicates the boundary mesh
  TriangleMesh* mmesh = new TriangleMesh(*bmesh);

  // Change the triangle vertices to use m_index
  for(int i = 0; i < mmesh->triangles.size(); i++)
    {
    for(int j = 0; j < 3; j++)
      {
      mmesh->triangles[i].vertices[j] = mIndex[bmesh->triangles[i].vertices[j]];
      }
    }

  // Compute boundary triangle edges, normals, and areas
  if(tc_bnd.min_dih_angle > 0 || tc_bnd.min_tri_angle > 0 || tc_bnd.min_tri_area > 0)
    {
    ComputeTriangleAndEdgeProperties(
          p, bmesh,
          X, NT_X, taX, tc_bnd.min_tri_area, tc_bnd.min_tri_angle > 0, TEL_X);
    }

  // Compute medial triangle edges, normals, and areas.
  // TODO: this is redundant, as it computes two sets of equal normals and
  // areas for each medial triangle, one facing each boundary triangle. We
  // could avoid this by checking for opposite triangles.
  if(tc_med.min_dih_angle > 0 || tc_med.min_tri_angle > 0 || tc_med.min_tri_area > 0)
    {
    ComputeTriangleAndEdgeProperties(
          p, mmesh,
          M, NT_M, taM, tc_med.min_tri_area, tc_med.min_tri_angle > 0, TEL_M);
    }

  // ------------------------------------------------------------------------
  // Export the medial atom mesh for debugging purposes.
  // ------------------------------------------------------------------------
  SaveMedialMesh("medial_before.vtk", p, bmesh, mIndex, mtbIndex, M, R, X, taM);

  // ------------------------------------------------------------------------
  // Common code for the regularization terms
  // ------------------------------------------------------------------------
  BigSum *objBasisResidual = new BigSum(p);

  // ------------------------------------------------------------------------
  // Define the subdivision basis regularization objective
  // ------------------------------------------------------------------------
  if(regOpts.Mode == RegularizationOptions::SUBDIVISION)
    {
    // Try to obtain a parent mesh
    BCMTemplate tmplParent;
    if(!ReverseEngineerSubdivisionMesh(tmpl, tmplParent, regOpts.SubdivisionLevel))
      {
      cerr << "Unable to deduce coarse-level mesh by reversing subdivision!" << endl;
      delete p;
      return -1;
      }

    // Save the parent-level mesh
    tmplParent.Save("reverse_eng.vtk");

    // TODO - do we want to regularize R as well???

    // Define the basis variables
    SubdivisionSurface::MeshLevel &pmesh = tmplParent.bmesh;
    VarVecArray XC(pmesh.nVertices, VarVec(3, NULL));
    for(int i = 0; i <  pmesh.nVertices; i++)
      {
      for(int j = 0; j < 3; j++)
        {
        sprintf(buffer, "XC[%d][%d]", i, j);
        XC[i][j] = p->AddVariable(buffer, tmplParent.x[i][j]);
        }
      }

    // Apply the weight matrix to compute the residual
    ImmutableSparseMatrix<double> &W = tmpl.bmesh.weights;

    // Define the residual objective
    for(int i = 0; i < nb; i++)
      {
      for(int j = 0; j < 3; j++)
        {
        // Use the sum generator for each child vertex
        WeightedSumGenerator wsg(p);
        wsg.AddTerm(X[i][j], -1.0);

        // Use the weighted matrix iterator
        for(ImmutableSparseMatrix<double>::RowIterator it = W.Row(i);
            !it.IsAtEnd(); ++it)
          {
          wsg.AddTerm(XC[it.Column()][j], it.Value());
          }

        // Add the square of this residual term to the residual objective
        objBasisResidual->AddSummand(new Square(p, wsg.GenerateSum()));
        }
      }
    }

  // ------------------------------------------------------------------------
  // Define the Laplace basis regularization objective
  // ------------------------------------------------------------------------
  if(regOpts.Mode == RegularizationOptions::SPECTRAL)
    {
    // Define a basis for the surface
    MeshBasisCoefficientMapping basismap_X(bmesh, regOpts.BasisSize, 3);

    // Define a basis for the medial axis
    // TODO: to what extent is this a valid basis - we need to visualize!
    // MeshBasisCoefficientMapping basismap_M(bmesh, nBasis, 4);

    // Create the coefficient variables
    VarVecArray XC(regOpts.BasisSize, VarVec(3, NULL));
    // VarVecArray MC(nBasis, VarVec(4, NULL));
    for(int i = 0; i < regOpts.BasisSize; i++)
      {
      for(int j = 0; j < 3; j++)
        {
        XC[i][j] = p->AddVariable("XC", 0.0);
        }
      // for(int j = 0; j < 4; j++)
      //   {
      //   MC[i][j] = p->AddVariable("MC", 0.0);
      //   }
      }

    // Define the objective on the basis
    for(int iBnd = 0; iBnd < nb; iBnd++)
      {
      SMLVec3d Xfixed = xInput[iBnd];

      for(int j = 0; j < 3; j++)
        {
        BigSum *xfit = new BigSum(p);

        xfit->AddSummand(new Constant(p, Xfixed[j]));
        for(int i = 0; i < regOpts.BasisSize; i++)
          {
          xfit->AddSummand(new ScalarProduct(p, XC[i][j],
                                             basismap_X.GetBasisComponent(i, iBnd)));
          }

        // Xfit is the approximation of X using the basis
        objBasisResidual->AddSummand(
              new Square(p, new BinaryDifference(p, xfit, X[iBnd][j])));
        }
      }

    // TODO: for the time being, I took out the medial residual computation,
    // because it is unreasonable to compute it as a difference from the initial
    // state, since the initial state is likely to be bogus. Need to revive this
    // later
    /*
    for(MedialAtomIterator it = model->GetAtomIterator();
        !it.IsAtEnd(); ++it)
      {
      int iatm = it.GetIndex();

      MedialAtom &a = model->GetAtomArray()[iatm];

      for(int j = 0; j < 4; j++)
        {
        double xfixed = (j < 3) ? a.X[j] : a.R;
        Expression *xdata = (j < 3) ? M[iatm][j] : R[iatm];

        BigSum *xfit = new BigSum(p);

        xfit->AddSummand(new Constant(p, xfixed));
        for(int i = 0; i < nBasis; i++)
          {
          xfit->AddSummand(new ScalarProduct(p, MC[i][j],
                                             basismap_M.GetBasisComponent(i, iatm)));
          }

        // Xfit is the approximation of X using the basis
        objBasisResidual->AddSummand(
              new Square(p, new BinaryDifference(p, xfit, xdata)));
        }
      }
    */
    }

  // ------------------------------------------------------------------------
  // Create a total volume objective -
  // ------------------------------------------------------------------------
  /*
  BigSum *objVolume = new BigSum(p);

  // Measure the wedge volume for each boundary triangle
  VarVec wedgeVol(bmesh->triangles.size(), NULL);

  // Add the Jacobian constraints - all tetrahedra must have positive volume
  for(int i = 0; i < bmesh->triangles.size(); i++)
    {
    Triangle &T = bmesh->triangles[i];

    // Get the boundary vertices
    VarVec x[] = { X[T.vertices[0]], X[T.vertices[1]], X[T.vertices[2]] };

    // And the medial vertices
    VarVec m[] = { M[mIndex[T.vertices[0]]],
                   M[mIndex[T.vertices[1]]],
                   M[mIndex[T.vertices[2]]] };

    // There are several ways to cut a wedge into tetras, we choose one
    Expression *c1 = TetraHedronVolume(p, m[2], x[0], x[1], x[2]);
    Expression *c2 = TetraHedronVolume(p, m[1], x[0], x[1], m[2]);
    Expression *c3 = TetraHedronVolume(p, m[2], x[0], m[0], m[1]);

    // printf("Tetra Vols: %f, %f, %f\n", c1->Evaluate(), c2->Evaluate(), c3->Evaluate());

    // Save the total wedge volume for use in integration
    wedgeVol[i] = new TernarySum(p, c1, c2, c3);

    // Each tetra should have positive volume
    /--*
    p->AddConstraint(c1, 0.1, 40);
    p->AddConstraint(c2, 0.1, 40);
    p->AddConstraint(c3, 0.1, 40);
    *--/

    // Total volume integral
    objVolume->AddSummand(wedgeVol[i]);
    }
    */


  // ------------------------------------------------------------------------
  // Create the medial/boundary Jacobian constraint (normals point in same direction)
  // ------------------------------------------------------------------------
  /*
  double constJacFact = 0.1;
  for(int i = 0; i < bmesh->triangles.size(); i++)
    {
    // For this constraint, we just want the medial triangle normal and
    // the boundary triangle normal to point in the same direction!
    Expression *dp = DotProduct(p, NT_X[i], NT_M[i]);

    // Depending on the side, constrain up or down
    if(dp->Evaluate() < constJacFact)
      std::cout << "Bad Jacobian constraint: " << dp->Evaluate()
                << " in triangle " << i << std::endl;

    // Add the constraint
    p->AddConstraint(dp, "Jac", constJacFact, ConstrainedNonLinearProblem::UBINF);
    }
    */


#ifdef OLD_CODE
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
#endif

  // ------------------------------------------------------------------------
  // Create the MIB constraint
  // ------------------------------------------------------------------------
  for(int iBnd = 0; iBnd < nb; iBnd++)
    {
    int iAtom = mIndex[iBnd];

    for(EdgeWalkAroundVertex walk(bmesh, iBnd); !walk.IsAtEnd(); ++walk)
      {
      int k = walk.MovingVertexId();
      Expression *distsq = DistanceSqr(p, M[iAtom], X[k]);
      p->AddConstraint(
            new BinaryDifference(p, distsq, new Square(p, R[iAtom])),
            "MIB", 0, ConstrainedNonLinearProblem::UBINF);
      }
    }



  // ------------------------------------------------------------------------
  // Construct the total surface area objective
  // ------------------------------------------------------------------------
  /*
  BigSum *objSurfArea = new BigSum(p);
  for(int i = 0; i < bmesh->triangles.size(); i++)
    {
    // Add the area to the objective
    objSurfArea->AddSummand(taX[i]);
    }
    */

  // ------------------------------------------------------------------------
  // Solve for the circumcenter and circumradius of each boundary triangle
  // ------------------------------------------------------------------------

#ifdef CIRCUMCENTER

  // Define arrays for circumcenter, circumradius, barycentric coords of the c.c.
  VarVecArray CC(model->GetNumberOfBoundaryTriangles(), VarVec(3, NULL));
  VarVecArray CCBC(model->GetNumberOfBoundaryTriangles(), VarVec(3, NULL));
  VarVec CR(model->GetNumberOfBoundaryTriangles(), NULL);

  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {

    // Get the corners of the triangle
    VarVec XT[3] = {
      X[trit.GetBoundaryIndex(0)],
      X[trit.GetBoundaryIndex(1)],
      X[trit.GetBoundaryIndex(2)]
    };



    // We define the barycentric coordinates of the circumcenter, the circum
    // center itself and the circumradius using implicit equations, arising
    // from solving the constrained problem R=argmin(R^2) subj to. |X-C|^2 = R,
    // |Y-C|^2 = R, |Z-C|^2 = R, using Lagrange multipliers
    VarVec &lambda = CCBC[trit.GetIndex()];
    VarVec &C = CC[trit.GetIndex()];

    // Expressions for squared edge lengths
    VarVec elen2(3, NULL);
    elen2[0] = DistanceSqr(p, XT[1], XT[2]);
    elen2[1] = DistanceSqr(p, XT[2], XT[0]);
    elen2[2] = DistanceSqr(p, XT[0], XT[1]);

    // Sum of the edge lengths squared
    double elen2sum = elen2[0]->Evaluate() + elen2[1]->Evaluate() + elen2[2]->Evaluate();

    // Compute what the barycenter coordinates should be (Wikipedia)
    vnl_vector_fixed<double, 3> lambdaVal;
    for(int j = 0; j < 3; j++)
      lambdaVal[j] = elen2[j]->Evaluate() * (elen2sum - 2 * elen2[j]->Evaluate());
    lambdaVal /= lambdaVal.sum();

    // We constrain lambdas to be positive - forces acute triangles! And
    // probably setting an upper limit would work to constrain the minimal
    // angle of the triangle
    for(int j = 0; j < 3; j++)
      {
      // Wikipedia expression for the lambdas
      lambda[j] = p->AddVariable("l", lambdaVal[j], 0);
      }

    // Now we define the seven constraints that tie these variables.

    // L1 + L2 + L3 = 1
    p->AddConstraint(new TernarySum(p, lambda[0], lambda[1], lambda[2]), 1.0, 1.0);

    // C = L1 X + L2 Y + L3 Z
    for(int j = 0; j < 3; j++)
      {
      // What C should equal
      Expression *rhs = new TernarySum(p,
                                       new BinaryProduct(p, lambda[0], XT[0][j]),
                                       new BinaryProduct(p, lambda[1], XT[1][j]),
                                       new BinaryProduct(p, lambda[2], XT[2][j]));
      // Add the variable and the constraint
      C[j] = p->AddExpressionAsConstrainedVariable(rhs);
      }

    // Compute squared distances from the vertices to the center
    VarVec XCd2(3, NULL);
    XCd2[0] = DistanceSqr(p, XT[0], C);
    XCd2[1] = DistanceSqr(p, XT[1], C);
    XCd2[2] = DistanceSqr(p, XT[2], C);

    // Create radius variable
    Expression *Rc = CR[trit.GetIndex()] = p->AddVariable("Rc", sqrt(XCd2[0]->Evaluate()), 0);

    // |C-X|^2 = R^2
    for(int j = 0; j < 3; j++)
      {
      p->AddConstraint(new BinaryDifference(p, XCd2[j], new Square(p,Rc)), 0.0, 0.0);
      }

    /*

    // Add variables for cotan alpha
    for(int j = 0; j < 3; j++)
      {
      // Edge length opposite vertex j
      Expression *elen = p->AddVariable("len", sqrt(elen2[j]->Evaluate()), 0);

      // Make it equal its square!
      Expression *conEdge = new BinaryDifference(p, new Square(p, elen), elen2[j]);
      p->AddConstraint(conEdge, 0, 0);

      // The distance from circumcenter to edge
      Expression *height = p->AddVariable("hgt", sqrt(Rc->Evaluate()*Rc->Evaluate() - 0.25 * elen->Evaluate() * elen->Evaluate()), 0);

      // Simple linking expressions: 0.25 * elen^2 + h^2 = r^2

      Expression *conPytha = new TernarySum(p,
                                            new ScalarProduct(p, new Square(p, elen), 0.25),
                                            new Square(p, height),
                                            new Negation(p, new Square(p, Rc)));

      p->AddConstraint(conPytha, 0, 0);

      // The cotan we care about!
      Expression *cotA = p->AddVariable("cotA", height->Evaluate() / elen->Evaluate(), 0);

      // cot-alpha * elen = height
      p->AddConstraint(new BinaryDifference(p,
                                            new BinaryProduct(p, cotA, elen),
                                            height), 0, 0);

      printf("CotA[%d] : %f\n", j, cotA->Evaluate());
      }
      */

    }

#endif

#define ASPECTRATIO 1

#ifdef ASPECTRATIO_OLD

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

/*
    // Create a variable for the cosine of each angle in the triangle
    for(int j = 0; j < 3; j++)
      {
      int k0 = trit.GetBoundaryIndex(j);
      int k1 = trit.GetBoundaryIndex((j + 1) % 3);
      int k2 = trit.GetBoundaryIndex((j + 2) % 3);

      // cos(theta) = dot(u,v) / (|u| * |v|)
      Expression *UdotV = DotProduct(
            p,
            VectorApplyPairwise<BinaryDifference>(p, X[k1], X[k0]),
            VectorApplyPairwise<BinaryDifference>(p, X[k2], X[k0]));

      Expression *lenU = edge[(j+1) % 3];
      Expression *lenV = edge[(j+2) % 3];

      // Add the cosine of alpha with constraints on range
      Expression *cosAlpha = p->AddVariable(
            "cos", UdotV->Evaluate() / (lenU->Evaluate() * lenV->Evaluate()),
            0, acos(vnl_math::pi * 15 / 180));

      // Create a constraint linking the cosine
      Expression *constr = new BinaryDifference(
            p, new TernaryProduct(p, lenU, lenV, cosAlpha), UdotV);

      p->AddConstraint(constr, 0, 0);
      }
*/

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

#endif

  // ------------------------------------------------------------------------
  // Add a constraint on minimal angle of boundary triangles
  // ------------------------------------------------------------------------
  if(tc_bnd.min_tri_angle > 0)
    {
    double min_angle = vnl_math::pi * tc_bnd.min_tri_angle / 180;
    double max_csc = 1.0 / sin(min_angle);

    for(int k = 0; k < bmesh->triangles.size(); k++)
      {
      for(int d = 0; d < 3; d++)
        {
        Expression *l1 = TEL_X[k][(d + 1) % 3];
        Expression *l2 = TEL_X[k][(d + 2) % 3];
        double v_csc_alpha = (l1->Evaluate() * l2->Evaluate()) / (2 * taX[k]->Evaluate());

        // Create a constrainted variable for the angle
        Expression *csc_alpha = p->AddVariable("cscAlpha", v_csc_alpha,
                                               ConstrainedNonLinearProblem::LBINF,
                                               max_csc);

        // Tie using the constraints, based on fmla 2*At*csc(alpha) = l1 * l2;
        Expression *con = new BinaryDifference(
              p,
              new ScalarProduct(p, new BinaryProduct(p, taX[k], csc_alpha), 2.0),
              new BinaryProduct(p, l1, l2));

        if(fabs(con->Evaluate()) > 1e-6)
          std::cout << "Con-CSC: " << con->Evaluate() << std::endl;

        p->AddConstraint(con, "CSC", 0, 0);
        }
      }
    }

  // ------------------------------------------------------------------------
  // Add a constraint on minimal angle of medial triangles
  // ------------------------------------------------------------------------
  if(tc_med.min_tri_angle > 0)
    {
    double min_angle_med = vnl_math::pi * tc_med.min_tri_angle / 180;
    double max_csc_med = 1.0 / sin(min_angle_med);

    for(int k = 0; k < mmesh->triangles.size(); k++)
      {
      for(int d = 0; d < 3; d++)
        {
        Expression *l1 = TEL_M[k][(d + 1) % 3];
        Expression *l2 = TEL_M[k][(d + 2) % 3];
        double v_csc_alpha = (l1->Evaluate() * l2->Evaluate()) / (2 * taM[k]->Evaluate());

        // Create a constrainted variable for the angle
        Expression *csc_alpha = p->AddVariable("cscAlpha", v_csc_alpha,
                                               ConstrainedNonLinearProblem::LBINF,
                                               max_csc_med);

        // Tie using the constraints, based on fmla 2*At*csc(alpha) = l1 * l2;
        Expression *con = new BinaryDifference(
              p,
              new ScalarProduct(p, new BinaryProduct(p, taM[k], csc_alpha), 2.0),
              new BinaryProduct(p, l1, l2));

        if(fabs(con->Evaluate()) > 1e-6)
          std::cout << "Con-Med-CSC: " << con->Evaluate() << std::endl;

        p->AddConstraint(con, "Med-CSC", 0, 0);
        }
      }
    }


  // ------------------------------------------------------------------------
  // Minimize the edge length
  // ------------------------------------------------------------------------
  // TODO: add a constraint on the angle between adjacent boundary edges, b/c
  // that seems to become quite small sometimes and setting a minimum on it might
  // work better than penalizing the total length
  BigSum *objEdgeLength = new BigSum(p);
  if(regOpts.EdgeLengthWeight > 0.0)
    {
    // Add all the crest edges
    for(int k = 0; k < bmesh->triangles.size(); k++)
      {
      for(int d = 0; d < 3; d++)
        {
        int v1 = bmesh->triangles[k].vertices[(d + 1) % 3];
        int v2 = bmesh->triangles[k].vertices[(d + 2) % 3];
        if(v1 < v2)
          {
          // Are these edge vertices?
          if(mtbIndex[mIndex[v1]].size() == 1 && mtbIndex[mIndex[v2]].size() == 1)
            {
            objEdgeLength->AddSummand(TEL_X[k][d]);
            }
          }
        }
      }
    }

  // ------------------------------------------------------------------------
  // Constrain the dihedral angle
  // ------------------------------------------------------------------------
  if(tc_bnd.min_dih_angle > 0)
    {
    double non_edge_thresh = cos(vnl_math::pi * (180. - tc_bnd.min_dih_angle) / 180.);
    double edge_thresh = -0.9;

    for(int k = 0; k < bmesh->triangles.size(); k++)
      {
      for(int d = 0; d < 3; d++)
        {
        int kopp = bmesh->triangles[k].neighbors[d];
        // One side only!
        if(kopp != NOID && k < kopp)
          {
          // Is this a boundary edge?
          int v1 = bmesh->triangles[k].vertices[(d+1)%3];
          int v2 = bmesh->triangles[k].vertices[(d+2)%3];
          bool isBnd = (mtbIndex[mIndex[v1]].size() == 1) && (mtbIndex[mIndex[v2]].size() == 1);

          Expression *da = DotProduct(p, NT_X[k], NT_X[kopp]);
          p->AddConstraint(da, "DA", isBnd ? edge_thresh : non_edge_thresh, ConstrainedNonLinearProblem::UBINF);

          }
        }
      }
    }

  // ------------------------------------------------------------------------
  // Constrain the dihedral angle of the medial axis
  // ------------------------------------------------------------------------
  if(tc_med.min_dih_angle > 0)
    {
    double non_edge_thresh = cos(vnl_math::pi * (180. - tc_med.min_dih_angle) / 180.);
    for(int k = 0; k < mmesh->triangles.size(); k++)
      {
      for(int d = 0; d < 3; d++)
        {
        int kopp = mmesh->triangles[k].neighbors[d];
        // One side only!
        if(kopp != NOID && k < kopp)
          {
          // Is this a boundary edge?
          int v1 = mmesh->triangles[k].vertices[(d+1)%3];
          int v2 = mmesh->triangles[k].vertices[(d+2)%3];
          bool isBnd = (mtbIndex[v1].size() == 1) && (mtbIndex[v2].size() == 1);

          Expression *da = DotProduct(p, NT_M[k], NT_M[kopp]);
          if(!isBnd)
            {
            if(da->Evaluate() < 0.9)
              printf("MDA: (%d, %d) via edge (%d, %d): %f\n", k, kopp, v1, v2, da->Evaluate());
            p->AddConstraint(da, "DA-Med", non_edge_thresh, ConstrainedNonLinearProblem::UBINF);
            }
          }
        }
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

#ifdef USE_DICE
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

#endif // USE_DICE

  // This simple objective is the displacement from the input boundary points.
  // It allows us to resolve all the constraints without fitting to new data.
  // It's a good sanity check to make sure the model is valid
  BigSum *objDisplacement = new BigSum(p);
  for(int i = 0; i < nb; i++)
    {
    objDisplacement->AddSummand(
          new TernaryGradientMagnitudeSqr(
            p,
            new BinaryDifference(p, X[i][0], new Constant(p, xLSQTarget[i][0])),
            new BinaryDifference(p, X[i][1], new Constant(p, xLSQTarget[i][1])),
            new BinaryDifference(p, X[i][2], new Constant(p, xLSQTarget[i][2]))));
    }

  // ------------------------------------------------------------------------
  // Derive the final objective
  // ------------------------------------------------------------------------

  // obj->AddSummand(objSqDist);
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, scaleRecip), objRecipSqDist));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 20), objImageMatch));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 0.01), objSurfArea));


  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 0.5), objSimpleBending));

  // Set of objectives for fitting the model to itself
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 1), objDisplacement));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 0.1), objSimpleBending));
  // obj->AddSummand(new BinaryProduct(p, new Constant(p, 0.02), objSurfArea));

  // Dice!
  /*
  obj->AddSummand(new ScalarProduct(p, objDice, 4 * (nb + nm)));
  obj->AddSummand(new ScalarProduct(p, objBasisResidual, 1));
  */


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

  // Save the true objective
  // SaveGradient(p, model, X, obj, "grad_obj_before.vtk");



#ifdef USE_DICE
  std::cout << "Dice objective: " << 1- objDice->Evaluate() << std::endl;
  std::cout << "Image match: " << objObjectIntegral->Evaluate() << std::endl;
  std::cout << "Volume integral: " << objVolumeIntegral->Evaluate() << std::endl;
  std::cout << "Target volume: " << volTarget << std::endl;

  // Save the sample and image values at samples
  SaveSamples(sampleX, sampleF, "samples_before.vtk");

  // Plot the gradient of the problem
  SaveGradient(p, model, X, objDice, "grad_obj_dice_before.vtk");
  SaveGradient(p, model, X, objSimpleBending, "grad_obj_bend_before.vtk");
#endif


  // Test some derivatives;
  // DerivativeTest(p, 1000);

  // Solve the problem
  SmartPtr<IPOptProblemInterface> ip = new IPOptProblemInterface(p);

  // Create a file for dumping the constraint info
  FILE *fnConstraintDump = fopen("constraint_dump.txt", "wt");
  ip->log_constraints(fnConstraintDump);

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
  app->Options()->SetStringValue("linear_solver", regOpts.Solver);
  // app->Options()->SetNumericValue("mu_init", 1e-3);
  // app->Options()->SetNumericValue("mu_target", 1e-5);
  // app->Options()->SetStringValue("mu_strategy", "adaptive");
  // app->Options()->SetStringValue("output_file", "ipopt.out");

  app->Options()->SetIntegerValue("max_iter", 200);
  if(!regOpts.UseHessian)
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

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

  // Combine the regularization objectives
  WeightedSumGenerator wsgObj(p);
  wsgObj.AddTerm(objBasisResidual, regOpts.Weight);
  wsgObj.AddTerm(objEdgeLength, regOpts.EdgeLengthWeight);
  Expression *objRegCombined = wsgObj.GenerateSum();

  // Try just fitting the boundary data
  ObjectiveCombiner ocfit(p);
  ocfit.AddObjectiveTerm(objBasisResidual, "Penalty (Basis Fit)", regOpts.Weight);
  ocfit.AddObjectiveTerm(objEdgeLength, "Penalty (Edge Length)", regOpts.EdgeLengthWeight);
  ocfit.AddObjectiveTerm(objDisplacement, "Displacement Objective", 1);
  Expression *obj = ocfit.GetTotalObjective();

  // Configure the problem
  p->SetObjective(obj);
  p->SetupProblem(true);

  // Evaluate the objective
  ocfit.PrintReport();

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(GetRawPtr(ip));

  // Evaluate the objective
  ocfit.PrintReport();

  // Remove the extension fromthe filename
  std::string fnOutBase = itksys::SystemTools::GetFilenameWithoutLastExtension(fnOutput);

  // Save the current state
  sprintf(buffer, "%s_fit2tmp_bnd.vtk", fnOutBase.c_str());
  SaveBoundaryMesh(buffer, p, bmesh, mIndex, mtbIndex, X, N, R);

  sprintf(buffer, "%s_fit2tmp_med.vtk", fnOutBase.c_str());
  SaveMedialMesh(buffer, p, bmesh, mIndex, mtbIndex, M, R, X, taM);

  // Continue if in ICP mode
  if(action == ACTION_FIT_TARGET)
    {
    // Create the closest point finder
    ClosestPointMatcher cpmatcher(target, 32);

    // The scale between the two distance functions
    double recip_scale = nb * 1.0 / cpmatcher.GetNumberOfTargetPointsUsed();


    // Repeat this several times
    Expression *objSqDist, *objRecipSqDist;
    for(int i = 0; i < 5; i++)
      {
      // ------------------------------------------------------------------------
      // Construct the first part of the objective function
      // ------------------------------------------------------------------------
      objSqDist = ComputeDistanceToMeshObjective(p, &cpmatcher, X);

      // ------------------------------------------------------------------------
      // Construct an opposite objective: distance from a selected set of points
      // to the closest point on the medial mesh
      // ------------------------------------------------------------------------
      objRecipSqDist = ComputeDistanceToModelObjective(p, &cpmatcher, bmesh, X);

      // ICP-style
      ObjectiveCombiner ocicp(p);
      ocicp.AddObjectiveTerm(objBasisResidual, "Penalty (Basis Fit)", regOpts.Weight);
      ocicp.AddObjectiveTerm(objEdgeLength, "Penalty (Edge Length)", regOpts.EdgeLengthWeight);
      ocicp.AddObjectiveTerm(objSqDist, "DistanceSqr to Target", 1);
      ocicp.AddObjectiveTerm(objRecipSqDist, "DistanceSqr to Model", recip_scale);
      obj = ocicp.GetTotalObjective();

      // Configure the problem
      p->SetObjective(obj);
      p->SetupProblem(true);

      // Evaluate the objective
      ocicp.PrintReport();

      // Ask Ipopt to solve the problem
      status = app->OptimizeTNLP(GetRawPtr(ip));

      // Evaluate the objective
      ocicp.PrintReport();

      // Save the current state
      sprintf(buffer, "%s_icp_%02d_bnd.vtk", fnOutBase.c_str(), i);
      SaveBoundaryMesh(buffer, p, bmesh, mIndex, mtbIndex, X, N, R);

      sprintf(buffer, "%s_icp_%02d_med.vtk", fnOutBase.c_str(), i);
      SaveMedialMesh(buffer, p, bmesh, mIndex, mtbIndex, M, R, X, taM);
      }

  #ifdef USE_DICE
    std::cout << "Image match: " << objObjectIntegral->Evaluate() << std::endl;
    std::cout << "Volume integral: " << objVolumeIntegral->Evaluate() << std::endl;
    std::cout << "Dice objective: " << 1- objDice->Evaluate() << std::endl;
    std::cout << "Target volume: " << volTarget << std::endl;
  #endif

    if (status == Solve_Succeeded) {
      printf("\n\n*** The problem solved!\n");
      }
    else {
      printf("\n\n*** The problem FAILED!\n");
      }

    // Test some derivatives;
    // DerivativeTest(p, 1000);

    // Save the result as a boundary mesh
    // SaveBoundaryMesh("result_bnd.vtk", p, bmesh, mIndex, mtbIndex, X, N, R);
    // SaveMedialMesh("result_med.vtk", p, bmesh, mIndex, M, R);

    SaveGradient(p, X, obj, "grad_obj_after.vtk");
    SaveGradient(p, X, objSqDist, "grad_obj_sqdist_after.vtk");
    SaveGradient(p, X, new ScalarProduct(p, objRecipSqDist, recip_scale), "grad_obj_recipsqdist_after.vtk");
    SaveGradient(p, X, new ScalarProduct(p, objBasisResidual, regOpts.Weight), "grad_obj_residual_after.vtk");

    }

#ifdef USE_DICE
  // Save the sample and image values at samples
  SaveSamples(sampleX, sampleF, "samples_after.vtk");

  // Plot the gradient of the problem
  SaveGradient(p, X, objDice, "grad_obj_dice_after.vtk");
  SaveGradient(p, X, objBasisResidual, "grad_obj_residual_after.vtk");
#endif


#ifdef CIRCUMCENTER
  SaveCircumcenterMesh(CC, CR, CCBC);
#endif

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.
  // delete p;

  fclose(fnConstraintDump);
  return (int) status;
}
