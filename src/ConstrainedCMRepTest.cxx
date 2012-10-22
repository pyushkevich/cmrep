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

#include "vtkPolyData.h"
#include "vtkCellLocator.h"
#include <vector>


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



int main(int argc, char *argv[])
{
  // The first parameter is the cm-rep to start from
  char *cmrepfile = argv[1];

  // The second is the VTK mesh for ICP
  char *targetmesh = argv[2];

  // Load and process the mrep
  MedialPDE mrep(cmrepfile);
  GenericMedialModel *model = mrep.GetMedialModel();

  // Load the target mesh
  vtkPolyData *target = ReadVTKMesh(targetmesh);

  // Create the optimization problem
  ConstrainedNonLinearProblem *p = new ConstrainedNonLinearProblem();

  // For each point on the model, find the closest target point
  std::vector<SMLVec3d> targetPoint = find_closest(model, target);

  // Add the variables to optimize over
  typedef std::vector<Variable *> VarVec;
  typedef std::vector<VarVec> VarVecArray;
  VarVecArray X(model->GetNumberOfBoundaryPoints(), VarVec(3, NULL));

  for(MedialBoundaryPointIterator it = model->GetBoundaryPointIterator();
      !it.IsAtEnd(); ++it)
    {
    SMLVec3d x = GetBoundaryPoint(it, model->GetAtomArray()).X;
    for(int j = 0; j < 3; j++)
      {
      ostringstream oss;
      oss << "X[" << it.GetIndex() << "," << j << "]";
      X[it.GetIndex()][j] = p->AddVariable(oss.str(), x[j], 0, 256);
      }
    }

  // Add a variable for each triangle
  VarVec A(model->GetNumberOfBoundaryTriangles(), NULL);
  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    ostringstream oss;
    oss << "A[" << trit.GetIndex() << "]";

    double a = TriangleArea(GetBoundaryPoint(trit, model->GetAtomArray(), 0).X,
                            GetBoundaryPoint(trit, model->GetAtomArray(), 1).X,
                            GetBoundaryPoint(trit, model->GetAtomArray(), 2).X);

    a += 0.2 * (rand() * 1.0 / RAND_MAX - 0.5);

    A[trit.GetIndex()] = p->AddVariable(oss.str(), a, 0);
    }

  // Construct the first part of the objective function
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

  // Construct the total surface area objective
  BigSum *objSurfArea = new BigSum(p);
  for(MedialBoundaryTriangleIterator trit = model->GetBoundaryTriangleIterator();
      !trit.IsAtEnd(); ++trit)
    {
    // Add the area to the objective
    objSurfArea->AddSummand(A[trit.GetIndex()]);

    // Create the corresponding constraint

    // Write an expression for the triangle area
    std::vector<Expression *> AB(3), AC(3);
    for(int j = 0; j < 3; j++)
      {
      AB[j] = new BinaryDifference(p,
                                   X[trit.GetBoundaryIndex(1)][j],
                                   X[trit.GetBoundaryIndex(0)][j]);
      AC[j] = new BinaryDifference(p,
                                   X[trit.GetBoundaryIndex(2)][j],
                                   X[trit.GetBoundaryIndex(0)][j]);
      }

    TernaryGradientMagnitudeSqr *tgm = new TernaryGradientMagnitudeSqr(
          p,
          new BinaryDifference(p,
                               new BinaryProduct(p, AB[1], AC[2]),
                               new BinaryProduct(p, AB[2], AC[1])),
          new BinaryDifference(p,
                               new BinaryProduct(p, AB[2], AC[0]),
                               new BinaryProduct(p, AB[0], AC[2])),
          new BinaryDifference(p,
                               new BinaryProduct(p, AB[0], AC[1]),
                               new BinaryProduct(p, AB[1], AC[0])));

    // Create the constraint (A^2 = (area)^2)
    Expression *constr =
        new BinaryDifference(p,
                             new Square(p, A[trit.GetIndex()]),
                             tgm);

    p->AddConstraint(constr, 0, 0);
    }

  // Derive the final objective
  Expression *obj =
      new BinarySum(p,
                    objSqDist,
                    new BinaryProduct(p, new Constant(p, 0.1), objSurfArea));

  // Evaluate the objective
  std::cout << "Initial MSD to target: " << objSqDist->Evaluate() << std::endl;
  std::cout << "Initial surface area: " << objSurfArea->Evaluate() << std::endl;
  std::cout << "Initial objective: " << obj->Evaluate() << std::endl;

  // Configure the problem
  p->SetObjective(obj);
  p->SetupProblem();

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
  // app->Options()->SetNumericValue("tol", 1e-9);
  // app->Options()->SetStringValue("mu_strategy", "adaptive");
  // app->Options()->SetStringValue("output_file", "ipopt.out");
  app->Options()->SetIntegerValue("max_iter", 1000);
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

  std::cout << "Final MSD to target: " << objSqDist->Evaluate() << std::endl;
  std::cout << "Final surface area: " << objSurfArea->Evaluate() << std::endl;
  std::cout << "Final objective: " << obj->Evaluate() << std::endl;

  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
  }
  else {
    printf("\n\n*** The problem FAILED!\n");
  }

  // Save the mesh for output
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



  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.
  delete p;

  return (int) status;
}
