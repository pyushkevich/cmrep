#include "IPOptProblemInterface.h"
#include "ScriptInterface.h"
#include "BasisFunctions2D.h"
#include "BruteForceSubdivisionMedialModel.h"
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
#include "IpIpoptAlg.hpp"
#include "MedialAtomGrid.h"
#include "VTKMeshBuilder.h"
#include "util/ReadWriteVTK.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"

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

#include "tetgen.h"

#include "vtkPolyData.h"
#include "vtkCellLocator.h"
#include "vtkSmartPointer.h"
#include <vtkPointLocator.h>

#include <vector>
#include <map>
#include <utility>
#include "itk_to_nifti_xform.h"
#include "itksys/SystemTools.hxx"

#include "QCQPProblem.h"
#include "IPOptQCQPProblemInterface.h"

#include <vnl_random.h>

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

  SMLVec3d FindClosestToTarget(const SMLVec3d &x);

  std::vector<MatchLocation> FindClosestToSource(TriangleMesh *mesh, std::vector<SMLVec3d> &X);

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

SMLVec3d ClosestPointMatcher::FindClosestToTarget(const SMLVec3d &x)
{
  double xs[3], d2, d;
  int subid;
  vtkIdType cellid;

  m_TargetLocator->FindClosestPoint(x.data_block(), xs, cellid, subid, d2);

  return SMLVec3d(xs);
}

std::vector<ClosestPointMatcher::MatchLocation>
ClosestPointMatcher::FindClosestToSource(TriangleMesh *mesh, std::vector<SMLVec3d> &X)
{
  // Create a VTK points object
  vtkSmartPointer<vtkPoints> out_pts = vtkSmartPointer<vtkPoints>::New();
  out_pts->Allocate(X.size());

  for(int i = 0; i < X.size(); i++)
    {
    out_pts->InsertNextPoint(X[i][0], X[i][1], X[i][2]);
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
    SMLVec3d A = X[T.vertices[0]];
    SMLVec3d B = X[T.vertices[1]];
    SMLVec3d C = X[T.vertices[2]];

    vnl_matrix<double> W(3, 2);
    SMLVec3d delBA = B - A, delCA = C - A;
    W.set_column(0, delBA.as_ref());
    W.set_column(1, delCA.as_ref());
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


/**
 * Optimal Mass Transport using Generalized Sinkhorn Iteration
 * Code based on Karlsson and Ringh paper (https://arxiv.org/pdf/1612.02273.pdf)
 *
 * This class is meant to be reused for multiple runs of Sinkhorn iteration
 * with the same target mesh but different starting points.
 */
class SinkhornIteration
{
public:
  friend class BCMRepQuadraticProblemBuilder;

  /**
   * The constructor takes a source mesh, a target mesh, and a set of parameters
   */
  SinkhornIteration(TriangleMesh *source, TriangleMesh *target, const std::vector<SMLVec3d> &Xtarget)
    : m_Source(source), m_Target(target)
    {
    unsigned int n0 = m_Source->triangles.size(), n1 = m_Target->triangles.size();

    // First, compute the z-transform parameters, so that the target shape is
    // centered on the origin and has radius 1
    SMLVec3d trg_sum(0.0), trg_sum_sq(0.0);
    for(auto xt : Xtarget)
      {
      trg_sum += xt;
      trg_sum_sq += element_product(xt, xt);
      }
    m_TrgCenter = trg_sum * (1.0 / Xtarget.size());
    SMLVec3d trg_var = trg_sum_sq * (1.0 / Xtarget.size()) - element_product(m_TrgCenter, m_TrgCenter);
    m_TrgScale = sqrt((trg_var[0]+trg_var[1]+trg_var[2]) / 3.0);

    // From the target mesh, extract triangle centers, normals, and areas
    m_ZTrgCenters.resize(n1);
    m_ZTrgNormals.resize(n1);
    m_ZTrgAreas.set_size(n1);
    m_Mu1.set_size(n1);
    m_U1.set_size(n1); m_U1.fill(1.0);
    m_ZTrgTotalArea = 0;
    for(int kt = 0; kt < n1; kt++)
      {
      Triangle &t = m_Target->triangles[kt];
      size_t *v = t.vertices;
      SMLVec3d X0 = (Xtarget[v[0]] - m_TrgCenter) / m_TrgScale;
      SMLVec3d X1 = (Xtarget[v[1]] - m_TrgCenter) / m_TrgScale;
      SMLVec3d X2 = (Xtarget[v[2]] - m_TrgCenter) / m_TrgScale;
      SMLVec3d Xu_cross_Xv = vnl_cross_3d(X1-X0, X2-X0);
      double mag_Xu_cross_Xv = Xu_cross_Xv.magnitude();

      m_ZTrgCenters[kt] = (X0 + X1 + X2) / 3.0;
      m_ZTrgNormals[kt] = Xu_cross_Xv / mag_Xu_cross_Xv;
      m_ZTrgAreas[kt] = 0.5 * mag_Xu_cross_Xv;
      m_ZTrgTotalArea += m_ZTrgAreas[kt];
      }

    // Target mass vector
    m_Mu1 = m_ZTrgAreas / m_ZTrgTotalArea;

    // Initialize the source arrays but without assignment
    m_ZSrcCenters.resize(n0);
    m_ZSrcNormals.resize(n0);
    m_ZSrcAreas.set_size(n0);
    m_Mu0.set_size(n0);
    m_U0.set_size(n0); m_U0.fill(1.0);

    // Initialize the matrix K
    m_C.set_size(n0, n1);
    m_K.set_size(n0, n1);
    }

  void PrintLoss(int iter)
    {
    unsigned int n0 = m_Source->triangles.size(), n1 = m_Target->triangles.size();

    // Check if the iterations are correct
    double loss = 0.0, entropy = 0.0;
    vnl_vector<double> mu0_test = m_Mu0, mu1_test = m_Mu1;
    for(unsigned int i = 0; i < n0; i++)
      {
      for(unsigned int j = 0; j < n1; j++)
        {
        // Compute the actual joint probability value for i,j
        double m_ij = m_U0(i) * m_K(i, j) * m_U1(j);
        loss += m_C(i,j) * m_ij;
        entropy -= m_ij < 1.e-15 ? 0.0 : m_ij * log(m_ij);
        mu0_test[i] -= m_ij;
        mu1_test[j] -= m_ij;
        }
      }

    printf("Sinkhorn Iter %03d   MC: %8g   Entropy: %8g   Err0: %8g   Err1: %8g\n",
           iter, loss, entropy, mu0_test.inf_norm(), mu1_test.inf_norm());
    }

  /**
   * Perform Sinkhorn iteration
   */
  void PerformIteration(const std::vector<SMLVec3d> &XSource,
                        double alpha = 1.0, double eps = 0.01, int niter = 100,
                        bool reset_starting_point = true)
    {
    // Store alpha
    m_Alpha = alpha;
    unsigned int n0 = m_Source->triangles.size(), n1 = m_Target->triangles.size();

    // Compute the source properties
    m_ZSrcTotalArea = 0;
    for(int kt = 0; kt < n0; kt++)
      {
      Triangle &t = m_Source->triangles[kt];
      size_t *v = t.vertices;
      SMLVec3d X0 = (XSource[v[0]] - m_TrgCenter) / m_TrgScale;
      SMLVec3d X1 = (XSource[v[1]] - m_TrgCenter) / m_TrgScale;
      SMLVec3d X2 = (XSource[v[2]] - m_TrgCenter) / m_TrgScale;
      SMLVec3d Xu_cross_Xv = vnl_cross_3d(X1-X0, X2-X0);
      double mag_Xu_cross_Xv = Xu_cross_Xv.magnitude();

      m_ZSrcCenters[kt] = (X0 + X1 + X2) / 3.0;
      m_ZSrcNormals[kt] = Xu_cross_Xv / mag_Xu_cross_Xv;
      m_ZSrcAreas[kt] = 0.5 * mag_Xu_cross_Xv;
      m_ZSrcTotalArea += m_ZSrcAreas[kt];
      }

    // Compute the mass vectors mu and the initial scaling vectors u
    m_Mu0 = m_ZSrcAreas / m_ZSrcTotalArea;

    // Compute the matrix K. TODO: this can be threaded
    for(int i = 0; i < n0; i++)
      {
      for(int j = 0; j < n1; j++)
        {
        double dist = 0.0;
        double cos_theta = 0.0;
        for(unsigned int d = 0; d < 3; d++)
          {
          double del = m_ZSrcCenters[i][d] - m_ZTrgCenters[j][d];
          dist += del * del;
          cos_theta += m_ZSrcNormals[i][d] * m_ZTrgNormals[j][d];
          }
        m_C(i,j) = dist * (1.0 + alpha * (1.0 - cos_theta));
        m_K(i,j) = exp(-m_C(i,j) / eps);
        }
      }

    // Not clear if we should reset the u vectors or reuse them
    if(reset_starting_point)
      {
      m_U0.fill(1.0);
      m_U1.fill(1.0);
      }

    // Perform the actual Sinkhorn iteration
    for(unsigned int iter = 0; iter < niter; iter++)
      {
      // Print the current loss
      if(iter % 20 == 0)
        PrintLoss(iter);

      // Compute u0 <- mu0 ./ K * u1
      for(unsigned int i = 0; i < n0; i++)
        {
        double f = 0;
        for(unsigned int j = 0; j < n1; j++)
          f += m_K.get(i,j) * m_U1[j];
        m_U0[i] = m_Mu0[i] / f;
        }

      // Compute u1 <- mu1 ./ K' * u0
      for(unsigned int j = 0; j < n1; j++)
        {
        double f = 0;
        for(unsigned int i = 0; i < n0; i++)
          {
          f += m_K.get(i,j) * m_U0[i];
          }
        m_U1[j] = m_Mu1[j] / f;
        }
      }

    PrintLoss(niter);
    }

  /** Compute the matched target locations for source faces */
  void ComputeTriangleCenterMatches(std::vector<SMLVec3d> &XMatches)
    {
    unsigned int n0 = m_Source->triangles.size(), n1 = m_Target->triangles.size();
    XMatches.resize(n0);
    for(unsigned int i = 0; i < n0; i++)
      {
      SMLVec3d match(0.0);
      double mass = 0.0;
      for(unsigned int j = 0; j < n1; j++)
        {
        // Compute the actual joint probability value for i,j
        double m_ij = m_U0(i) * m_K(i, j) * m_U1(j);

        // Use as weight against target centers
        match += m_ij * m_ZTrgCenters[j];
        mass += m_ij;
        }

      XMatches[i] = match * (m_TrgScale / mass) + m_TrgCenter;
      }
    }

protected:
  TriangleMesh *m_Source, *m_Target;

  // Target center and scale, used for z-transforming data
  SMLVec3d m_TrgCenter;
  double m_TrgScale;

  // Source and target triangle centers and normals, z-transformed
  std::vector<SMLVec3d> m_ZTrgCenters, m_ZTrgNormals, m_ZSrcCenters, m_ZSrcNormals;

  // Source and target triangle areas, z-transformed
  vnl_vector<double> m_ZTrgAreas, m_ZSrcAreas;

  // Total surface area
  double m_ZTrgTotalArea, m_ZSrcTotalArea;

  // The matrix K
  vnl_matrix<double> m_C, m_K;

  // The mass vectors and the scaling vectors
  vnl_vector<double> m_Mu0, m_Mu1, m_U0, m_U1;

  // The alpha used for the cost computation
  double m_Alpha = 1.0;
};




struct TriangleProperties
{
  qcqp::VariableRef<2> triangle_normals;
  qcqp::VariableRef<1> triangle_areas;
};

TriangleProperties
ComputeTriangleProperties(
    qcqp::Problem &p, TriangleMesh *mesh, qcqp::VariableRef<2> &X,
    const std::string &name_prefix, double minArea)
{
  // Number of triangles
  int nt = mesh->triangles.size();

  // Array names
  std::string nm_nt = name_prefix + "NT", nm_at = name_prefix + "AT";
  std::string nm_con_ntna = name_prefix + "NTNA", nm_con_ntnt = name_prefix + "NTNT";

  // Create the triangle properties
  qcqp::VariableRef<2> tri_normals = p.AddVariableTensor<2>(nm_nt.c_str(), {nt, 3});
  qcqp::VariableRef<1> tri_areas = p.AddVariableTensor<1>(nm_at.c_str(), {nt}, minArea);

  // Iterate over all the triangles in this mesh
  for(int it = 0; it < nt; it++)
    {
    // Here is a triangle
    Triangle &t = mesh->triangles[it];
    int v0 = t.vertices[0], v1 = t.vertices[1], v2 = t.vertices[2];

    // First compute the actual normal and actual area
    SMLVec3d Xu = X(v1).as_vector<3>() - X(v0).as_vector<3>();
    SMLVec3d Xv = X(v2).as_vector<3>() - X(v0).as_vector<3>();
    SMLVec3d Xu_cross_Xv = vnl_cross_3d(Xu, Xv);
    double v_area = Xu_cross_Xv.magnitude();
    SMLVec3d v_normal = Xu_cross_Xv / v_area;
    v_area *= 0.5;

    // Assign the variables initial values
    tri_normals(it).from_vector(v_normal);
    tri_areas(it) = v_area;

    // Constraint that ties the area and the cross-product of the edges.
    //   (a x b) + (b x c) + (c x a) = 2 * A * N
    // where a,b,c are the triangle vertices
    for(int j = 0; j < 3; j++)
      {
      auto &con_an = p.AddConstraint(nm_con_ntna.c_str(), 0.0, 0.0);

      // For each pair of vertices, we want to add the j-th component of their
      // cross product
      int j1 = (j + 1) % 3, j2 = (j + 2) % 3;
      con_an.A(X(v0,j1), X(v1,j2)) = -1.0;
      con_an.A(X(v0,j2), X(v1,j1)) = 1.0;
      con_an.A(X(v1,j1), X(v2,j2)) = -1.0;
      con_an.A(X(v1,j2), X(v2,j1)) = 1.0;
      con_an.A(X(v2,j1), X(v0,j2)) = -1.0;
      con_an.A(X(v2,j2), X(v0,j1)) = 1.0;

      // Place 2 * A * N in the right hand side
      con_an.A(tri_normals(it,j), tri_areas(it)) = 2.0;
      }

    // Constraint that makes the normal have unit length
    auto &con_nn = p.AddConstraint(nm_con_ntnt.c_str(), 1.0, 1.0);
    for(int j = 0; j < 3; j++)
      con_nn.A(tri_normals(it,j), tri_normals(it,j)) = 1.0;
    }

  return { tri_normals, tri_areas };
}

struct EdgeProperties
{
  // Edge length variables
  qcqp::VariableRef<1> edge_lengths;

  // Maps each (triangle,vertex) pair to an edge
  std::vector< std::tuple<int, int> > edge_triangle_index;

  // Maps each edge to a (triangle,vertex) pair
  vnl_matrix<int> tri_edges;
};

EdgeProperties
ComputeEdgeProperties(
    qcqp::Problem &p,
    TriangleMesh *mesh,
    qcqp::VariableRef<2> X,
    const std::string &name_prefix,
    double min_edge_length = 0.0)
{
  // Number of triangles
  int nt = mesh->triangles.size();

  // Array names
  std::string nm_el = name_prefix + "EL", nm_con_el = name_prefix + "EL";

  // First pass is to compute all edges and mark unique edges
  int n_edges = 0;
  vnl_matrix<int> tri_edges(nt, 3, -1);
  std::vector< std::tuple<int, int> > edge_tri_idx;
  for(int it = 0; it < nt; it++)
    {
    Triangle &t = mesh->triangles[it];
    for(int d = 0; d < 3; d++)
      {
      if(tri_edges(it, d) < 0)
        {
        tri_edges(it, d) = n_edges;
        if(t.neighbors[d] != NOID)
          tri_edges(t.neighbors[d],t.nedges[d]) = n_edges;
        n_edges++;
        edge_tri_idx.push_back({it, d});
        }
      }
    }

  // Create the edge lengths variable
  qcqp::VariableRef<1> edge_lengths = p.AddVariableTensor<1>(nm_el.c_str(), {n_edges}, 0.0);

  // Now we can compute the edge lengths and set constraints
  for(int ie = 0; ie < n_edges; ie++)
    {
    int it = std::get<0>(edge_tri_idx[ie]);
    int d = std::get<1>(edge_tri_idx[ie]);

    // The opposite vertices
    Triangle &t = mesh->triangles[it];
    int v1 = t.vertices[(d + 1) % 3], v2 = t.vertices[(d + 2) % 3];

    // Compute the edge length and assign as initial value
    double e_len = (X(v1).as_vector<3>() - X(v2).as_vector<3>()).magnitude();
    edge_lengths(ie) = e_len;

    // Create a constraint: square of the edge length is equal to the squared
    // distance between the points, el * el = sum_j (v2[j] - v1[j])^2
    auto &con_el = p.AddConstraint(nm_con_el.c_str(), 0.0, 0.0);
    con_el.A(edge_lengths(ie), edge_lengths(ie)) = 1.0;
    for(unsigned int j = 0; j < 3; j++)
      {
      con_el.A(X(v1,j), X(v1,j)) = -1.0;
      con_el.A(X(v2,j), X(v2,j)) = -1.0;
      con_el.A(X(v1,j), X(v2,j)) = 2.0;
      }
    }

  return { edge_lengths, edge_tri_idx, tri_edges };
}

// Create dihedral angle constraints
void CreateTriangleMinAngleConstraints(
    qcqp::Problem &qp, TriangleMesh *mesh,
    TriangleProperties &tp, EdgeProperties &ep,
    double min_tri_angle, std::string name_prefix)
{
  // Names
  std::string nm_ca = name_prefix + "csc_alpha";
  std::string nm_con_ca = name_prefix + "csc_alpha";

  int nt = (int) mesh->triangles.size();
  double min_angle = vnl_math::pi * min_tri_angle / 180;
  double max_csc = 1.0 / sin(min_angle);

  // Create the csc_alpha variable
  qcqp::VariableRef<2> q_csc_alpha =
      qp.AddVariableTensor<2>(nm_ca.c_str(), {nt, 3}, qcqp::LBINF, max_csc);

  for(int k = 0; k < nt; k++)
    {
    for(int d = 0; d < 3; d++)
      {
      // Index of edges around vertex d in triangle k
      int e1 = ep.tri_edges(k, (d + 1) % 3);
      int e2 = ep.tri_edges(k, (d + 2) % 3);

      qcqp::ElementRef l1 = ep.edge_lengths(e1);
      qcqp::ElementRef l2 = ep.edge_lengths(e2);
      qcqp::ElementRef area = tp.triangle_areas(k);
      qcqp::ElementRef csc_alpha = q_csc_alpha(k,d);

      // Compute the actual csc_alpha
      csc_alpha = l1 * l2 / (2. * area);

      // Create a constraint of the form 2*At*csc(alpha) = l1 * l2;
      auto &con_csc_alpha = qp.AddConstraint(nm_con_ca.c_str(), 0, 0);
      con_csc_alpha.A(area, csc_alpha) = 2.0;
      con_csc_alpha.A(l1, l2) = -1.0;
      }
    }
}


void CreateTriangleDihedralAngleConstraints(
    qcqp::Problem &qp, TriangleMesh *mesh, TriangleProperties &tp,
    vnl_matrix<double> &edgewise_min_angle, std::string name_prefix)
{
  std::string nm_con_da = name_prefix + "DA";
  for(size_t k = 0; k < mesh->triangles.size(); k++)
    {
    for(int d = 0; d < 3; d++)
      {
      size_t kopp = mesh->triangles[k].neighbors[d];
      double min_da = edgewise_min_angle(k,d);

      // One side only!
      if(kopp != NOID && k < kopp && min_da > 0.0)
        {
        // The dot product between triangle normals should be greater than this
        double min_cos_da = cos(vnl_math::pi * (180. - min_da) / 180.);

        // Create the constraint
        auto &con_cos_da = qp.AddConstraint(nm_con_da.c_str(), min_cos_da, qcqp::UBINF);
        for(int j = 0; j < 3; j++)
          {
          con_cos_da.A(tp.triangle_normals(k,j), tp.triangle_normals(kopp,j)) = 1.0;
          }
        }
      }
    }
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

  // Subdivision depth (optional info)
  std::vector<int> subDepth;

  void Load(const char *file);
  void Save(const char *file);

  void Subdivide(bool edge_only, bool flat_mode);
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

  vtkSmartPointer<vtkDataArray> sdepth = pd->GetPointData()->GetArray("SubdivisionDepth");

  // Get the coordinates and the medial index
  x.resize(pd->GetNumberOfPoints());
  mIndex.resize(pd->GetNumberOfPoints());
  subDepth.resize(pd->GetNumberOfPoints(), 0);
  for(int i = 0; i < pd->GetNumberOfPoints(); i++)
    {
    x[i].set(pd->GetPoint(i));
    mIndex[i] = (int) floor(0.5 + mix->GetTuple1(i));
    if(sdepth)
      subDepth[i] = (int) floor(0.5 + sdepth->GetTuple1(i));
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

  vtkSmartPointer<vtkIntArray> sdepth = vtkSmartPointer<vtkIntArray>::New();
  sdepth->SetNumberOfComponents(1);
  sdepth->Allocate(bmesh.nVertices);
  sdepth->SetName("SubdivisionDepth");

  for(int i = 0; i < bmesh.nVertices; i++)
    {
    pts->InsertNextPoint(x[i].data_block());
    mix->InsertNextTuple1(mIndex[i]);
    if(subDepth.size())
      sdepth->InsertNextTuple1(subDepth[i]);
    }

  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  pd->Allocate(bmesh.triangles.size());
  pd->SetPoints(pts);
  pd->GetPointData()->AddArray(mix);
  if(subDepth.size())
    pd->GetPointData()->AddArray(sdepth);
  
  if(Nx.size())
    {
    vtkSmartPointer<vtkFloatArray> nrm = vtkSmartPointer<vtkFloatArray>::New();
    nrm->SetNumberOfComponents(3);
    nrm->Allocate(bmesh.nVertices);

    for(int i = 0; i < bmesh.nVertices; i++)
      nrm->InsertNextTuple(Nx[i].data_block());

    pd->GetPointData()->SetNormals(nrm);
    }

  if(R.size())
    {
    vtkSmartPointer<vtkFloatArray> rad = vtkSmartPointer<vtkFloatArray>::New();
    rad->SetNumberOfComponents(1);
    rad->Allocate(bmesh.nVertices);
    rad->SetName("Radius");

    for(int i = 0; i < bmesh.nVertices; i++)
      rad->InsertNextTuple1(R[i]);

    pd->GetPointData()->AddArray(rad);
    }

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

template <class T>
void apply_subdivision(std::vector<T> &z, ImmutableSparseMatrix<double> &W)
{
  std::vector<T> zsrc = z;
  z.clear();

  for(int i = 0; i < W.GetNumberOfRows(); i++)
    {
    T zi(0.0);
    for(typename ImmutableSparseMatrix<double>::RowIterator it = W.Row(i); !it.IsAtEnd(); ++it)
      zi += zsrc[it.Column()] * it.Value();
    z.push_back(zi);
    }
}

void BCMTemplate::Subdivide(bool edges_only, bool flat_mode)
{
  // Make a copy of the boundary mesh
  bmesh.SetAsRoot();
  SubdivisionSurface::MeshLevel src = bmesh;
  
  // If subdividing edges only, mark those triangles
  if(edges_only)
    {
    // Triangle selection
    std::set<size_t> tsel;
    
    // Iterate over triangles and find the ones on the edge
    for(size_t t = 0; t < src.triangles.size(); t++)
      {
      auto &tri = src.triangles[t];
      for(int a = 0; a < 3; a++)
        if(tri.neighbors[a] == NOID)
          tsel.insert(t);
      }
    
    // Only subdivide those triangles
    // TODO: This code is broken. If you run this, the code below where the mIndex is assigned
    // will crash. This can be overcome but needs additional coding
    SubdivisionSurface::SubdivideSelected(&src, &bmesh, tsel);
    }
  else
    {
    // Subdivide into the current mesh
    SubdivisionSurface::Subdivide(&src, &bmesh, flat_mode);
    }

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
  apply_subdivision(x, W);
  if(Nx.size())
    apply_subdivision(Nx, W);
  if(R.size())
    apply_subdivision(R, W);

  // Compute the m-indices
  mIndex.resize(bmesh.nVertices, -1);
  subDepth.resize(bmesh.nVertices, -1);
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

      if(subDepth[vnew[j]] == -1)
        {
        int max_depth = std::max(subDepth[vold[(j+1) % 3]], subDepth[vold[(j+2) % 3]]);
        subDepth[vnew[j]] = max_depth + 1;
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

void merge_vertices(TriangleMesh &tm, vnl_matrix<double> &X, unsigned int v1, unsigned int v2)
{
  // Get rid of references to v2
  for(auto &tr : tm.triangles)
    {
    for(unsigned int i = 0; i < 3; i++)
      {
      if(tr.vertices[i] == v2)
        tr.vertices[i] = v1;
      }
    }

  // Merge the vertices themselves
  X.set_row(v1, (X.get_row(v1) + X.get_row(v2)) * 0.5);
}

template<class T>
unsigned int count_nnz(vnl_sparse_matrix<T> &mat)
{
  unsigned int nnz = 0;
  for(unsigned int i = 0; i < mat.rows(); i++)
    {
    auto &r = mat.get_row(i);
    for(unsigned int j = 0; j < r.size(); j++)
      {
      if(r[j].second != 0.0)
        nnz++;
      }
    }
  return nnz;
}

void InflateMedialModelWithBranches(const char *fn_input, const char *fn_output, double rad, int edge_label)
{
  // This inflation code accepts non-mesh medial surfaces, i.e., medial surfaces with branches. We can no
  // longer use the half-edge construct.

  // Load a triangular mesh from file
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(fn_input);
  reader->Update();

  // Convert it into a triangle mesh
  vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();
  vtkSmartPointer<vtkDataArray> alab = pd->GetPointData()->GetArray("Label");
  unsigned int nv = (unsigned int) pd->GetNumberOfPoints();

  // An edge is a pair of vertices, always stored in sorted order
  typedef std::pair<unsigned int, unsigned int> Edge;

  // A reference to a triangle edge (tri index, edge index, forward/backward)
  typedef std::tuple<unsigned int, unsigned int, bool> TriEdgeRef;

  // Each edge is associated with some number of triangles
  typedef std::map<Edge, std::list<TriEdgeRef> > EdgeTriMap;

  // Edge-triangle map
  EdgeTriMap etm;

  // List of duplicate triangles
  std::vector<Triangle> tdup;

  // List of triangle normals
  typedef vnl_vector_fixed<double, 3> Vec3;
  std::vector<Vec3> tnorm;

  // Find all the edges in the mesh
  for(unsigned int i = 0; i < pd->GetNumberOfCells(); i++)
    {
    // Read the cell
    vtkCell *c = pd->GetCell(i);
    if(c->GetNumberOfPoints() != 3)
      throw MedialModelException("Bad cell in input");

    // Create duplicates of the triangle with opposite windings
    Triangle t1;
    t1.vertices[0] = c->GetPointId(0);
    t1.vertices[1] = c->GetPointId(1);
    t1.vertices[2] = c->GetPointId(2);
    tdup.push_back(t1);

    Triangle t2;
    t2.vertices[0] = c->GetPointId(2);
    t2.vertices[1] = c->GetPointId(1);
    t2.vertices[2] = c->GetPointId(0);
    tdup.push_back(t2);

    // Compute the normals of the two triangles
    Vec3 A(pd->GetPoint(c->GetPointId(0)));
    Vec3 B(pd->GetPoint(c->GetPointId(1)));
    Vec3 C(pd->GetPoint(c->GetPointId(2)));
    Vec3 N = vnl_cross_3d(B-A, C-A).normalize();

    tnorm.push_back(N);
    tnorm.push_back(-N);
    }

  // Find edges across the duplicate triangles
  for(unsigned int i = 0; i < tdup.size(); i++)
    {
    for(unsigned int k = 0; k < 3; k++)
      {
      size_t v1 = tdup[i].vertices[(k+1) % 3];
      size_t v2 = tdup[i].vertices[(k+2) % 3];
      Edge ek = make_pair(min(v1,v2), max(v1,v2));

      // Add the triangle to this edge, marking it as either forward-traversed
      // or backward traversed.
      etm[ek].push_back(std::make_tuple(i, k, v1 > v2));
      }
    }

  // For each edge, find the triangles that are adjacent across the edge. Adjacent
  // triangles must traverse the edge in opposite order.
  for(auto &eit : etm)
    {
    // Get the edge vector direction
    Vec3 e_X1(pd->GetPoint(eit.first.first));
    Vec3 e_X2(pd->GetPoint(eit.first.second));

    for(auto &tref : eit.second)
      {
      // Get the normal of the current triangle
      unsigned int i_tri = std::get<0>(tref);
      unsigned int i_tri_edge_idx = std::get<1>(tref);
      bool winding = std::get<2>(tref);

      const Vec3 &N = tnorm[i_tri];
      Vec3 Z = (e_X2 - e_X1).normalize();
      if(!winding)
        Z = -Z;
      Vec3 X = vnl_cross_3d(Z, N);

      // Find the triangle that is closest by converting each opposite-winded triangle
      // to an angle and selecting the one with the minimum angle
      unsigned int opp_tri = NOID, opp_tri_edge_idx = -1;
      double min_angle = 0.0;
      for(auto &tref_test : eit.second)
        {
        // Only consider opposite winding
        if(std::get<2>(tref_test) != std::get<2>(tref))
          {
          // Find the 'X' of the test triangle
          unsigned int i_tri_test = std::get<0>(tref_test);
          const Vec3 &N_test = -tnorm[i_tri_test];
          Vec3 X_test = vnl_cross_3d(Z, N_test);

          // Find the angle with the current triangle.
          double a_test = (i_tri / 2 == i_tri_test / 2)
                          ? vnl_math::twopi
                          : atan2(dot_product(X_test, N), dot_product(X_test, X));
          if(a_test <= 0.0)
            a_test += vnl_math::twopi;

          printf("Angle of triangle %d with triangle %d over edge (%d,%d) is %f\n",
                 i_tri, i_tri_test, eit.first.first, eit.first.second, a_test);

          // Is this the best match
          if(opp_tri == NOID || a_test < min_angle)
            {
            opp_tri = i_tri_test;
            opp_tri_edge_idx = std::get<1>(tref_test);
            min_angle = a_test;
            }
          }
        }

      // We can now mark the neighbor of the triangle across this edge
      tdup[i_tri].neighbors[i_tri_edge_idx] = opp_tri;
      tdup[i_tri].nedges[i_tri_edge_idx] = (short) opp_tri_edge_idx;

      printf("Triangle %d matched to triangle %d\n", i_tri, opp_tri);
      }
    }

  // Create a vertex adjacency matrix. The rows/columns refer to the triangle vertices
  // which are at this point all considered to be disjoint points. When the adjacency
  // matrix contains 1, this means that the two vertices are actually the same point
  vnl_sparse_matrix<int> tv_adj(tdup.size() * 3, tdup.size() * 3);

  // Visit each edge in each triangle and match the vertices with the opposite edge
  // in the opposite triangle
  for(unsigned int i = 0; i < tdup.size(); i++)
    {
    for(unsigned int k = 0; k < 3; k++)
      {
      // Add identity element to matrix
      tv_adj(i * 3  + k, i * 3 + k) = 1;

      // Take triangle that's opposite
      unsigned int i_opp = tdup[i].neighbors[k];
      if(i_opp == NOID)
        throw MedialModelException("Triangle missing neighbor");

      // Set the matches
      unsigned int k_opp = tdup[i].nedges[k];
      unsigned int v1 = (k + 1) % 3, v2 = (k + 2) % 3;
      unsigned int v1_opp = (k_opp + 1) % 3, v2_opp = (k_opp + 2) % 3;

      tv_adj(i * 3 + v1, i_opp * 3 + v2_opp) = 1;
      tv_adj(i * 3 + v2, i_opp * 3 + v1_opp) = 1;
      }
    }

  // Find the connected components in the adjacency matrix. A lazy way to do this is to take powers of the
  // matrix until it converges.
  unsigned int nnz_last = count_nnz(tv_adj);
  printf("Adjacency matrix, nnz = %d\n", nnz_last);
  vnl_sparse_matrix<int> tv_adj_pow = tv_adj * tv_adj;
  while(count_nnz(tv_adj_pow) > nnz_last)
    {
    nnz_last = count_nnz(tv_adj_pow);
    tv_adj_pow = tv_adj_pow * tv_adj;
    printf("Adjacency multiplication, nnz = %d\n", nnz_last);
    }

  // Go through and remap the disjoint vertices to new vertices
  std::vector<unsigned int> vnew(tdup.size() * 3, NOID);
  unsigned int vcurr = 0;
  for(unsigned int i = 0; i < tdup.size() * 3; i++)
    {
    if(vnew[i] == NOID)
      {
      // Assign a new vertex ID to this vertex
      vnew[i] = vcurr;

      // Assign it to every other vertex in its row
      auto &row = tv_adj_pow.get_row(i);
      for(unsigned int j = 0; j < row.size(); j++)
        {
        if(vnew[row[j].first] != NOID && vnew[row[j].first] != vcurr)
          throw MedialModelException("Vertex traversal logic violation");

        vnew[row[j].first] = vcurr;
        }
      vcurr++;
      }
    }

  // Now we have a valid mesh structure in place. We can store this into a proper
  // triangle array
  vnl_matrix<unsigned int> m_tri(tdup.size(), 3);

  // We also need to compute the positions of the new vertices, i.e., by pushing them out
  // along the outward normals. We initialize each point to its original mesh location and
  // then add to the vertex all the normals of all the triangles that contain it
  vnl_matrix<double> m_pt(vcurr, 3), m_pt_offset(vcurr, 3);
  std::vector<unsigned int> valence(vcurr, 0);

  // Create the medial index array - this is just the original medial vertex
  vnl_vector<int> m_mindex(vcurr);

  // First pass through triangles, assigning new vertices and vertex coordinates
  for(unsigned int i = 0; i < tdup.size(); i++)
    {
    // Compute the triangle normal and center
    Vec3 P[] = {
      Vec3(pd->GetPoint(tdup[i].vertices[0])),
      Vec3(pd->GetPoint(tdup[i].vertices[1])),
      Vec3(pd->GetPoint(tdup[i].vertices[2])) };

    Vec3 N = vnl_cross_3d(P[1]-P[0], P[2]-P[0]).normalize();
    Vec3 C = (P[0] + P[1] + P[2])/3.0;

    for(unsigned int k = 0; k < 3; k++)
      {
      // Assign the new vertex
      m_tri(i,k) = vnew[i * 3 + k];

      // Get the coordinate of this vertex
      m_pt.set_row(m_tri(i,k), P[k]);

      // Add up the valence of this vertex
      valence[m_tri(i,k)]++;

      // Add up to the shift vector
      m_pt_offset.set_row(m_tri(i,k), m_pt_offset.get_row(m_tri(i,k)) + N);

      // Set the medial index (original index before inflation)
      m_mindex[m_tri(i,k)] = tdup[i].vertices[k];
      }
    }

  // Offset the vertices
  for(unsigned int j = 0; j < vcurr; j++)
    m_pt.set_row(j, m_pt.get_row(j) + rad * m_pt_offset.get_row(j) / valence[j]);


  // Generate an output triangle mesh
  VTKMeshBuilder<vtkPolyData> vmb;
  vmb.SetTriangles(m_tri);
  vmb.SetPoints(m_pt);
  vmb.AddArray(m_mindex, "MedialIndex");
  vmb.Save(fn_output);
}


void InflateMedialModel(const char *fn_input, const char *fn_output, double rad, int edge_label)
{
  typedef vnl_vector_fixed<double, 3> Vec3;
  typedef vnl_matrix_fixed<double, 3, 3> Mat3;

  // Load a triangular mesh from file
  vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(fn_input);
  reader->Update();

  // Convert it into a triangle mesh
  vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();
  vtkSmartPointer<vtkDataArray> alab = pd->GetPointData()->GetArray("Label");
  unsigned int nv = pd->GetNumberOfPoints();

  // On the first pass, we just duplicate each triangle and each vertex. Each vertex
  // i is matched to the corresponding vertex nv + i
  TriangleMesh model;
  TriangleMeshGenerator tmg(&model, 2 * nv);

  // Constants for pushing out vertices along the normals and laterally
  double r_norm = rad, r_lat = 0.5 * rad;

  // The coordinates of the vertices. These are constructed additively by pushing
  // each triangle vertex along its normal and outwards
  vnl_matrix<double> X(nv * 2, 3, 0.0);

  // The valence of the vertices - used for coordinate computation
  vnl_vector<unsigned int> valence(nv * 2, 0);

  // For each triangle in the mesh, generate another triangle that is its opposite. 
  for(int i = 0; i < pd->GetNumberOfCells(); i++)
    {
    // Read the cell
    vtkCell *c = pd->GetCell(i);
    if(c->GetNumberOfPoints() != 3)
      throw MedialModelException("Bad cell in input");

    // The vertices of the triangle
    vnl_vector_fixed<unsigned int, 3> v, w;
    Mat3 x;
    for(unsigned int a = 0; a < 3; a++)
      {
      v[a] = c->GetPointId(a);
      w[a] = v[a] + nv;
      for(unsigned int b = 0; b < 3; b++)
        x(a,b) = pd->GetPoint(v[a])[b];
      }

    // Add the two triangles
    tmg.AddTriangle(v[0], v[1], v[2]);
    tmg.AddTriangle(w[2], w[1], w[0]);

    // Get triangle normal and center
    Vec3 N = vnl_cross_3d(x.get_row(1)-x.get_row(0), x.get_row(2)-x.get_row(0)).normalize();
    Vec3 C = (x.get_row(0) + x.get_row(1) + x.get_row(2)) / 3.0;

    // Update the coordinates and valences
    for(unsigned int a = 0; a < 3; a++)
      {
      Vec3 y1 = x.get_row(a) * (1 + r_lat) - C * r_lat + N * r_norm;
      Vec3 y2 = x.get_row(a) * (1 + r_lat) - C * r_lat - N * r_norm;
      X.set_row(v[a], X.get_row(v[a]) + y1);
      X.set_row(w[a], X.get_row(w[a]) + y2);
      valence[v[a]]++; valence[w[a]]++;
      }
    }

  // Scale the X by valence
  for(int j = 0; j < 2 * nv; j++)
    X.set_row(j, X.get_row(j) * 1.0 / valence[j]);

  // Generate the mesh
  tmg.GenerateMesh();

  // Write the mesh
  /*
  VTKMeshBuilder<vtkPolyData> vmb;
  vmb.SetPoints(X);
  vmb.SetTriangles(model);
  vmb.Save("/tmp/lifted.vtk");
  */

  // Find all the edge vertices and merge them
  std::vector<int> remap(X.rows(), 0);
  for(unsigned int j = 0; j < nv; j++)
    {
    printf("label : %d\n", (alab ? (int) alab->GetTuple1(j) : -1));
    if(!model.IsVertexInternal(j) && (!alab || edge_label < 0 || (int) alab->GetTuple1(j) == edge_label))
      {
      merge_vertices(model, X, j, j + nv);
      remap[j+nv] = -1;
      }
    }

  // Create a vertex index mapping
  unsigned int k = 0;
  for(unsigned int j = 0; j < remap.size(); j++)
    {
    if(remap[j] == 0)
      remap[j] = k++;
    }

  // Apply the mapping
  TriangleMesh model_merged;
  TriangleMeshGenerator tmg_merged(&model_merged, k);
  for(auto &tr : model.triangles)
    tmg_merged.AddTriangle(remap[tr.vertices[0]], remap[tr.vertices[1]], remap[tr.vertices[2]]);
  tmg_merged.GenerateMesh();

  // Medial index
  vnl_vector<int> mindex(k, 0);

  // Apply to X
  vnl_matrix<double> X_merged(k, 3);
  for(unsigned int j = 0; j < remap.size(); j++)
    {
    if(remap[j] >= 0)
      {
      X_merged.set_row(remap[j], X.get_row(j));
      mindex[remap[j]] = j < nv ? j : j - nv;
      }
    }

  // Write the mesh
  VTKMeshBuilder<vtkPolyData> vmb_merged;
  vmb_merged.SetPoints(X_merged);
  vmb_merged.SetTriangles(model_merged);
  vmb_merged.AddArray(mindex,"MedialIndex");
  vmb_merged.Save(fn_output);
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
    bool edges_only, bool flat_mode,
    std::string fnOutput)
{
  // Load the template!
  BCMTemplate tmpl;
  tmpl.Load(fnTemplate.c_str());

  // Create the output template
  for(int i = 0; i < subdivisionLevel; i++)
    tmpl.Subdivide(edges_only, flat_mode);

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
  tParent.subDepth.resize(tParent.bmesh.nVertices);
  for(int i = 0; i < tParent.bmesh.nVertices; i++)
    {
    tParent.mIndex[i] = tChild.mIndex[i];
    tParent.subDepth[i] = tChild.subDepth[i];
    }


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
  ACTION_FIT_TARGET_OMT,
  ACTION_FIT_SELF,
  ACTION_SUBDIVIDE,
  ACTION_SUBDIVIDE_EDGE,
  ACTION_CONVERT_CMREP,
  ACTION_INFLATE_CMREP,
  ACTION_FIX_CMREP,
  ACTION_NONE
};

int usage()
{
  const char *usage =
      "bcmrep_main - Boundary-constrained medial representation utility\n"
      "usage:\n"
      "  bcmrep_main [options]\n"
      "options that select an action (use one):\n"
      "  -omt template.vtk target.vtk   : fit template to target with optimal mass transport\n"
      "  -icp template.vtk target.vtk   : fit template to target with iterative closest point loss\n"
      "  -lsq template.vtk target.vtk   : fit template to target least squares loss (matching vertices)\n"
      "  -sub template.vtk factor       : subdivide template by a factor\n"
      "  -cmr input.cmrep               : import a cm-rep template\n"
      "  -cfx input.cmrep               : fix a cm-rep template with bad triangles\n"
      "  -inf medial.vtk radius edgelab : inflate a medial model created with the GUI tool.\n"
      "                                   edgelab specifies labels of edge vertices, or -1\n"
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
      "  -max-iter N                    : maximum iterations for IpOpt (200)\n"
      "  -icp-iter N                    : Number of ICP iterations (5)\n"
      "  -no-fit-self                   : For ICP/OMT, skip the fit-to-self step\n"
      "triangle shape constraints:\n"
      "  -bc alpha beta min_area        : Set boundary triangle constraints \n"
      "                                     alpha:    minimum triangle angle (def: 12)\n"
      "                                     beta:     minimum dihedral angle (def: 120)\n"
      "                                     min_area: minimum triangle area  (def: 0.1)\n"
      "  -mc alpha beta min_area        : Set medial triangle constraints (see above) \n"
      "                                     by default all three are set to 0 (no constraint)\n"
      "subdivision-related options:\n"
      "  -sub-edge                      : only subdivide trianges on free edges and branches\n"
      "  -sub-flat                      : flat (as opposed to Loop) subdivision\n"
      "IPOpt options:\n"
      "  -hsllib <path>                 : path to the HSL dynamic library to load at runtime\n"
      "  -solver NAME                   : select which solver to use (see Coin-Or HSL docs. Def: ma86)\n"
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
  std::string Solver, HSLDylibPath, PardisoDylibPath;
  bool UseHessian;

  double EdgeLengthWeight;
  unsigned int MaxIpOptIter, MaxICPIter;

  bool SubdivisionFlat, SubdivisionEdgeOnly;
  bool SkipSelfFit;

  RegularizationOptions()
    : SubdivisionLevel(0), BasisSize(20),
      Weight(4.0), Mode(SPECTRAL), Solver("ma86"),
      EdgeLengthWeight(0.0), UseHessian(true), MaxIpOptIter(200), MaxICPIter(5),
      SubdivisionFlat(false), SubdivisionEdgeOnly(false), SkipSelfFit(false) {}
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

struct Statistics
{
  int n = 0;
  double sum = 0.0;
  double sum_sq = 0.0;
  double minx = qcqp::UBINF, maxx = qcqp::LBINF;

  void add(double x)
  {
    n++;
    sum += x; sum_sq += x*x;
    minx = std::min(minx, x);
    maxx = std::max(maxx, x);
  }

  double mean() { return sum / n; }
  double sd() { return sqrt(sum_sq / n - sum * sum / (n * n)); }
  double min() { return minx; }
  double max() { return maxx; }
};

std::string qp_var_name(qcqp::Problem &qp, int i)
{
  for(int q = 0; q < qp.GetNumberOfVariableTensors(); q++)
    {
    auto *vt = qp.GetVariableTensor(q);
    if(vt->GetStartIndex() <= i && vt->GetStartIndex()+vt->GetFlatSize() > i)
      return vt->GetName();
    }

  return "";
}

class BCMRepQuadraticProblemBuilder
{
public:
  typedef std::vector<std::vector<int> > MedialToBoundaryIndex;

  BCMRepQuadraticProblemBuilder(
      BCMTemplate &tmpl,
      RegularizationOptions &regOpts,
      TriangleConstraint tc_bnd, TriangleConstraint tc_med);

  ~BCMRepQuadraticProblemBuilder()
    {
    delete mmesh;
    }

  // Build the geometry without attachment terms
  void BuildGeometry(const std::vector<SMLVec3d> &x);

  // Build the regularization terms
  void BuildRegularization(const std::vector<SMLVec3d> &x = std::vector<SMLVec3d>());

  // Build the least squares objective
  void BuildLeastSquaresObjective(const std::vector<SMLVec3d> &xLSQTarget);

  // Build the least squares objective for triangle face centers
  void BuildFacesLeastSquaresObjective(const std::vector<SMLVec3d> &xLSQTarget);

  // Build an optimal mass transport objective
  void BuildOMTObjective(SinkhornIteration &si);

  // Build the ICP objective
  void BuildICPObjective(ClosestPointMatcher &cpm, double w_model_to_target, double w_target_to_model);

  // Save the current medial mesh
  void SaveBoundaryMesh(const char *fn);

  // Save the current medial mesh
  void SaveMedialMesh(const char *fn);

  // Get the wrapped problem
  qcqp::Problem &GetProblem() { return qp; }


public:
  // Parameters
  RegularizationOptions &regOpts;
  TriangleConstraint tc_bnd, tc_med;

  // Template
  BCMTemplate &tmpl;

  // Input meshes
  TriangleMesh *bmesh, *mmesh;

  // Medial link encoding
  std::vector<int> mIndex;
  MedialToBoundaryIndex mtbIndex;

  // Problem dimensions
  int nb, nm;

  // The wrapped problem
  qcqp::Problem qp;

  // The boundary positions and normals, medial positions
  qcqp::VariableRef<2> qX, qN, qM;

  // The radius values
  qcqp::VariableRef<1> qR;

  // Differential geometry at crests
  qcqp::VariableRef<3> qFF1, qSO;
  qcqp::VariableRef<1> qKappa1;

  // Triangle and edge geometry for medial and boundary triangles
  TriangleProperties tp_bnd, tp_med;
  EdgeProperties ep_bnd, ep_med;

  // Regularlization vars
  qcqp::VariableRef<2> qXC;
};

BCMRepQuadraticProblemBuilder::BCMRepQuadraticProblemBuilder(
    BCMTemplate &tmpl, RegularizationOptions &regOpts,
    TriangleConstraint tc_bnd, TriangleConstraint tc_med)
  : tmpl(tmpl), regOpts(regOpts), tc_bnd(tc_bnd), tc_med(tc_med)
{
  // Get the properties from the template
  bmesh = &tmpl.bmesh;
  mIndex = tmpl.mIndex;

  // Set the vertex numbers
  nb = bmesh->nVertices;

  // Create a medial mesh that duplicates the boundary mesh
  mmesh = new TriangleMesh(*bmesh);

  // Change the triangle vertices to use m_index
  for(int i = 0; i < mmesh->triangles.size(); i++)
    for(int j = 0; j < 3; j++)
      mmesh->triangles[i].vertices[j] = mIndex[bmesh->triangles[i].vertices[j]];

  // Get the number of medial points. TODO: we need to make sure that every
  // index between 0 and nm is represented.
  nm = 1 + *std::max_element(mIndex.begin(), mIndex.end());

  // Create a list of boundary atoms for each medial atom
  mtbIndex.resize(nm);
  for(int i = 0; i < nb; i++)
    mtbIndex[mIndex[i]].push_back(i);
}

void BCMRepQuadraticProblemBuilder::BuildGeometry(const std::vector<SMLVec3d> &x)
{
  // The boundary positions
  qX = qp.AddVariableTensor<2>("X", {nb, 3});

  // The medial positions
  qM = qp.AddVariableTensor<2>("M", {nm, 3});

  // The radius values
  qR = qp.AddVariableTensor<1>("R", {nm}, 0.);

  // The boundary normal vectors
  qN = qp.AddVariableTensor<2>("N", {nb, 3});

  // A list of all the boundary vertices that are crest points
  std::vector<int> crestVertices;

  // Index into the crest array for each boundary vertex
  std::vector<int> btcIndex(nb);

  // Assign the X variables
  for(int i = 0; i < nb; i++)
    for(int j = 0; j < 3; j++)
      qX(i,j) = x[i][j];

  // Compute crest indices
  for(int i = 0; i < nb; i++)
    {
    if(mtbIndex[mIndex[i]].size() == 1)
      {
      btcIndex[i] = crestVertices.size();
      crestVertices.push_back(i);
      }
    else
      {
      btcIndex[i] = -1;
      }
    }

  int nCrest = crestVertices.size();

  // Create the variables that are stored at the crest
  qFF1 = qp.AddVariableTensor<3>("FF1", {nCrest, 2, 2});
  qSO = qp.AddVariableTensor<3>("SO", {nCrest, 2, 2});
  qKappa1 = qp.AddVariableTensor<1>("Kappa1", {nCrest});

  // For each vertex, we need to temporarily store an array of neighbor weights
  // for partial derivative computation and the computed partial derivatives.
  struct VertexPartials
  {
    struct NeighborWeight { int index; double weight[2]; };
    std::vector<NeighborWeight> neighbor_weights;
    vnl_vector_fixed<double, 3> Xd[2];
  };

  std::vector<VertexPartials> vertex_partials(nb);

  // Create a LoopScheme for specifying normal vector constraints
  LoopTangentScheme lts;
  lts.SetMesh(bmesh);

  // The first pass computes the first partial derivatives and normals at all
  // boundary vertices, and sets the constraints linking N, M, R and X
  for(int i = 0; i < nb; i++)
    {
    // Medial atom index
    int iAtom = mIndex[i];

    // Walk around the vertex neighbors
    EdgeWalkAroundVertex walk(bmesh, i);

    // Collect the weights for this vertex for the u and v derivatives
    auto &vp_i = vertex_partials[i];

    // Take our own weight
    vp_i.neighbor_weights.push_back({ i, { lts.GetOwnWeight(0, i), lts.GetOwnWeight(1, i) } });

    // Take neighbor weights
    for(EdgeWalkAroundVertex walk(bmesh, i); !walk.IsAtEnd(); ++walk)
      {
      int i_nbr = (int) walk.MovingVertexId();
      vp_i.neighbor_weights.push_back(
            { i_nbr, { lts.GetNeighborWeight(0, walk), lts.GetNeighborWeight(1, walk) } });
      }

    // Compute the partial derivatives of X with respect to u/v and the
    // constraints of the normal vector
    for(int d = 0; d < 2; d++)
      {
      // Constraint that dot(N, Xd) == 0
      auto &con_norm_xd = qp.AddConstraint("N.Xu", 0.0, 0.0);
      vp_i.Xd[d].fill(0.0);
      for(auto w : vp_i.neighbor_weights)
        {
        if(w.weight[d] != 0.0)
          {
          // Add to Xu/Xv
          vp_i.Xd[d] += w.weight[d] * qX(w.index).as_vector<3u>();

          // Set coefficient for the constraint
          for(unsigned int j = 0; j < 3; j++)
            con_norm_xd.A(qN(i,j), qX(w.index,j)) = w.weight[d];
          }
        }
      }

    // Assign the initial value for the normal
    vnl_vector_fixed<double, 3> q_v_n = vnl_cross_3d(vp_i.Xd[0], vp_i.Xd[1]).normalize();
    qN(i).from_vector(q_v_n);

    // Constraint that dot(N, N) == 1
    auto &con_norm_len = qp.AddConstraint("N.N", 1.0, 1.0);
    for(unsigned int j = 0; j < 3; j++)
      con_norm_len.A(qN(i,j), qN(i,j)) = 1.0;

    // Now setup the medial constraints, X - r * N - M = 0
    for(int j = 0; j < 3; j++)
      {
      auto &con_medial = qp.AddConstraint("X-rNM", 0.0, 0.0);
      con_medial.b(qX(i,j)) = 1.0;
      con_medial.b(qM(iAtom,j)) = -1.0;
      con_medial.A(qN(i,j), qR(iAtom)) = -1.0;
      }
    }

  // The second pass is for second derivative properties, and requires all
  // the normals to have been computed already
  for(int i = 0; i < nb; i++)
    {
    int iAtom = mIndex[i];
    auto &vp_i = vertex_partials[i];

    // The rest of the constraints are in the single-tangency case
    int i_crest = btcIndex[i];
    if(i_crest >= 0)
      {
      // Compute the partial derivatives of N with respect to u/v and the
      // constraints of the first fundamental form
      vnl_vector_fixed<double, 3> q_Ndi[2];
      vnl_matrix_fixed<double, 2, 2> m_FF1, m_FF2, m_SO;
      for(int d = 0; d < 2; d++)
        {
        q_Ndi[d].fill(0.0);
        for(auto w : vp_i.neighbor_weights)
          {
          if(w.weight[d] != 0.0)
            {
            // Add to Nu/Nv
            q_Ndi[d] += w.weight[d] * qN(w.index).as_vector<3u>();
            }
          }

        // Set up first fundamental form constraints
        for(int r = 0; r < 2; r++)
          {
          // Set the value of the first fundamental form
          qFF1(i_crest, d, r) = m_FF1(d,r) = dot_product(vp_i.Xd[d], vp_i.Xd[r]);

          // Constraint that FF1[d][r] == dot(Xd[d], Xd[r])
          auto &con_ff1_dr = qp.AddConstraint("FF1", 0.0, 0.0);
          con_ff1_dr.b(qFF1(i_crest, d, r)) = 1;

          for(auto w1 : vp_i.neighbor_weights)
            {
            for(auto w2 : vp_i.neighbor_weights)
              {
              if(w1.weight[d] * w2.weight[r] != 0.0)
                {
                for(unsigned int j = 0; j < 3; j++)
                  con_ff1_dr.A(qX(w1.index, j), qX(w2.index, j)) -= w1.weight[d] * w2.weight[r];
                }
              }
            }
          }
        }

      // Compute the second fundamental form and shape operator
      for(int d = 0; d < 2; d++)
        for(int r = 0; r < 2; r++)
          m_FF2[d][r] = dot_product(vp_i.Xd[d], q_Ndi[r]);

      // Numerically solve for the shape operator
      m_SO = - vnl_inverse(m_FF1) * m_FF2;

      // Set up the constraints for the shape operator
      for(int d = 0; d < 2; d++)
        {
        for(int r = 0; r < 2; r++)
          {
          // Set the initial value for the shape operator variable
          qSO(i_crest, d, r) = m_SO(d,r);

          // Create the constraint on the shape operator, that FF1 * SO = SFF
          auto &con_so_dr = qp.AddConstraint("SO", 0.0, 0.0);

          // The left hand side of the constraint is (FF1 * SO)[d,r], which is to say
          // FF1[d][0] * SO[0][r] + FF1[d][1] * SO[1][r]
          con_so_dr.A(qFF1(i_crest, d, 0), qSO(i_crest, 0, r)) = 1.0;
          con_so_dr.A(qFF1(i_crest, d, 1), qSO(i_crest, 1, r)) = 1.0;

          // The right hand side is the second fundamental form, i.e., dot product of
          // the derivative of X and derivative of N
          for(auto w1 : vp_i.neighbor_weights)
            {
            for(auto w2 : vp_i.neighbor_weights)
              {
              if(w1.weight[d] * w2.weight[r] != 0.0)
                {
                for(unsigned int j = 0; j < 3; j++)
                  con_so_dr.A(qX(w1.index, j), qN(w2.index, j)) = w1.weight[d] * w2.weight[r];
                }
              }
            }
          }
        }

      // Numerically solve for the principal curvature k1
      double mH = vnl_trace(m_SO) / 2;
      double mK = vnl_det(m_SO);
      double mk1 = mH - sqrt(mH*mH - mK);

      // Assign the initial value for the principal curvature
      qKappa1(i_crest) = mk1;

      // The equality constraint on kappa1 is to ensure it is an eigenvalue, i.e,
      // det(SO-k1*I) = 0. Written out, this has the form
      // SO[0][0] * SO[1][1] - SO[0][1] * SO[1][0] - k1 * SO[0][0] - k1 * SO[1][1] + k1^2 = 0
      auto &con_kappa1 = qp.AddConstraint("Kappa1_eq", 0.0, 0.0);
      con_kappa1.A(qSO(i_crest,0,0), qSO(i_crest,1,1)) = 1.0;
      con_kappa1.A(qSO(i_crest,0,1), qSO(i_crest,1,0)) = -1.0;
      con_kappa1.A(qSO(i_crest,0,0), qKappa1(i_crest)) = -1.0;
      con_kappa1.A(qSO(i_crest,1,1), qKappa1(i_crest)) = -1.0;
      con_kappa1.A(qKappa1(i_crest), qKappa1(i_crest)) = 1.0;

      // The inequality constraint on kappa1 is to ensure it is less than the mean curvature
      // (trace of SO/2), i.e., k1 - 0.5 (S[0][0]+S[1][1]) < 0
      auto &con_kappa1_ineq = qp.AddConstraint("Kappa1_ineq", qcqp::LBINF, 0.0);
      con_kappa1_ineq.b(qKappa1(i_crest)) = 1.0;
      con_kappa1_ineq.b(qSO(i_crest,0,0)) = -0.5;
      con_kappa1_ineq.b(qSO(i_crest,1,1)) = -0.5;

      // Set the value of the radius at this crest atom
      qR(iAtom) = -1.0 / mk1;

      // Create constraint linking R and kappa1, R * kappa1 = -1
      auto &con_R_kappa1 = qp.AddConstraint("R_Kappa1", -1.0, -1.0);
      con_R_kappa1.A(qR(iAtom), qKappa1(i_crest)) = 1.0;

      // Set the value of the M variable based on the single R
      vnl_vector_fixed<double, 3> v_m = qX(i).as_vector<3>() - qN(i).as_vector<3>() * (double) qR(iAtom);
      qM(iAtom).from_vector(v_m);
      }
    }

  // Set initial values of M and R for the bitangent or greater medial atoms
  for(int i = 0; i < nm; i++)
    {
    int k = mtbIndex[i].size();
    if(k > 1)
      {
      // MULTIPLE TANGENCY CASE

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
        SMLVec3d Xq = qX(iq).as_vector<3>(), Nq = qN(iq).as_vector<3>();
        for(int p = 0; p < k; p++)
          {
          int ip = mtbIndex[i][p];
          SMLVec3d Xp = qX(ip).as_vector<3>(), Np = qN(ip).as_vector<3>();

          numerator += dot_product(Xq - Xp, Nq);
          denominator -= dot_product(Np, Nq);
          }

        sumX += Xq; sumN += Nq;
        }

      // Compute the best fit r and m
      double v_r = numerator / denominator;
      SMLVec3d v_m = (sumX - sumN * v_r) / ((double) k);

      // Store those values
      qR(i) = v_r;
      qM(i).from_vector(v_m);
      }
    }

  // ------------------------------------------------------------------------
  // Create the QCQP MIB constraint (NEW Code)
  // ------------------------------------------------------------------------
  for(int iBnd = 0; iBnd < nb; iBnd++)
    {
    int iAtom = mIndex[iBnd];

    for(EdgeWalkAroundVertex walk(bmesh, iBnd); !walk.IsAtEnd(); ++walk)
      {
      int k = walk.MovingVertexId();

      // Constraint is in the form |X[k] - M|^2 - R^2 > 0
      auto &con_mib = qp.AddConstraint("MIB", 0, qcqp::UBINF);
      for(int j = 0; j < 3; j++)
        {
        con_mib.A(qX(k,j), qX(k,j)) = 1.0;
        con_mib.A(qX(k,j), qM(iAtom,j)) = -2.0;
        con_mib.A(qM(iAtom,j), qM(iAtom,j)) = 1.0;
        con_mib.A(qR(iAtom), qR(iAtom)) = -1.0;
        }
      }
    }


  // ------------------------------------------------------------------------
  // Boundary QCQP Triangle and Edge Geometry Constraints (NEW CODE)
  // ------------------------------------------------------------------------
  if(tc_bnd.min_dih_angle > 0 || tc_bnd.min_tri_angle > 0 || tc_bnd.min_tri_area > 0)
    tp_bnd = ComputeTriangleProperties(qp, bmesh, qX, "bnd_", tc_bnd.min_tri_area);

  if(tc_bnd.min_tri_angle > 0 || regOpts.EdgeLengthWeight > 0.0)
    ep_bnd = ComputeEdgeProperties(qp, bmesh, qX, "bnd_", 0.0);

  // Boundary minimal angle constraint
  if(tc_bnd.min_tri_angle > 0)
    CreateTriangleMinAngleConstraints(
          qp, bmesh, tp_bnd, ep_bnd, tc_bnd.min_tri_angle, "bnd_");

  // Boundary minimal dihedral angle constraint
  if(tc_bnd.min_dih_angle > 0)
    {
    // First we need to find the minimal angle for each edge, since crest edges
    // should be allowed a larger dihedral angle than plain edges
    vnl_matrix<double> edgewise_min_angle(bmesh->triangles.size(), 3);
    for(unsigned int k = 0; k < bmesh->triangles.size(); k++)
      {
      for(int d = 0; d < 3; d++)
        {
        int v1 = bmesh->triangles[k].vertices[(d+1)%3];
        int v2 = bmesh->triangles[k].vertices[(d+2)%3];
        bool isBnd = (mtbIndex[mIndex[v1]].size() == 1) && (mtbIndex[mIndex[v2]].size() == 1);

        // TODO: why this hard-coded value?
        edgewise_min_angle(k, d) = isBnd ? 25.84 : tc_bnd.min_dih_angle;
        }
      }

    // Create the actual constraints
    CreateTriangleDihedralAngleConstraints(
          qp, bmesh, tp_bnd, edgewise_min_angle, "bnd_");
    }


  // ------------------------------------------------------------------------
  // Medial QCQP Triangle and Edge Geometry Constraints (NEW CODE)
  // ------------------------------------------------------------------------
  // TODO: this is redundant, as it computes two sets of equal normals and
  // areas for each medial triangle, one facing each boundary triangle. We
  // could avoid this by checking for opposite triangles.
  if(tc_med.min_dih_angle > 0 || tc_med.min_tri_angle > 0 || tc_med.min_tri_area > 0)
    tp_med = ComputeTriangleProperties(qp, mmesh, qX, "med_", tc_med.min_tri_area);

  if(tc_med.min_tri_angle > 0)
    ep_med = ComputeEdgeProperties(qp, mmesh, qX, "med_", 0.0);

  // Boundary dihedral angle constraint
  if(tc_med.min_tri_angle > 0)
    CreateTriangleMinAngleConstraints(
          qp, mmesh, tp_med, ep_med, tc_med.min_tri_angle, "med_");

  if(tc_med.min_dih_angle > 0)
    {
    // We first need to compute the multiplicity of each edge in the medial mesh so that we
    // can exclude the edges that are on seam and edge curves of the medial model. Those edges
    // are shared by a number of triangles other than two.
    typedef std::pair<size_t, size_t> Edge;
    std::map<Edge, unsigned int> edge_mult;

    // First pass computes the edge multiplicities, visiting each edge twice
    for(unsigned int k = 0; k < mmesh->triangles.size(); k++)
      {
      for(int d = 0; d < 3; d++)
        {
        unsigned int kopp = mmesh->triangles[k].neighbors[d];
        if(kopp != NOID && k < kopp)
          {
          // Is this a boundary edge?
          int v1 = mmesh->triangles[k].vertices[(d+1)%3];
          int v2 = mmesh->triangles[k].vertices[(d+2)%3];
          Edge edge(std::min(v1,v2), std::max(v1,v2));
          auto it = edge_mult.find(edge);
          if(it == edge_mult.end())
            edge_mult.insert(make_pair(edge, 1));
          else
            it->second++;
          }
        }
      }

    // Minimum dihedral angle for each edge (represented by triangle/vertex)
    vnl_matrix<double> edgewise_min_angle(bmesh->triangles.size(), 3);

    // Second pass assigns minimum DA to the edges
    for(unsigned int k = 0; k < mmesh->triangles.size(); k++)
      {
      for(int d = 0; d < 3; d++)
        {
        unsigned int kopp = mmesh->triangles[k].neighbors[d];
        if(kopp != NOID && k < kopp)
          {
          // Is this a boundary edge?
          int v1 = mmesh->triangles[k].vertices[(d+1)%3];
          int v2 = mmesh->triangles[k].vertices[(d+2)%3];
          Edge edge(std::min(v1,v2), std::max(v1,v2));
          edgewise_min_angle(k,d) = (edge_mult[edge] == 2) ? tc_med.min_dih_angle : 0;
          }
        }
      }

    // Finally, create the constraints
    CreateTriangleDihedralAngleConstraints(
          qp, mmesh, tp_med, edgewise_min_angle, "med_");
    }
}

void BCMRepQuadraticProblemBuilder::BuildRegularization(const std::vector<SMLVec3d> &x)
{
  // ------------------------------------------------------------------------
  // Define the QCQP regularization objective (NEW)
  // ------------------------------------------------------------------------
  auto &q_loss_basis_residual = qp.AddLoss("BasisResidual", 1.0);

  if(regOpts.Mode == RegularizationOptions::SUBDIVISION)
    {
    // Try to obtain a parent mesh
    BCMTemplate tmplParent;
    if(!ReverseEngineerSubdivisionMesh(tmpl, tmplParent, regOpts.SubdivisionLevel))
      throw MedialModelException("Unable to deduce coarse-level mesh by reversing subdivision.");

    // Save the parent-level mesh
    // tmplParent.Save("reverse_eng.vtk");

    // TODO - do we want to regularize R as well???

    // Create control points, i.e., vertices in the parent mesh
    SubdivisionSurface::MeshLevel &pmesh = tmplParent.bmesh;
    qXC = qp.AddVariableTensor<2>("XC", { (int) pmesh.nVertices, 3});

    // Assign the control point initial values
    for(int i = 0; i <  pmesh.nVertices; i++)
      qXC(i).from_vector<3>(x.size() ? x[i] : tmplParent.x[i]);

    // Apply the weight matrix to compute the residual
    ImmutableSparseMatrix<double> &W = tmpl.bmesh.weights;

    // Define the residual objective
    for(int i = 0; i < nb; i++)
      {
      for(int j = 0; j < 3; j++)
        {
        // A linear expression for the distance between interpolated and actual vertices
        qcqp::LinearExpression delta;
        delta.b(qX(i,j)) = -1.0;

        // Use the weighted matrix iterator
        for(ImmutableSparseMatrix<double>::RowIterator it = W.Row(i); !it.IsAtEnd(); ++it)
          delta.b(qXC(it.Column(), j)) += it.Value();

        // Add the squared distance to the loss
        q_loss_basis_residual.add_dotted(delta, delta);
        }
      }
    }

  if(regOpts.Mode == RegularizationOptions::SPECTRAL)
    {
    // Define a basis for the surface
    MeshBasisCoefficientMapping basismap_X(bmesh, regOpts.BasisSize, 3);

    // Define a basis for the medial axis
    // TODO: to what extent is this a valid basis - we need to visualize!
    // MeshBasisCoefficientMapping basismap_M(bmesh, nBasis, 4);

    // Create the coefficient variables
    qXC = qp.AddVariableTensor<2>("XC", { (int) regOpts.BasisSize, 3});
    for(int i = 0; i < regOpts.BasisSize; i++)
      for(int j = 0; j < 3; j++)
        qXC(i,j) = x.size() ? x[i][j] : 0.0;

    // Define the objective on the basis
    for(int iBnd = 0; iBnd < nb; iBnd++)
      {
      SMLVec3d Xfixed = tmpl.x[iBnd];

      for(int j = 0; j < 3; j++)
        {
        // Compute expression (X_init - X) - interpolate(XC)
        qcqp::LinearExpression delta;
        delta.c() = Xfixed[j];
        delta.b(qX(iBnd,j)) = -1.0;
        for(int i = 0; i < regOpts.BasisSize; i++)
          delta.b(qXC(i,j)) += basismap_X.GetBasisComponent(i, iBnd);

        // Add square of this expression to the loss
        q_loss_basis_residual.add_dotted(delta, delta);
        }
      }
    }

  // ------------------------------------------------------------------------
  // QCQP Boundary Total Edge Length Loss (NEW CODE)
  // ------------------------------------------------------------------------
  // TODO: add a constraint on the angle between adjacent boundary edges, b/c
  // that seems to become quite small sometimes and setting a minimum on it might
  // work better than penalizing the total length
  if(regOpts.EdgeLengthWeight > 0.0)
    {
    auto &loss_el = qp.AddLoss("edge_length", regOpts.EdgeLengthWeight);

    // Add all the crest edges
    for(unsigned int k = 0; k < bmesh->triangles.size(); k++)
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
            int i_edge = ep_bnd.tri_edges(k,d);
            loss_el.b(ep_bnd.edge_lengths(i_edge)) = 1.0;
            }
          }
        }
      }
    }
}

void BCMRepQuadraticProblemBuilder::BuildLeastSquaresObjective(const std::vector<SMLVec3d> &xLSQTarget)
{
  // ------------------------------------------------------------------------
  // QCQP displacement objective (NEW CODE)
  // ------------------------------------------------------------------------
  auto &loss_disp = qp.AddLoss("disp", 1.0);
  for(int i = 0; i < nb; i++)
    {
    for(unsigned int j = 0; j < 3; j++)
      {
      // Loss is in the form sum [ |X_i - xtarg_i|^2 ]
      double xtarg = xLSQTarget[i][j];
      loss_disp.A(qX(i,j), qX(i,j)) = 1.0;
      loss_disp.b(qX(i,j)) = -2.0 * xtarg;
      loss_disp.c() += xtarg * xtarg;
      }
    }
  }

void BCMRepQuadraticProblemBuilder::BuildFacesLeastSquaresObjective(const std::vector<SMLVec3d> &xFacesLSQTarget)
{
  // ------------------------------------------------------------------------
  // QCQP displacement objective (NEW CODE)
  // ------------------------------------------------------------------------
  auto &loss_disp = qp.AddLoss("disp_faces", 1.0);
  for(int i = 0; i < bmesh->triangles.size(); i++)
    {
    auto &t = bmesh->triangles[i];

    for(unsigned int j = 0; j < 3; j++)
      {
      // Loss is in the form sum [ |X_i - xtarg_i|^2 ]
      qcqp::ElementRef a = qX(t.vertices[0],j);
      qcqp::ElementRef b = qX(t.vertices[1],j);
      qcqp::ElementRef c = qX(t.vertices[2],j);
      double z = xFacesLSQTarget[i][j];

      /*
      qcqp::LinearExpression delta;
      delta.b(qX(t.vertices[0],j)) = 1.0 / 3.0;
      delta.b(qX(t.vertices[1],j)) = 1.0 / 3.0;
      delta.b(qX(t.vertices[2],j)) = 1.0 / 3.0;
      loss_disp.c() = -xtarg;
      loss_disp.add_dotted(delta, delta);
      */
      loss_disp.A(a,a) += 1. / 9.;
      loss_disp.A(b,b) += 1. / 9.;
      loss_disp.A(c,c) += 1. / 9.;
      loss_disp.A(a,b) += 2. / 9.;
      loss_disp.A(b,c) += 2. / 9.;
      loss_disp.A(c,a) += 2. / 9.;
      loss_disp.b(a) += -2. * z / 3.;
      loss_disp.b(b) += -2. * z / 3.;
      loss_disp.b(c) += -2. * z / 3.;
      loss_disp.c() += z * z;
      }
    }
}

void BCMRepQuadraticProblemBuilder::BuildOMTObjective(SinkhornIteration &si)
{
  int n0 = bmesh->triangles.size();
  int n1 = si.m_U1.size();

  // The objective is of the form |X_i - Y_j|^2 * (1 + alpha - alpha <N_i,M_j>)
  // where X,N are the model triangle centers and normals, and Y,M are the target
  // triangle centers and normals. To make this quadratic requires the introduction
  // of a new variable X_i * N_i
  auto &loss = qp.AddLoss("omp_loss", 1.0);

  // Add the variable representing the center of the triangle times the normal of
  // the triangle. We could also add the variable for the vertices instead, but I
  // am not sure that will make any difference in problem size or efficiency
  qcqp::VariableRef<3> qTXN = qp.AddVariableTensor<3>("TXN", {n0, 3, 3});

  // Iterate over the model mesh triangles
  for(int i = 0; i < n0; i++)
    {
    auto &t = bmesh->triangles[i];
    int a = t.vertices[0], b = t.vertices[1], c = t.vertices[2];

    // Initialize qTXN variable and set its constraint
    for(unsigned int d = 0; d < 3; d++)
      {
      for(unsigned int q = 0; q < 3; q++)
        {
        // Initialize the qTXN variable
        qTXN(i,d,q) = tp_bnd.triangle_normals(i,q) * (qX(a,d) + qX(b,d) + qX(c,d)) / 3.0;

        // Create a constraint for the qTXN variable
        auto &con_txn = qp.AddConstraint("TXN", 0.0, 0.0);
        con_txn.A(tp_bnd.triangle_normals(i,q), qX(a,d)) = 1./3.;
        con_txn.A(tp_bnd.triangle_normals(i,q), qX(b,d)) = 1./3.;
        con_txn.A(tp_bnd.triangle_normals(i,q), qX(c,d)) = 1./3.;
        con_txn.b(qTXN(i,d,q)) = -1.0;
        }
      }
    }

  // TODO: this is just for checking the computation
  auto x_curr = qp.GetVariableValues();

  // Iterate over the model mesh triangles
  for(int i = 0; i < n0; i++)
    {
    auto &t = bmesh->triangles[i];
    int a = t.vertices[0], b = t.vertices[1], c = t.vertices[2];

    // Create a temporary quadratic expression for this triangle
    qcqp::QuadraticExpression loss_i;

    // Expected value of the loss - obtained by adding up C * gamma in sinkhorn
    double loss_exp_value = 0.0;
    double loss_exp_value_2 = 0.0;

    // Running sums of coefficients to compute
    double w_1 = 0.0;
    vnl_vector_fixed<double, 3> w_Xd(0.0), w_XXd(0.0), w_Nq(0.0);
    vnl_matrix_fixed<double, 3, 3> w_XXd_Nq(0.0), w_Xd_Nq(0.0);

    // Compute these sums going over the target triangles
    double mass_i = 0.0;
    for(unsigned int j = 0; j < n1; j++)
      {
      // Get the transport probability between i and j
      double gamma_ij = si.m_U0(i) * si.m_K(i, j) * si.m_U1(j);
      mass_i += gamma_ij;

      // Get the target triangle center in native coordinates
      SMLVec3d Yj = si.m_ZTrgCenters[j] * si.m_TrgScale + si.m_TrgCenter;
      SMLVec3d Mj = si.m_ZTrgNormals[j];

      // Add up the weights
      for(unsigned int d = 0; d < 3; d++)
        {
        w_1       += Yj[d] * Yj[d] * (1+si.m_Alpha) * gamma_ij;
        w_XXd[d]  += (1+si.m_Alpha) * gamma_ij;
        w_Xd[d]   += -2 * Yj[d] * (1+si.m_Alpha) * gamma_ij;
        for(unsigned int q = 0; q < 3; q++)
          {
          w_XXd_Nq(d,q) += -si.m_Alpha * Mj[q] * gamma_ij;
          w_Xd_Nq(d,q)  += 2 * si.m_Alpha * Yj[d] * Mj[q] * gamma_ij;
          w_Nq[q]       += -si.m_Alpha * Yj[d] * Yj[d] * Mj[q] * gamma_ij;
          }
        }

      // Add the expected value
      loss_exp_value += gamma_ij * si.m_C(i,j) * si.m_TrgScale * si.m_TrgScale;

      for(unsigned int d = 0; d < 3; d++)
        {
        double xi = (qX(a,d) + qX(b,d) + qX(c,d))/3;
        double ca = dot_product(tp_bnd.triangle_normals(i).as_vector<3>(), Mj);
        loss_exp_value_2 += gamma_ij * (xi - Yj[d]) * (xi - Yj[d]) * (1 + si.m_Alpha - si.m_Alpha * ca);
        }
      }

    // Time to assign weights to all the variables. However, since the weights for the X's above
    // refer to the triangle centers, each weight has to be distributed among the vertices
    const double w19 = 1.0/9.0, w29 = 2.0/9.0, w13 = 1.0/3.0;
    for(unsigned int d = 0; d < 3; d++)
      {
      // All the elements involved
      auto Xa = qX(a,d), Xb = qX(b,d), Xc = qX(c,d), N = tp_bnd.triangle_normals(i,d);

      // Assign the weight for X^2 term
      loss_i.A(Xa, Xa) += w19 * w_XXd[d];
      loss_i.A(Xb, Xb) += w19 * w_XXd[d];
      loss_i.A(Xc, Xc) += w19 * w_XXd[d];
      loss_i.A(Xa, Xb) += w29 * w_XXd[d];
      loss_i.A(Xb, Xc) += w29 * w_XXd[d];
      loss_i.A(Xc, Xa) += w29 * w_XXd[d];

      // Assign the weight for X term
      loss_i.b(Xa) += w13 * w_Xd[d];
      loss_i.b(Xb) += w13 * w_Xd[d];
      loss_i.b(Xc) += w13 * w_Xd[d];

      for(unsigned int q = 0; q < 3; q++)
        {
        auto  XN = qTXN(i,d,q);

        // Assign weights for the term X * N
        loss_i.A(Xa, N) += w13 * w_Xd_Nq(d,q);
        loss_i.A(Xb, N) += w13 * w_Xd_Nq(d,q);
        loss_i.A(Xc, N) += w13 * w_Xd_Nq(d,q);

        // Assign weights for the term X * X * N
        loss_i.A(Xa, XN) += w13 * w_XXd_Nq(d,q);
        loss_i.A(Xb, XN) += w13 * w_XXd_Nq(d,q);
        loss_i.A(Xc, XN) += w13 * w_XXd_Nq(d,q);
        }

      // Assign weights for the term N (d here refers to q)
      loss_i.b(N) += w_Nq[d];
      }

    // Add constant term
    loss_i.c() += w_1;

    // Compare the loss with expected value
    double loss_value = loss_i.Evaluate(x_curr.data_block());
    printf("Triangle %5d:  OMT Loss expected: %12.8f  %12.8f    actual: %12.8f\n", i, loss_exp_value, loss_exp_value_2, loss_value);

    // Add the loss
    loss.add(loss_i);
    }
}

void BCMRepQuadraticProblemBuilder::BuildICPObjective(
    ClosestPointMatcher &cpm, double w_model_to_target, double w_target_to_model)
{
  // We need a vector with the X locations for later
  std::vector<SMLVec3d> x(nb);

  // Create the distance from model to target loss
  auto &loss_to_target = qp.AddLoss("model_to_target", w_model_to_target);
  for(int i = 0; i < nb; i++)
    {
    x[i] = qX(i).as_vector<3>();
    SMLVec3d xtarg = cpm.FindClosestToTarget(x[i]);
    for(unsigned int j = 0; j < 3; j++)
      {
      // Loss is in the form sum [ |X_i - xtarg_i|^2 ]
      loss_to_target.A(qX(i,j), qX(i,j)) = 1.0;
      loss_to_target.b(qX(i,j)) = -2.0 * xtarg[j];
      loss_to_target.c() += xtarg[j] * xtarg[j];
      }
    }

  // Create the distance from target to model loss
  auto &loss_to_model = qp.AddLoss("target_to_model", w_target_to_model);

  // Compute the match points
  auto x_match = cpm.FindClosestToSource(bmesh, x);

  for(int i = 0; i < x_match.size(); i++)
    {
    // Get the match triangle and vertices
    auto &loc = x_match[i];
    size_t *v = bmesh->triangles[loc.iTriangle].vertices;

    // Iterate over dimension
    for(int d = 0; d < 3; d++)
      {
      // Use linear expression to describe the distance long dimension d
      qcqp::LinearExpression delta;

      // Use the barycentric coordinates to express match location as weighted sum
      for(unsigned int j = 0; j < 3; j++)
        delta.b(qX(v[j], d)) = loc.xBary[j];
      delta.c() = -loc.xTarget[d];

      // Add quadratic expression
      loss_to_model.add_dotted(delta, delta);
      }
    }
}

void BCMRepQuadraticProblemBuilder::SaveBoundaryMesh(const char *fn)
{
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  pts->Allocate(nb);

  vtkSmartPointer<vtkFloatArray> rad = vtkSmartPointer<vtkFloatArray>::New();
  rad->SetNumberOfComponents(1);
  rad->Allocate(nb);
  rad->SetName("Radius");

  vtkSmartPointer<vtkIntArray> mix = vtkSmartPointer<vtkIntArray>::New();
  mix->SetNumberOfComponents(1);
  mix->Allocate(nb);
  mix->SetName("MedialIndex");

  vtkSmartPointer<vtkIntArray> sdepth = vtkSmartPointer<vtkIntArray>::New();
  sdepth->SetNumberOfComponents(1);
  sdepth->Allocate(nb);
  sdepth->SetName("SubdivisionDepth");

  vtkSmartPointer<vtkIntArray> mult = vtkSmartPointer<vtkIntArray>::New();
  mult->SetNumberOfComponents(1);
  mult->Allocate(nb);
  mult->SetName("Tangency");

  vtkSmartPointer<vtkFloatArray> norm = vtkSmartPointer<vtkFloatArray>::New();
  norm->SetNumberOfComponents(3);
  norm->Allocate(nb);

  for(int i = 0; i < nb; i++)
    {
    int j = mIndex[i];
    pts->InsertNextPoint(qX(i,0), qX(i,1), qX(i,2));
    norm->InsertNextTuple3(qN(i,0), qN(i,1), qN(i,2));
    rad->InsertNextTuple1(qR(j));
    mix->InsertNextTuple1(j);
    mult->InsertNextTuple1(mtbIndex[j].size());
    if(tmpl.subDepth.size())
      sdepth->InsertNextTuple1(tmpl.subDepth[i]);
    }

  vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  pd->Allocate(bmesh->triangles.size());
  pd->SetPoints(pts);
  pd->GetPointData()->SetNormals(norm);
  pd->GetPointData()->AddArray(mix);
  pd->GetPointData()->AddArray(mult);
  pd->GetPointData()->AddArray(rad);

  if(tmpl.subDepth.size())
    pd->GetPointData()->AddArray(sdepth);

  for(int i = 0; i < bmesh->triangles.size(); i++)
    {
    vtkIdType vtx[3];
    for(int j = 0; j < 3; j++)
      vtx[j] = bmesh->triangles[i].vertices[j];
    pd->InsertNextCell(VTK_TRIANGLE, 3, vtx);
    }

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetInputData(pd);
  writer->SetFileName(fn);
  writer->Update();
}

void BCMRepQuadraticProblemBuilder::SaveMedialMesh(const char *fn)
{
  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  pts->Allocate(nm);

  vtkSmartPointer<vtkFloatArray> rad = vtkSmartPointer<vtkFloatArray>::New();
  rad->SetNumberOfComponents(1);
  rad->Allocate(nm);
  rad->SetName("Radius Function");

  vtkSmartPointer<vtkFloatArray> spoke[3];

  spoke[0] = vtkSmartPointer<vtkFloatArray>::New();
  spoke[0]->SetNumberOfComponents(3);
  spoke[0]->Allocate(nm);
  spoke[0]->SetName("Spoke1");

  spoke[1] = vtkSmartPointer<vtkFloatArray>::New();
  spoke[1]->SetNumberOfComponents(3);
  spoke[1]->Allocate(nm);
  spoke[1]->SetName("Spoke2");

  spoke[2] = vtkSmartPointer<vtkFloatArray>::New();
  spoke[2]->SetNumberOfComponents(3);
  spoke[2]->Allocate(nm);
  spoke[2]->SetName("Spoke3");

  for(int i = 0; i < nm; i++)
    {
    vnl_vector_fixed<double, 3> xm;
    for(int j = 0; j < 3; j++)
      xm[j] = qM(i,j);

    pts->InsertNextPoint(xm[0], xm[1], xm[2]);
    rad->InsertNextTuple1(qR(i));

    std::vector<int> &ix = mtbIndex[i];

    vnl_vector_fixed<double, 3> xs;
    for(int k = 0; k < std::min((int) ix.size(), 3); k++)
      {
      for(int j = 0; j < 3; j++)
        xs[j] = qX(ix[k],j) - xm[j];
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

  if(tp_med.triangle_areas.GetFlatSize())
    {
    vtkSmartPointer<vtkFloatArray> t_area = vtkSmartPointer<vtkFloatArray>::New();
    t_area->SetNumberOfComponents(1);
    t_area->Allocate(bmesh->triangles.size());
    t_area->SetName("Triangle area");

    for(int i = 0; i < bmesh->triangles.size(); i++)
      t_area->InsertNextTuple1(tp_med.triangle_areas(i));

    pd->GetCellData()->AddArray(t_area);
    }

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetInputData(pd);
  writer->SetFileName(fn);
  writer->Update();
}

void PrintLossReport(qcqp::Problem::LossReport &l_init, qcqp::Problem::LossReport &l_curr)
{
  printf("%30s   %10s %10s   %10s   %10s\n", "TERM", "RawValue", "Weight", "WgtValue", "Delta");
  double total = 0.0, total_init = 0.0;
  for(auto it: l_curr)
    {
    double val = std::get<1>(it.second);
    double weight = std::get<0>(it.second);
    double val_init = std::get<1>(l_init[it.first]);

    printf("%30s   %10.4f %10.4f   %10.4f   %10.4f\n",
           it.first.c_str(), val, weight, val * weight, val * weight - val_init * weight);

    total += val * weight;
    total_init += val_init * weight;
    }
  printf("%30s   %10.4s %10.4s   %10.4f   %10.4f\n",
         "TOTAL OBJECTIVE", "", "", total, total - total_init);
}



int main(int argc, char *argv[])
{
  // Usage help
  if(argc < 2) return usage();

  // The action specified for the program
  ProgramAction action = ACTION_NONE;
  std::string fnTemplate, fnTarget, fnOutput, fnImportSource;
  int subdivisionLevel = 0;
  double infl_radius;
  int infl_edge_label = -1;

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
    if(cmd == "-omt")
      {
      action = ACTION_FIT_TARGET_OMT;
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
    else if(cmd == "-sub-edge")
      {
      regOpts.SubdivisionEdgeOnly = true;
      }
    else if(cmd == "-sub-flat")
      {
      regOpts.SubdivisionFlat = true;
      }
    else if(cmd == "-cmr")
      {
      action = ACTION_CONVERT_CMREP;
      fnImportSource = argv[++p];
      }
    else if(cmd == "-inf")
      {
      action = ACTION_INFLATE_CMREP;
      fnImportSource = argv[++p];
      infl_radius = atof(argv[++p]);
      infl_edge_label = atoi(argv[++p]);
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
    else if(cmd == "-max-iter")
      {
      regOpts.MaxIpOptIter = atoi(argv[++p]);
      }
    else if(cmd == "-icp-iter")
      {
      regOpts.MaxICPIter = atoi(argv[++p]);
      }
    else if(cmd == "-hsllib")
      {
      regOpts.HSLDylibPath = argv[++p];
      }
    else if(cmd == "-pardisolib")
      {
      regOpts.PardisoDylibPath = argv[++p];
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
    else if(cmd == "-no-fit-self")
      {
      regOpts.SkipSelfFit = true;
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

  else if(action == ACTION_INFLATE_CMREP)
    {
    InflateMedialModelWithBranches(fnImportSource.c_str(), fnOutput.c_str(), infl_radius, infl_edge_label);
    return 0;
    }

  else if(action == ACTION_FIX_CMREP)
    {
    FixCmrepMedialMesh(fnImportSource, fnOutput);
    return 0;
    }

  else if(action == ACTION_SUBDIVIDE)
    {
    return SubdivideBoundaryRepresentation(fnTemplate, subdivisionLevel, regOpts.SubdivisionEdgeOnly, regOpts.SubdivisionFlat, fnOutput);
    }

  else if(action == ACTION_NONE)
    {
    std::cerr << "No action specified" << std::endl;
    return -1;
    }
  
  // Get the parts of the output filename
  std::string fnOutputFull = itksys::SystemTools::CollapseFullPath(fnOutput);
  std::string fnOutBase = itksys::SystemTools::GetFilenameWithoutLastExtension(fnOutputFull);
  std::string fnOutDir = itksys::SystemTools::GetParentDirectory(fnOutputFull);
  std::cout << "Output base: " << fnOutDir << "/" << fnOutBase << std::endl;

  // Load the template
  BCMTemplate tmpl;
  tmpl.Load(fnTemplate.c_str());

  // Get the number of boundary points
  int nb = tmpl.bmesh.nVertices;

  // Load the target vertex coordintates if requesting LSQ fit
  std::vector<SMLVec3d> xLSQTarget = tmpl.x;
  if(action == ACTION_FIT_SELF)
    {
    BCMTemplate lsq_target;
    lsq_target.Load(fnTarget.c_str());
    xLSQTarget = lsq_target.x;
    }

  // Load the target mesh
  vtkSmartPointer<vtkPolyData> target = ReadVTKMesh(fnTarget.c_str());

  // A buffer for making variable names
  char buffer[1024];

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
  app->Options()->SetBoolValue("print_timing_statistics", true);
  app->Options()->SetStringValue("linear_solver", regOpts.Solver);
  if(regOpts.HSLDylibPath.size())
    app->Options()->SetStringValue("hsllib", regOpts.HSLDylibPath);
  if(regOpts.PardisoDylibPath.size())
    app->Options()->SetStringValue("pardisolib", regOpts.PardisoDylibPath);

  // app->Options()->SetNumericValue("mu_init", 1e-3);
  // app->Options()->SetNumericValue("mu_target", 1e-5);
  // app->Options()->SetStringValue("mu_strategy", "adaptive");
  // app->Options()->SetStringValue("output_file", "ipopt.out");

  app->Options()->SetIntegerValue("max_iter", regOpts.MaxIpOptIter);
  if(!regOpts.UseHessian)
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");

  // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  // app->Options()->SetStringValue("derivative_test", "second-order");
  // app->Options()->SetStringValue("derivative_test_print_all", "yes");

  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded)
    {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
    }

  // Print a message describing IpOpt
  app->PrintCopyrightMessage();

  // Keep track of the running optimal values and optimal absis coefficients
  std::vector<SMLVec3d> x_opt = tmpl.x, xc_opt;

  // Do we need to perform the initial self-fit?
  if(regOpts.SkipSelfFit == false || action == ACTION_FIT_TARGET)
    {
    // Build the new efficient quadratic problem
    BCMRepQuadraticProblemBuilder qp_builder(tmpl, regOpts, tc_bnd, tc_med);
    qp_builder.BuildGeometry(tmpl.x);
    qp_builder.BuildRegularization();
    qp_builder.BuildLeastSquaresObjective(xLSQTarget);
    qcqp::Problem &qp = qp_builder.GetProblem();
    SmartPtr<IPOptQCQPProblemInterface> q_ip = new IPOptQCQPProblemInterface(qp);

    // TODO: Create a file for dumping the constraint info
    snprintf(buffer, 1024, "%s/%s_qcqp_constraint_dump.txt", fnOutDir.c_str(), fnOutBase.c_str());
    FILE *fnQCQPConstraintDump = fopen(buffer, "wt");
    q_ip->log_constraints(fnQCQPConstraintDump);

    // Print the QP objective
    vnl_vector<double> val_init = qp.GetVariableValues();
    qcqp::Problem::LossReport lr_init = qp.GetLossReport(val_init.data_block());
    PrintLossReport(lr_init, lr_init);

    // Now the QP
    qp.SetupProblem();
    printf("QP Jacobian size (%zu, %zu), NNZ %zu\n",
           qp.GetConstraintsJacobian().GetNumberOfRows(),
           qp.GetConstraintsJacobian().GetNumberOfColumns(),
           qp.GetConstraintsJacobian().GetNumberOfSparseValues());

    printf("QP Hessian size (%zu, %zu), NNZ %zu\n",
           qp.GetHessianOfLagrangean().GetNumberOfRows(),
           qp.GetHessianOfLagrangean().GetNumberOfColumns(),
           qp.GetHessianOfLagrangean().GetNumberOfSparseValues());

 #ifdef __PRINT_STATS_

    printf("QP varible stats\n");
    int p_count = 0;
    vnl_vector<double> x_qp(qp.GetNumberOfVariables());
    for(unsigned int i = 0; i < qp.GetNumberOfVariableTensors(); i++)
      {
      auto *v = qp.GetVariableTensor(i);
      Statistics stat;
      for(int j = v->GetStartIndex(); j < v->GetStartIndex() + v->GetFlatSize(); j++)
        {
        x_qp[j] = qp.GetVariableValue(j);
        stat.add(x_qp[j]);
        }

      printf("%16s: %-4d \t %8.4f ± %8.4f \t [ %8.4f, %8.4f ]\n", v->GetName(),
             stat.n, stat.mean(), stat.sd(), stat.min(), stat.max());
      p_count += stat.n;
      }
    printf("Total: %d\n", p_count);

    std::map<std::string, Statistics> qp_con_count;
    for(unsigned int i = 0; i < qp.GetNumberOfConstraints(); i++)
      {
      qp_con_count[qp.GetConstraintName(i)].add(
            qp.EvaluateConstraint(x_qp.data_block(), i));
      }

    printf("QP varible stats\n");
    p_count = 0;
    for(auto &it : qp_con_count)
      {
      printf("%16s: %-4d \t %8.4f ± %8.4f \t [ %8.4f, %8.4f ]\n", it.first.c_str(),
             it.second.n, it.second.mean(), it.second.sd(), it.second.min(), it.second.max());
      p_count += it.second.n;
      }
    printf("Total: %d\n", p_count);

    printf("QP objective: %12.8f\n", qp.EvaluateLoss(x_qp.data_block()));

    /* Compute the grad-f statistics, this has to be done by variable and constraint type */
    std::map<std::string, Statistics> qp_gradf_count;
    vnl_vector<double> gradf_qp(qp.GetNumberOfVariables());
    qp.EvaluateLossGradient(x_qp.data_block(), gradf_qp.data_block());
    for(int i = 0; i < qp.GetNumberOfVariables(); i++)
      {
      std::string n_var = qp_var_name(qp, i);
      qp_gradf_count[n_var].add(gradf_qp[i]);
      }

    printf("QP Grad Objective stats\n");
    p_count = 0;
    for(auto &it : qp_gradf_count)
      {
      printf("%22s: %-4d \t %8.4f ± %8.4f \t [ %8.4f, %8.4f ]\n", it.first.c_str(),
             it.second.n, it.second.mean(), it.second.sd(), it.second.min(), it.second.max());
      p_count += it.second.n;
      }
    printf("Total: %d\n", p_count);

    /* Compute the Jacobian statistics, this has to be done by variable and constraint type */
    auto &jqp = qp.GetConstraintsJacobian();
    std::map<std::string, Statistics> qp_cv_count;
    for(int i = 0; i < jqp.GetNumberOfRows(); i++)
      {
      std::string n_con = qp.GetConstraintName(i);
      for(auto ri = jqp.Row(i); !ri.IsAtEnd(); ++ri)
        {
        std::string n_var = qp_var_name(qp, ri.Column());
        std::string name = n_con + "/" + n_var;
        auto &wz = ri.Value();
        double val = wz.z;
        for(auto it_w : wz.w)
          val += qp.GetVariableValue(std::get<0>(it_w)) * std::get<1>(it_w);

        qp_cv_count[name].add(val);
        }
      }

    printf("QP jacobian stats\n");
    p_count = 0;
    for(auto &it : qp_cv_count)
      {
      printf("%22s: %-4d \t %8.4f ± %8.4f \t [ %8.4f, %8.4f ]\n", it.first.c_str(),
             it.second.n, it.second.mean(), it.second.sd(), it.second.min(), it.second.max());
      p_count += it.second.n;
      }
    printf("Total: %d\n", p_count);

    // Compute qp hessian stats
    std::map<std::string, Statistics> qp_hess_count;
    auto &qp_hol = qp.GetHessianOfLagrangean();

    for(int i = 0; i < qp_hol.GetNumberOfRows(); i++)
      {
      std::string nm_i = qp_var_name(qp, i);
      for(auto it = qp_hol.Row(i); !it.IsAtEnd(); ++it)
        {
        std::string nm_j = qp_var_name(qp, it.Column());
        auto &wz = it.Value();
        double hval = wz.z;
        for(auto qtw : wz.w)
          hval += std::get<1>(qtw);
        qp_hess_count[nm_i+","+nm_j].add(hval);
        }
      }

    printf("QP Hessian stats\n");
    p_count = 0;
    for(auto &it : qp_hess_count)
      {
      printf("%22s: %-4d \t %8.4f ± %8.4f \t [ %8.4f, %8.4f ]\n", it.first.c_str(),
             it.second.n, it.second.mean(), it.second.sd(), it.second.min(), it.second.max());
      p_count += it.second.n;
      }
    printf("Total: %d\n", p_count);

#endif // __PRINT_STATS

    // Ask Ipopt to solve the problem
    app->Options()->SetStringValue("derivative_test", "none");
    status = app->OptimizeTNLP(GetRawPtr(q_ip));
    if(status < 0)
      {
      printf("\n\n*** Error %d during optimization!\n", (int) status);
      return -1;
      }

    // Get the current optimal X values
    for(unsigned int i = 0; i < nb; i++)
      x_opt[i] = qp_builder.qX(i).as_vector<3>();

    xc_opt.resize(qp_builder.qXC.GetFlatSize() / 3);
    for(unsigned int i = 0; i < xc_opt.size(); i++)
      xc_opt[i] = qp_builder.qXC(i).as_vector<3>();

    // Print the QP objective
    vnl_vector<double> val_final = qp.GetVariableValues();
    qcqp::Problem::LossReport lr_final = qp.GetLossReport(val_final.data_block());
    PrintLossReport(lr_init, lr_final);

    // Save the current state
    snprintf(buffer, 1024, "%s/%s_fit2tmp_bnd.vtk", fnOutDir.c_str(), fnOutBase.c_str());
    qp_builder.SaveBoundaryMesh(buffer);

    snprintf(buffer, 1024, "%s/%s_fit2tmp_med.vtk", fnOutDir.c_str(), fnOutBase.c_str());
    qp_builder.SaveMedialMesh(buffer);
    }

  // Continue if in ICP or OMT mode
  if(action == ACTION_FIT_TARGET || action == ACTION_FIT_TARGET_OMT)
    {
    // Create the closest point finder
    std::shared_ptr<ClosestPointMatcher> cpmatcher;
    if(ACTION_FIT_TARGET)
      {
      cpmatcher = std::shared_ptr<ClosestPointMatcher>(new ClosestPointMatcher(target, 32));
      }

    std::shared_ptr<SinkhornIteration> sinkhorn;
    std::shared_ptr<TriangleMesh> target_mesh;
    if(ACTION_FIT_TARGET_OMT)
      {
      // Generate the target mesh
      target_mesh = std::shared_ptr<TriangleMesh>(new TriangleMesh());
      TriangleMeshGenerator m_target_gen(target_mesh.get(), target->GetNumberOfPoints());
      for(int i = 0; i < target->GetNumberOfCells(); i++)
        {
        vtkCell *c = target->GetCell(i);
        if(c->GetNumberOfPoints() != 3)
          throw MedialModelException("Non-triangle cell in OMT target mesh");
        m_target_gen.AddTriangle(c->GetPointId(0), c->GetPointId(1), c->GetPointId(2));
        }
      m_target_gen.GenerateMesh();

      // Get the target coordinates
      std::vector<SMLVec3d> x_target(target->GetNumberOfPoints());
      for(int i = 0; i < target->GetNumberOfPoints(); i++)
        x_target[i].set(target->GetPoint(i));

      // Generate the Sinkhorn object
      sinkhorn = std::shared_ptr<SinkhornIteration>(
            new SinkhornIteration(&tmpl.bmesh, target_mesh.get(), x_target));
      }

    // Repeat this several times
    for(int i = 0; i < regOpts.MaxICPIter; i++)
      {
      // Create a new QP to solve
      BCMRepQuadraticProblemBuilder qp_builder_icp(tmpl, regOpts, tc_bnd, tc_med);
      qp_builder_icp.BuildGeometry(x_opt);
      qp_builder_icp.BuildRegularization(xc_opt);

      if(ACTION_FIT_TARGET)
        {
        // The scale between the two distance functions
        double recip_scale = nb * 1.0 / cpmatcher->GetNumberOfTargetPointsUsed();
        qp_builder_icp.BuildICPObjective(*cpmatcher, 1.0, recip_scale);
        }
      else
        {
        // Perform Sinkhorn iterations
        printf("Performing Sinkhorn iterations\n");
        sinkhorn->PerformIteration(x_opt, 1.0, 0.04, 100, true);
        printf("done\n");

        // Get the target locations
        std::vector<SMLVec3d> xTriangleCenterTargets(tmpl.bmesh.triangles.size());
        sinkhorn->ComputeTriangleCenterMatches(xTriangleCenterTargets);

        qp_builder_icp.BuildOMTObjective(*sinkhorn);

        // Write target locations
        vtkSmartPointer<vtkPolyData> pd = ReadVTKMesh(fnTemplate.c_str());
        vtkSmartPointer<vtkFloatArray> vMatch = vtkSmartPointer<vtkFloatArray>::New();
        vMatch->SetNumberOfComponents(3);
        vMatch->Allocate(3 * pd->GetNumberOfPoints());
        vMatch->SetName("ClosestPoint");
        for(unsigned int i = 0; i < xTriangleCenterTargets.size(); i++)
          vMatch->InsertNextTuple3(xTriangleCenterTargets[i][0],xTriangleCenterTargets[i][1],xTriangleCenterTargets[i][2]);
        pd->GetCellData()->AddArray(vMatch);
        vtkSmartPointer<vtkPolyDataWriter> wr =
            vtkSmartPointer<vtkPolyDataWriter>::New();
        wr->SetInputData(pd);
        snprintf(buffer, 1024, "%s/%s_omt_%02d_facetarget.vtk", fnOutDir.c_str(), fnOutBase.c_str(), i);
        wr->SetFileName(buffer);
        wr->Update();
        }

      qcqp::Problem &qp = qp_builder_icp.GetProblem();
      SmartPtr<IPOptQCQPProblemInterface> q_ip = new IPOptQCQPProblemInterface(qp);

      // Now the QP
      qp.SetupProblem();
      printf("QP Jacobian size (%zu, %zu), NNZ %zu\n",
             qp.GetConstraintsJacobian().GetNumberOfRows(),
             qp.GetConstraintsJacobian().GetNumberOfColumns(),
             qp.GetConstraintsJacobian().GetNumberOfSparseValues());

      printf("QP Hessian size (%zu, %zu), NNZ %zu\n",
             qp.GetHessianOfLagrangean().GetNumberOfRows(),
             qp.GetHessianOfLagrangean().GetNumberOfColumns(),
             qp.GetHessianOfLagrangean().GetNumberOfSparseValues());

      // Print the QP objective
      vnl_vector<double> val_init = qp.GetVariableValues();
      qcqp::Problem::LossReport lr_init = qp.GetLossReport(val_init.data_block());
      PrintLossReport(lr_init, lr_init);

      // Solve this qp
      status = app->OptimizeTNLP(GetRawPtr(q_ip));
      if(status < 0)
        {
        printf("\n\n*** Error %d during optimization!\n", (int) status);
        return -1;
        }

      // Get the current optimal X values
      for(unsigned int i = 0; i < nb; i++)
        x_opt[i] = qp_builder_icp.qX(i).as_vector<3>();

      xc_opt.resize(qp_builder_icp.qXC.GetFlatSize() / 3);
      for(unsigned int i = 0; i < xc_opt.size(); i++)
        xc_opt[i] = qp_builder_icp.qXC(i).as_vector<3>();

      // Print the QP objective
      vnl_vector<double> val_final = qp.GetVariableValues();
      qcqp::Problem::LossReport lr_final = qp.GetLossReport(val_final.data_block());
      PrintLossReport(lr_init, lr_final);
      
      // Save the current state
      snprintf(buffer, 1024, "%s/%s_icp_%02d_bnd.vtk", fnOutDir.c_str(), fnOutBase.c_str(), i);
      qp_builder_icp.SaveBoundaryMesh(buffer);

      snprintf(buffer, 1024, "%s/%s_icp_%02d_med.vtk", fnOutDir.c_str(), fnOutBase.c_str(), i);
      qp_builder_icp.SaveMedialMesh(buffer);
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

    //SaveGradient(p, X, obj, "grad_obj_after.vtk");
    //SaveGradient(p, X, objSqDist, "grad_obj_sqdist_after.vtk");
    //SaveGradient(p, X, new ScalarProduct(p, objRecipSqDist, recip_scale), "grad_obj_recipsqdist_after.vtk");
    //SaveGradient(p, X, new ScalarProduct(p, objBasisResidual, regOpts.Weight), "grad_obj_residual_after.vtk");

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

  // fclose(fnConstraintDump);
  return (int) status;
}
