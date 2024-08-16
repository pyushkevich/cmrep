#include "VCGTriMesh.h"
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkCellArrayIterator.h>
#include <vcg/complex/algorithms/clean.h>


using namespace vcg;
using namespace std;

void VCGTriMesh::ImportFromVTK(vtkPolyData *pd)
{
  int np = pd->GetNumberOfPoints();

  // Remove existing elements
  m_Mesh.Clear();

  // Create adder functions
  Mesh::VertexIterator vi = vcg::tri::Allocator<Mesh>::AddVertices(m_Mesh, np);
  Mesh::FaceIterator fi = vcg::tri::Allocator<Mesh>::AddFaces(m_Mesh, pd->GetPolys()->GetNumberOfCells());

  // Container for created vertex pointers
  vector<Mesh::VertexPointer> ivp(np);
  for(unsigned int i = 0; i < np; i++)
    {
    ivp[i] = &*vi;
    auto *p = pd->GetPoint(i);
    vi->P() = Mesh::CoordType(p[0], p[1], p[2]);
    ++vi;
    }

  // Add the faces
  int bad_faces = 0;
  auto iter = vtk::TakeSmartPointer(pd->GetPolys()->NewIterator());
  for (iter->GoToFirstCell(); !iter->IsDoneWithTraversal(); iter->GoToNextCell())
  {
    // do work with iter
    const vtkIdType *ids;
    vtkIdType npts;
    iter->GetCurrentCell(npts, ids);
    if(npts != 3)
      {
      ++bad_faces;
      continue;
      }
    fi->V(0)=ivp[ids[0]];
    fi->V(1)=ivp[ids[1]];
    fi->V(2)=ivp[ids[2]];
    fi++;
  }

  if(bad_faces > 0)
    std::cerr << "Conversion from VTK to VCG encountered " << bad_faces << " non-triangular faces" << std::endl;
  }

#include <vtkPointData.h>
#include <vtkFloatArray.h>

void VCGTriMesh::ExportToVTK(vtkPolyData *pd)
{
  // Point array
  vtkNew<vtkPoints> pts;

  // Normal array
  vtkNew<vtkFloatArray> normals;
  normals->SetNumberOfComponents(3);

  std::vector<int> VertexId(m_Mesh.vert.size());
  int numvert = 0;
  for(auto vi = m_Mesh.vert.begin(); vi != m_Mesh.vert.end(); ++vi)
    {
    if( !(*vi).IsD() )
      {
      const auto &p = vi->cP();
      const auto &n = vi->cN();
      pts->InsertNextPoint(p[0], p[1], p[2]);
      normals->InsertNextTuple3(n[0], n[1], n[2]);
      VertexId[vi - m_Mesh.vert.begin()] = numvert++;
      }
    }

  pd->SetPoints(pts);
  pd->GetPointData()->SetNormals(normals);

  // Face array
  vtkNew<vtkCellArray> faces;
  for(auto fi = m_Mesh.face.begin(); fi != m_Mesh.face.end(); ++fi)
    {
    if( !(*fi).IsD() )
      {
      faces->InsertNextCell(fi->VN());
      for(int k=0; k< fi->VN(); k++)
        faces->InsertCellPoint(VertexId[tri::Index(m_Mesh, fi->V(k))]);
      }
    }

  pd->SetPolys(faces);
}

void VCGTriMesh::CleanMesh()
{
  // some cleaning to get rid of bad file formats like stl that duplicate vertexes..
  int dup = tri::Clean<Mesh>::RemoveDuplicateVertex(m_Mesh);
  int unref = tri::Clean<Mesh>::RemoveUnreferencedVertex(m_Mesh);
  printf("Removed %i duplicate and %i unreferenced vertices from mesh\n", dup, unref);

  // Remove degenerate and zero area faces
  int dupf = tri::Clean<Mesh>::RemoveDuplicateFace(m_Mesh);
  int degf = tri::Clean<Mesh>::RemoveDegenerateFace(m_Mesh);
  int zaf = tri::Clean<Mesh>::RemoveZeroAreaFace(m_Mesh);
  printf("Removed %i duplicate, %i degenerate, and %i zero area faces from mesh\n", dupf, degf, zaf);

  tri::Allocator<Mesh>::CompactEveryVector(m_Mesh);
  tri::UpdateTopology<Mesh>::VertexFace(m_Mesh);
  tri::UpdateBounding<Mesh>::Box(m_Mesh);
  tri::UpdateFlags<Mesh>::FaceBorderFromVF(m_Mesh);
  tri::UpdateFlags<Mesh>::VertexBorderFromFaceBorder(m_Mesh);
}

void VCGTriMesh::RecomputeNormals()
{
  // updateBoxAndNormals() in meshlab
  tri::UpdateBounding<Mesh>::Box(m_Mesh);
  if(m_Mesh.fn > 0)
    {
    tri::UpdateNormal<Mesh>::PerFaceNormalized(m_Mesh);
    tri::UpdateNormal<Mesh>::PerVertexAngleWeighted(m_Mesh);
    }

  // Normal cleaning
  tri::UpdateNormal<Mesh>::NormalizePerFace(m_Mesh);
  tri::UpdateNormal<Mesh>::PerVertexFromCurrentFaceNormal(m_Mesh);
  tri::UpdateNormal<Mesh>::NormalizePerVertex(m_Mesh);
  }

tri::TriEdgeCollapseQuadricParameter VCGTriMesh::GetDefaultQuadricEdgeCollapseRemeshingParameters()
{
  vcg::tri::TriEdgeCollapseQuadricParameter qparams;
  qparams.QualityThr = 0.3;
  qparams.PreserveBoundary=true;
  qparams.BoundaryQuadricWeight=1.0;
  qparams.PreserveTopology=true;
  qparams.QualityWeight=false;
  qparams.NormalCheck=true;
  qparams.OptimalPlacement=true;
  qparams.QualityQuadric=true;
  qparams.QualityQuadricWeight=0.001;

  return qparams;
  }

void VCGTriMesh::QuadricEdgeCollapseRemeshing(
    double reduction_factor,
    vcg::tri::TriEdgeCollapseQuadricParameter qparams)
{
  // decimator initialization
  vcg::LocalOptimization<VCGTriMesh::Mesh> DeciSession(m_Mesh, &qparams);

  // specialization
  typedef vcg::tri::BasicVertexPair<VCGTriMesh::Vertex> VertexPair;
  class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric<
      VCGTriMesh::Mesh, VertexPair,
      MyTriEdgeCollapse,
      vcg::tri::QInfoStandard<VCGTriMesh::Vertex>  >
  {
  public:
    typedef vcg::tri::TriEdgeCollapseQuadric<
        VCGTriMesh::Mesh, VertexPair,
        MyTriEdgeCollapse,
        vcg::tri::QInfoStandard<VCGTriMesh::Vertex>  > TECQ;
    typedef  VCGTriMesh::Mesh::VertexType::EdgeType EdgeType;
    inline MyTriEdgeCollapse(const VertexPair &p, int i, vcg::BaseParameterClass *pp) :TECQ(p,i,pp){}
  };

  // Target number of vertices
  int n_target = reduction_factor < 1.0
      ? (int) (m_Mesh.FN() * reduction_factor)
      : (int) reduction_factor;

  int t1=clock();
  DeciSession.Init<MyTriEdgeCollapse>();
  int t2=clock();
  printf("Simplifying to reduce from %7i to %7i faces\n", m_Mesh.FN(), n_target);

  DeciSession.SetTargetSimplices(n_target);
  DeciSession.SetTimeBudget(2.0f);
  // DeciSession.SetTargetOperations(100000);
  // if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

  double TargetError = std::numeric_limits<double>::max();
  while( DeciSession.DoOptimization() && m_Mesh.FN() > n_target && DeciSession.currMetric < TargetError)
    printf("Current Mesh size %7i heap sz %9i err %9g \n", m_Mesh.FN(), int(DeciSession.h.size()),DeciSession.currMetric);

  int t3=clock();
  printf("mesh (%d,%d) Error %g \n", m_Mesh.VN(), m_Mesh.FN(), DeciSession.currMetric);
  printf("\nCompleted in (%5.3f+%5.3f) sec\n", float(t2-t1)/CLOCKS_PER_SEC, float(t3-t2)/CLOCKS_PER_SEC);

  DeciSession.Finalize<MyTriEdgeCollapse>();
}
