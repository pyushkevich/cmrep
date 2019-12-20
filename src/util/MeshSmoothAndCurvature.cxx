#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/quality.h>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include<vcg/complex/algorithms/clustering.h>

// input output
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/export_stl.h>

// include VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkCell.h>
#include <vtkCellArray.h>

using namespace vcg;
using namespace std;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<  Use<MyVertex>::AsVertexType, Use<MyFace>::AsFaceType>{};
class MyVertex  : public Vertex< MyUsedTypes, vertex::VFAdj, vertex::Coord3f, vertex::Normal3f, vertex::CurvatureDirf, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::VFAdj, face::FFAdj, face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh<vector<MyVertex>, vector<MyFace> > {};

template <class TMesh>
class VTK_VCG_Interface
{
public:
  static void ReadPolyData(TMesh &m, const char *file)
    {
    // Import a mesh from VTK
    vtkSmartPointer<vtkPolyDataReader> vReader = vtkPolyDataReader::New();
    vReader->SetFileName(file);
    vReader->Update();

    // Perform triangulation
    vtkSmartPointer<vtkTriangleFilter> vTriFilter = vtkTriangleFilter::New();
    vTriFilter->SetInputConnection(vReader->GetOutputPort());
    vTriFilter->Update();

    // Get the polydata
    vtkSmartPointer<vtkPolyData> vMesh = vTriFilter->GetOutput();
    vtkSmartPointer<vtkPointData> vPointData = vMesh->GetPointData();
    vtkSmartPointer<vtkDataArray> vNormals = vPointData->GetNormals();

    // Add all the vertices
    for(int i = 0; i < vMesh->GetNumberOfPoints(); i++)
      {
      double *vP = vMesh->GetPoint(i);

      if(vNormals)
        {
        double *vN = vNormals->GetTuple(i);
        tri::Allocator<TMesh>::AddVertex(m, 
          typename TMesh::CoordType(vP[0], vP[1], vP[2]),
          typename TMesh::CoordType(vN[0], vN[1], vN[2])
        );
        }
      else
        {
        tri::Allocator<TMesh>::AddVertex(m, typename TMesh::CoordType(vP[0], vP[1], vP[2]));
        }
      }

    // Add all the faces
    vtkIdType ids[3];
    for(int i = 0; i < vMesh->GetNumberOfCells(); i++)
      {
      vtkCell *cell = vMesh->GetCell(i);
      tri::Allocator<TMesh>::AddFace(m, cell->GetPointId(0), cell->GetPointId(1), cell->GetPointId(2));
      }
    }

  enum Attr {
    A_NORMAL=0x0001,
    A_CURVATURE=0x0010,
  };

  static vtkDataArray* MakeDoubleArray(vtkPolyData *pd, int n_comp, const char *name)
    {
    vtkSmartPointer<vtkDoubleArray> arr = vtkDoubleArray::New();
    arr->SetName(name);
    arr->SetNumberOfComponents(n_comp);
    pd->GetPointData()->AddArray(arr);
    return pd->GetPointData()->GetArray(name);
    }

  static void WritePolyData(TMesh &m, const char *filename, int attr = 0)
    {
    // Create a new polydata
    vtkSmartPointer<vtkPolyData> pd = vtkPolyData::New();
    vtkSmartPointer<vtkPoints> pts = vtkPoints::New();
    pd->SetPoints(pts);

    // Add the points
    typename TMesh::VertexIterator vi;
    std::map<typename TMesh::VertexPointer, int> index;
    int ind = 0;

    // Add arrays to the points
    vtkDataArray *arrNormal = 0, 
                 *arrCurvMean = 0, *arrCurvGauss = 0, *arrCurviness = 0, *arrCurvSI = 0, *arrCurvK1 = 0, *arrCurvK2 = 0;

    if(attr & A_NORMAL)
      {
      arrNormal = MakeDoubleArray(pd, 3, "normals");
      }
    if(attr & A_CURVATURE)
      {
      arrCurvMean = MakeDoubleArray(pd, 1, "MeanCurvature");
      arrCurvGauss = MakeDoubleArray(pd, 1, "GaussianCurvature");
      arrCurviness = MakeDoubleArray(pd, 1, "Curviness");
      arrCurvSI = MakeDoubleArray(pd, 1, "ShapeIndex");
      arrCurvK1 = MakeDoubleArray(pd, 1, "Kappa1");
      arrCurvK2 = MakeDoubleArray(pd, 1, "Kappa2");
      }

    // Add the vertices and their properties
    for(vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      {
      if(!vi->IsD()) 
        {
        pts->InsertNextPoint(vi->P()[0], vi->P()[1], vi->P()[2]);
        index[&*vi] = ind++;

        if((attr & A_NORMAL)) // && vi->HasNormal())
          {
          arrNormal->InsertNextTuple3(vi->N()[0], vi->N()[1], vi->N()[2]);
          }
        if((attr & A_CURVATURE)) // && vi->HasCurvature())
          {
          // double H = vi->Kh();
          // double K = vi->Kg();
          double k1 = vi->K1(); // H + sqrt(H * H - K);
          double k2 = vi->K2(); // H - sqrt(H * H - K);
          double H = 0.5 * (k1 + k2), K = k1 * k2;
          double C = sqrt(k1 * k1 + k2 * k2);
          double S = atan2((k2 + k1) , (k2 - k1)) * 2.0 / 3.14159265359;
          arrCurvMean->InsertNextTuple1(H);
          arrCurvGauss->InsertNextTuple1(K);
          arrCurviness->InsertNextTuple1(C);
          arrCurvSI->InsertNextTuple1(S);
          arrCurvK1->InsertNextTuple1(k1);
          arrCurvK2->InsertNextTuple1(k2);
          }
        }
      }
    printf(".\n");


    // Add the faces
    vtkSmartPointer<vtkCellArray> cells = vtkCellArray::New();
    typename TMesh::FaceIterator fi;
    for(fi = m.face.begin(); fi != m.face.end(); ++fi)
      {
      if(!fi->IsD())
        {
        // Get the entries
        vtkIdType ids[3];
        ids[0] = index[fi->V(0)];
        ids[1] = index[fi->V(1)];
        ids[2] = index[fi->V(2)];
        cells->InsertNextCell(3, ids);
        }
      }

    pd->SetPolys(cells);

    // Write
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkPolyDataWriter::New();
    writer->SetInputData(pd);
    writer->SetFileName(filename);
    writer->Update();
    }
};

int usage()
{
  printf("vcg_smooth: Taubin mesh smoothing and curvature computation\n");
  printf("Usage:\n");
  printf("  vcg_smooth [options] in_mesh.vtk out_mesh.vtk\n");
  printf("Options:\n");
  printf("  -mu X.XX       : Taubin mu parameter (-0.51) \n"); 
  printf("  -lambda X.XX   : Taubin lambda parameter (0.5) \n"); 
  printf("  -iter N        : Taubin iter parameter (1000) \n"); 
  printf("  -r X.XX        : Radius for the PCA curvature method\n");
  return -1;
}

int main(int argc,char ** argv)
{
  if(argc<3) return usage();

  const char *in_mesh = argv[argc-2];
  const char *out_mesh = argv[argc-1];

  double smooth_mu = -0.51;
  double smooth_lambda = 0.5;
  int smooth_iter = 1000;
  double pca_R = 0.0;
  double curve_N = 20.0;

  // Parse arguments
  for (int i = 1; i < argc - 2; i++)
    {
    std::string arg = argv[i];
    if(arg == "-mu")
      {
      smooth_mu = atof(argv[++i]);
      }
    else if(arg == "-lambda")
      {
      smooth_lambda = atof(argv[++i]);
      }
    else if(arg == "-iter")
      {
      smooth_iter = atoi(argv[++i]);
      }
    else if(arg == "-r")
      {
      pca_R = atof(argv[++i]);
      }
    else
      {
      std::cerr << "Unknown option " << arg << std::endl;
      return usage();
      }
    }

  MyMesh m;

  // Read the VTK mesh
  printf("Reading VTK polydata from %s\n", in_mesh);
  VTK_VCG_Interface<MyMesh>::ReadPolyData(m, in_mesh);

  // some cleaning to get rid of bad file formats like stl that duplicate vertexes..
  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
  int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
  if(dup > 0 || unref > 0)
    printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[1]);

  tri::Clean<MyMesh>::FlipMesh(m);

  tri::Allocator<MyMesh>::CompactEveryVector(m);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  tri::UpdateBounding<MyMesh>::Box(m);

  tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
  int cnt = tri::UpdateSelection<MyMesh>::VertexFromFaceStrict(m);

  // Perform the smoothing
  printf("Smoothing mesh with mu=%f, lambda=%f, iter=%d\n", smooth_mu, smooth_lambda, smooth_iter);
  tri::Smooth<MyMesh>::VertexCoordTaubin(m, smooth_iter, smooth_lambda, smooth_mu, cnt>0);
  printf("Done\n");

  // Fix up the mesh
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  vcg::tri::Allocator<MyMesh>::CompactVertexVector(m);
  vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);

  int nonman = tri::Clean<MyMesh>::CountNonManifoldEdgeFF(m);
  if ( nonman > 0 ) 
    {
    printf("Mesh has %d not 2-manifold faces, cannot compute PC curvature", nonman); // text
    exit(-1);
    }

  int delvert=tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
  tri::Allocator<MyMesh>::CompactVertexVector(m);
  if(delvert > 0)
    printf( "Removed %d unreferenced vertices\n", delvert);

  // Compute PCA curvature
  if(pca_R == 0)
    pca_R = m.bbox.Diag() / curve_N;
  printf("Computing curvatures using PCA, radius %f\n", pca_R);
  tri::UpdateCurvature<MyMesh>::PrincipalDirectionsPCA(m, pca_R,  true);
  printf("Done\n");

  // Write mesh
  printf("Writing VTK polydata to %s\n", out_mesh);
  typedef VTK_VCG_Interface<MyMesh> VTKIO;
  VTKIO::WritePolyData(m, out_mesh, VTKIO::A_NORMAL | VTKIO::A_CURVATURE);

  return 0;
}
