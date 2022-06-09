#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/normal.h>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/clustering.h>

// input output
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

using namespace vcg;
using namespace std;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<  Use<MyVertex>::AsVertexType, Use<MyFace>::AsFaceType>{};
class MyVertex  : public Vertex< MyUsedTypes, vertex::VFAdj, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::VFAdj, face::Normal3f, face::VertexRef, face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh<vector<MyVertex>, vector<MyFace> > {};

int main(int argc,char ** argv)
{
  if(argc<4)
  {
    printf("Usage: poisson_disk <input> <output> <samples>\n");
    return 0;
  }

  int nSamples = atoi(argv[3]);

  MyMesh m;

  //open a mesh
  int err = tri::io::Importer<MyMesh>::Open(m,argv[1]);
  if(err) { // all the importers return 0 in case of success
    printf("Error in reading %s: '%s'\n",argv[1], tri::io::Importer<MyMesh>::ErrorMsg(err));
    exit(-1);
  }

  // some cleaning to get rid of bad file formats like stl that duplicate vertexes..
  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
  int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
  printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[1]);

  tri::Allocator<MyMesh>::CompactEveryVector(m);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  tri::UpdateBounding<MyMesh>::Box(m);
  tri::UpdateFlags<MyMesh>::FaceBorderFromVF(m);
  tri::UpdateFlags<MyMesh>::VertexBorderFromFaceBorder(m);

  typedef tri::SurfaceSampling<MyMesh,tri::TrivialSampler<MyMesh> > MySampling;
  MyMesh::ScalarType radius = MySampling::ComputePoissonDiskRadius(m,nSamples);

  cout << "Radius = " << radius << endl;

  tri::TrivialSampler<MyMesh> mps;
  MySampling::PoissonDiskParam pp;

  pp.radiusVariance = 1.0;
  pp.preGenFlag=false;
  pp.geodesicDistanceFlag=false;
  pp.bestSampleChoiceFlag=true;
  pp.bestSamplePoolSize=10;

  MySampling::PoissonDiskPruning(mps, m, radius, pp);

  MyMesh sampleMesh;
  tri::BuildMeshFromCoordVector(sampleMesh,mps.SampleVec()); 

  cout << "Pruning done, new vertices: " << sampleMesh.vn << endl;

  //LaplacianSmooth(m,atoi(argv[2]));
  tri::io::ExporterPLY<MyMesh>::Save(sampleMesh,argv[2],false);

  return 0;
}
