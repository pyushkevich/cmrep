#ifndef VCG_WRAP_H
#define VCG_WRAP_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

class vtkPolyData;

class VCGTriMesh
{
public:
  // Define the mesh type
  class Face;
  class Edge;
  class Vertex;
  struct UsedTypes : public vcg::UsedTypes<vcg::Use<Vertex>::AsVertexType, vcg::Use<Edge>::AsEdgeType, vcg::Use<Face>::AsFaceType>{};
  class Vertex     : public vcg::Vertex<UsedTypes,
                                        vcg::vertex::VFAdj, vcg::vertex::Coord3f, vcg::vertex::Mark,
                                        vcg::vertex::Qualityf, vcg::vertex::Normal3f, vcg::vertex::BitFlags  >
    {
    public:
      vcg::math::Quadric<double> &Qd() { return q; }
    private:
      vcg::math::Quadric<double> q;
     };

  class Edge       : public vcg::Edge<UsedTypes> {};
  class Face       : public vcg::Face<UsedTypes, vcg::face::VFAdj, vcg::face::Normal3f, vcg::face::VertexRef, vcg::face::BitFlags > {};
  class Mesh       : public vcg::tri::TriMesh<std::vector<Vertex>, std::vector<Face> > {};

  // Import a VTK polydata into the mesh structure
  void ImportFromVTK(vtkPolyData *pd);

  // Import a VTK polydata into the mesh structure
  void ExportToVTK(vtkPolyData *pd);

  // Clean the mesh and compute internal quantities
  void CleanMesh();

  // Recompute normals
  void RecomputeNormals();

  // Access the  mesh
  Mesh &GetMesh() { return m_Mesh; }


protected:

  Mesh m_Mesh;

};

#endif // VCG_WRAP_H
