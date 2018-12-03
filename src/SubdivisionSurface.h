#ifndef __SubdivisionSurface_h_
#define __SubdivisionSurface_h_

#include "MeshTraversal.h"

class vtkPolyData;

using namespace std;


/**
 * Subdivision surface representation
 */
class SubdivisionSurface
{
public:

  struct MeshLevel : public TriangleMesh
  {

    // Parent mesh level
    const MeshLevel *parent;

    // Sparse matrix mapping the mesh level to the parent level
    ImmutableSparseMatrix<double> weights;

    // Set the mesh level as the 'root'. This means that its weight matrix
    // becomes identity and its parent is null
    void SetAsRoot()
      { parent = NULL; weights.SetIdentity(nVertices); }

    // Constructor
    MeshLevel() : TriangleMesh()
      { parent = NULL; }

    // Constructor that makes a copy of a triangle mesh
    MeshLevel(const TriangleMesh &other)
      {
      *((TriangleMesh *)(this)) = other;
      this->SetAsRoot();
      }
  };

  /** Subdivide a mesh level once */
  static void Subdivide(const MeshLevel *src, MeshLevel *dst, bool flat_mode = false);

  /** Subdivide selected triangles on the mesh */
  static void SubdivideSelected(
    const MeshLevel *parent, MeshLevel *child, std::set<size_t> tsel, bool flat_mode = false);

  /**
   * Subdivide a mesh level n times. The intermediate levels will be
   * discarded.
   */
  static void RecursiveSubdivide(const MeshLevel *src, MeshLevel *dst, size_t n, bool flat_mode = false);

  
  /**
   * Recursive subdivision, but applied only to the triangles along the boundary
   * of the mesh.
   */
  static void RecursiveSubdivideBoundary(const MeshLevel *src, MeshLevel *dst, size_t n, bool flat_mode = false);

  /** Import a mesh from a VTK mesh */
  static void ImportLevelFromVTK(vtkPolyData *, MeshLevel &dest);

  /** Export a mesh to VTK (only sets the cells, not the points) */
  static void ExportLevelToVTK(const MeshLevel &src, vtkPolyData *mesh);

  /** Apply a subdivision to a vtk mesh */
  static void ApplySubdivision(
    vtkPolyData *src, vtkPolyData *target, MeshLevel &m);

  /**
   * Apply subdivision to arbitrary double data (in strides of nComp
   * components)
   */
  static void ApplySubdivision(const double *xsrc, double *xdst, size_t nComp, MeshLevel &m);

  /** Apply subdivision to raw data */
  static void ApplySubdivision(
    SMLVec3d *xSrc, double *rhoSrc,
    SMLVec3d *xTrg, double *rhoTrg, MeshLevel &m);

  /** Load a mesh level from a registry */
  static void LoadMeshLevel(Registry &registry);

  /** Test the correctness of a mesh level */
  static bool CheckMeshLevel(MeshLevel *mesh);

  /** Get a sub-mesh in which triangles have a certain label */
  static void PickTrianglesWithLabel(const MeshLevel &src, size_t label,
                                     MeshLevel &dst, std::vector<size_t> &vtx_full_to_sub);

private:
  // Mutable sparse matrix
  typedef vnl_sparse_matrix<double> MutableSparseMatrix;

  // Vertex assignment function (visit vertices to assign labels)
  static void RecursiveAssignVertexLabel(MeshLevel *mesh, size_t t, size_t v, size_t id);

  // Set the weights for an even vertex
  static void SetEvenVertexWeights(
    MutableSparseMatrix &W, size_t ivc, 
    const MeshLevel *parent, size_t ivp, bool flat_mode);

  // Set the weights for an odd vertex
  static void SetOddVertexWeights(
    MutableSparseMatrix &W, size_t ivc,
    const MeshLevel *parent, size_t t, size_t v, bool flat_mode);

  // Get boundary triangles
  static std::set<size_t> GetBoundaryTriangles(const MeshLevel *src);
};

#endif
