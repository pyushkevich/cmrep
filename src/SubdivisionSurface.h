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
  };

  /** Subdivide a mesh level once */
  static void Subdivide(const MeshLevel *src, MeshLevel *dst);

  /**
   * Subdivide a mesh level n times. The intermediate levels will be
   * discarded.
   */
  static void RecursiveSubdivide(const MeshLevel *src, MeshLevel *dst, size_t n);

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
  static bool CheckMeshLevel(MeshLevel &mesh);

private:
  // Mutable sparse matrix
  typedef vnl_sparse_matrix<double> MutableSparseMatrix;

  // Vertex assignment function (visit vertices to assign labels)
  static void RecursiveAssignVertexLabel(MeshLevel *mesh, size_t t, size_t v, size_t id);

  // Set the weights for an even vertex
  static void SetEvenVertexWeights(MutableSparseMatrix &W,
    const MeshLevel *parent, MeshLevel *child, size_t t, size_t v);

  // Set the weights for an odd vertex
  static void SetOddVertexWeights(MutableSparseMatrix &W,
    const MeshLevel *parent, MeshLevel *child, size_t t, size_t v);

};

#endif
