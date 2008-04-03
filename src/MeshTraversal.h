#ifndef __MeshTraversal_h_
#define __MeshTraversal_h_


/**
 * A representation of a triangle in a mesh
 */
struct Triangle
{
  // Index of the triangle's vertices
  size_t vertices[3];

  // Index of the neighbors
  size_t neighbors[3];

  // Each edge is associated with an index of the vertex opposite to it.
  // This value tells us for each edge its index in the adjacent triangle
  short nedges[3];

  // Initializes to dummy values
  Triangle();
};

/**
 * Each vertex is associated with a walk. A walk starts at some triangle and
 * ends at a different triangle. We store the starting and ending points of
 * the walk, as well as the number of steps in the walk
 */
struct Vertex
{
  // First triangle in the walk and this vertex's index in it
  size_t t0; short v0;

  // Second triangle in the walk and this vertex's index in it
  size_t t1; short v1;

  // The size of the walk (number of triangles) - equal to valence for
  // internal vertices but not for boundary ones!!!
  size_t n;

  // Whether the vertex is on the boundary.
  bool bnd;

  // Is this a boundary vertex
  bool IsBoundary() const { return bnd; }

  // Get the valence of the vertex. For internal vertices it is equal to the
  // number of triangles that share the vertex, but for boundary vertices it
  // is equal to 1 plus that (when the mesh has a disk topology)
  size_t Valence() const { return bnd ? n + 1 : n; }

  // Constructors
  Vertex(size_t it0, short iv0, size_t it1, short iv1, size_t in, bool ibnd)
    : t0(it0), v0(iv0), t1(it1), v1(iv1), n(in), bnd(ibnd) {}
  Vertex() : t0(NOID), t1(NOID), v0(-1), v1(-1), n(0), bnd(false) {}
};

// Information describing vertex neighborhood relationship
struct NeighborInfo
{
  // Triangle in front and behind
  size_t tFront, tBack;

  // Index of the vertex relative to these triangles
  short vFront, vBack;

  // Constructor and destructor
  NeighborInfo() : tFront(NOID), tBack(NOID), vFront(-1), vBack(-1) {}
  NeighborInfo(size_t tf, short vf, size_t tb, short vb)
    : tFront(tf), vFront(vf), tBack(tb), vBack(vb) {}
};

struct TriangleMesh
{
  typedef vector<Triangle>::iterator TriangleIt;

  // List of triangles in this mesh level
  vector<Triangle> triangles;

  // Number of vertices at this mesh level
  size_t nVertices;

  // Sparse matrix representing the neighborhood relationship between
  // vertices in the mesh. This information allows easy walk-line iteration
  // around vertices.
  ImmutableSparseArray<NeighborInfo> nbr;

  // Return the valence of a vertex
  size_t GetVertexValence(size_t ivtx)
    { return nbr.GetRowIndex()[ivtx+1] - nbr.GetRowIndex()[ivtx]; }

  // Check if the vertex is on the boundary
  bool IsVertexInternal(size_t ivtx) const
    { return nbr.GetSparseData()[nbr.GetRowIndex()[ivtx]].tBack != NOID; }

  // Constructor
  MeshLevel()
    { parent = NULL; nVertices = 0; }
};



#endif
