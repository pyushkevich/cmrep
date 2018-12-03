#ifndef __MeshTraversal_h_
#define __MeshTraversal_h_

#include <vector>
#include <set>
#include "smlmath.h"
#include "Registry.h"
#include "SparseMatrix.h"

using namespace std;

// Useful global-level inline routines.
// TODO: are there faster ops for this?
inline short ror(short x) { return (x + 1) % 3; }
inline short rol(short x) { return (x + 2) % 3; }

// Constant used for unreasonable size_t objects
#define NOID 0xffffffff


/**
 * A representation of a triangle in a mesh
 */
struct Triangle
{
  // Index of the triangle's vertices
  size_t vertices[3];

  // Index of the neighbors
  size_t neighbors[3];

  // Optional label of the triangle. Will be propagated to children
  size_t label;

  // Each edge is associated with an index of the vertex opposite to it.
  // This value tells us for each edge its index in the adjacent triangle
  short nedges[3];

  // Initializes to dummy values
  Triangle()
    {
    vertices[0] = vertices[1] = vertices[2] = NOID;
    neighbors[0] = neighbors[1] = neighbors[2] = NOID;
    nedges[0] = nedges[1] = nedges[2] = -1;
    label = NOID;
    }
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
  Vertex() : t0(NOID), v0(-1), t1(NOID), v1(-1), n(0), bnd(false) {}
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
    : tFront(tf), tBack(tb), vFront(vf), vBack(vb) {}
};

class TriangleMesh
{

public:

  typedef vector<Triangle>::iterator TriangleIt;

  // List of triangles in this mesh level
  vector<Triangle> triangles;

  // Number of vertices at this mesh level
  size_t nVertices;

  // Sparse matrix representing the neighborhood relationship between
  // vertices in the mesh. This information allows easy walk-line iteration
  // around vertices.
  typedef ImmutableSparseArray<NeighborInfo> NeighborMatrix;
  
  // Return the valence of a vertex
  size_t GetVertexValence(size_t ivtx)
    { return nbr.GetRowIndex()[ivtx+1] - nbr.GetRowIndex()[ivtx]; }

  // Check if the vertex is on the boundary
  bool IsVertexInternal(size_t ivtx) const
    { return nbr.GetSparseData()[nbr.GetRowIndex()[ivtx]].tBack != NOID; }

  // Constructor
  TriangleMesh()
    { nVertices = 0; }

  /** 
   * Compute the walks in a mesh level. This associates each vertex with a 
   * triangle, which makes subsequent processing a lot easier
   */
  void ComputeWalks();

  /**
   * This method is going to make a mesh Delaunay by flipping edges. It takes
   * the coordinates of the vertices in the mesh as parameters.
   */
  void MakeDelaunay(vnl_vector_fixed<double,3> *X);

  /** Return a constant reference to the sparse array structure */
  const NeighborMatrix &GetNeighborMatrix()
    { return nbr; }

private:
  NeighborMatrix nbr;

  friend class TriangleMeshGenerator;
  friend class EdgeWalkAroundVertex;
};




/**
 * A simple class for generating triangle mesh objects from lists of
 * vertices. Simply pass in the mesh that you want to initialize, then
 * add all the triangles one by one, and call Generate();
 */
class TriangleMeshGenerator
{
public:

  // Constructor, pass in the mesh you want to 'generate'
  TriangleMeshGenerator(TriangleMesh *target, size_t nVertices);
    
  // Add a triangle to the generator
  void AddTriangle(size_t v0, size_t v1, size_t v2, size_t label = NOID);

  // Populate the triangle mesh
  void GenerateMesh();

private:

  typedef pair<size_t, size_t> HalfEdge;
  typedef pair<size_t, short> TriangleRep;
  typedef map<HalfEdge, TriangleRep> TriangleMap;
  TriangleMap tmap;

  TriangleMesh *xTargetMesh;
};


/**
 * This is an iterator that represents a walk around the vertex where at each
 * point in time you are standing on one edge. To your right is a vertex that
 * you are walking around, to your left is a neighbor vertex. In front of you
 * or behind of you may be a triangle (unless you are at the beginning). The
 * number of steps is equal to the valence of the vertex
 */
class EdgeWalkAroundVertex
{
public:
  /** Constructor: takes a vertex in a fully initialized mesh level */
  EdgeWalkAroundVertex(const TriangleMesh *mesh, size_t iVertex)
    {
    this->mesh = mesh;
    this->iVertex = iVertex;
    this->pos = mesh->nbr.GetRowIndex()[iVertex];
    this->end = mesh->nbr.GetRowIndex()[iVertex + 1];
    }

  /** Check if this walk is closed or open (boundary vertex) */
  bool IsOpen()
    { return mesh->nbr.GetSparseData()[end - 1].tFront == NOID; }

  /** Return the valence of the vertex in the center */
  size_t Valence()
    { 
    return mesh->nbr.GetRowIndex()[iVertex + 1] 
      - mesh->nbr.GetRowIndex()[iVertex]; 
    }

  /**
   * Are we at the end of the walk? The end is reached when you have made a
   * full circle or moved past the end of the mesh. No operations should be
   * performed once you are at the end of the walk
   */
  bool IsAtEnd() const
    { return pos >= end; }

  /** Increment operator */
  EdgeWalkAroundVertex &operator++ ()
    { pos++; return *this; }

  /** Go to the last edge in the walk (next step puts you at end) */
  void GoToLastEdge()
    { pos = end - 1; }

  /** Go to the first edge in the walk (next step puts you at end) */
  void GoToFirstEdge()
    { pos = mesh->nbr.GetRowIndex()[iVertex]; }

  /** Get the id of the 'fixed' vertex. This is always the same id */
  size_t FixedVertexId()
    { return iVertex; }

  /** Get the moving vertex id */
  size_t MovingVertexId()
    { return mesh->nbr.GetColIndex()[pos]; }

  /** Get the triangle ahead */
  size_t TriangleAhead()
    { return mesh->nbr.GetSparseData()[pos].tFront; }

  /** Get the triangle behind */
  size_t TriangleBehind()
    { return mesh->nbr.GetSparseData()[pos].tBack; }

  /** Get the index of the fixed vertex in the triangle ahead */
  short FixedVertexIndexInTriangleAhead()
    { return mesh->nbr.GetSparseData()[pos].vFront; }

  /** Get the index of the fixed vertex in the triangle behind */
  short FixedVertexIndexInTriangleBehind()
    { return mesh->nbr.GetSparseData()[pos].vBack; }

  /** Get the index of the moving vertex in the triangle ahead */
  short MovingVertexIndexInTriangleAhead()
    { return ror(FixedVertexIndexInTriangleAhead()); }

  /** Get the index of the moving vertex in the triangle behind */
  short MovingVertexIndexInTriangleBehind()
    { return rol(FixedVertexIndexInTriangleBehind()); }

  /** Get the index of the face-connected vertex in triangle ahead */
  short OppositeVertexIndexInTriangleAhead()
    { return rol(FixedVertexIndexInTriangleAhead()); }

  /** Get the index of the face-connected vertex in triangle behind */
  short OppositeVertexIndexInTriangleBehind()
    { return ror(FixedVertexIndexInTriangleBehind()); }

  /** Get the id of the face-connected vertex ahead */
  size_t VertexIdAhead()
    {
    if(pos == end - 1)
      if(mesh->nbr.GetSparseData()[pos].tFront == NOID)
        return NOID;
      else
        return mesh->nbr.GetColIndex()[mesh->nbr.GetRowIndex()[iVertex]];
    else
      return mesh->nbr.GetColIndex()[pos + 1];
    }

  /** Get the id of the face-connected vertex behind */
  size_t VertexIdBehind()
    {
    if(pos == mesh->nbr.GetRowIndex()[iVertex])
      if(mesh->nbr.GetSparseData()[pos].tBack == NOID)
        return NOID;
      else
        return mesh->nbr.GetColIndex()[end - 1];
    else
      return mesh->nbr.GetColIndex()[pos - 1];
    }

  /** Get the index of the iterator's position in the sparse array */
  size_t GetPositionInMeshSparseArray() const
    { return pos; }

private:
  // The mesh being iterated around
  const TriangleMesh *mesh;

  // The index of the vertex walked around by this iterator
  size_t iVertex;

  // The position in the walk, points to the column entry in the mesh's nbr
  // sparse matrix
  size_t pos, end;
};

/**
 * A scheme for computing tangent vectors from mesh levels. It also can be used
 * to compute the limit surface using the index 2
 */
class LoopTangentScheme
{
public:

  enum Index { TANGENT_U = 0, TANGENT_V = 1, LIMIT_SURFACE = 2 };

  // Define a pair of weights (wu, wv, w_limit)
  struct Weight { double w[3]; };

  // A sparse matrix of weights (computed using Loop fmla)
  typedef ImmutableSparseArray<Weight> WeightMatrix;

  /** Pass a mesh to the tangent scheme */
  void SetMesh(TriangleMesh *mesh);

  /** Direct access to the weights. This method returns the contribution
   * of a neighbor pointed by a walk iterator to the d-tangent at a vertex */
  double GetNeighborWeight(size_t d, EdgeWalkAroundVertex &walk)
    {
    size_t i = walk.FixedVertexId();
    size_t k = walk.GetPositionInMeshSparseArray() + i + 1;
    return W.GetSparseData()[k].w[d];
    }

  /** Direct access to the weights. This method returns the contribution
   * of the vertex itself to a d-tangent. This is zero for vertices with 
   * regular valence */
  double GetOwnWeight(size_t d, size_t i)
    {   
    size_t k = W.GetRowIndex()[i];
    return W.GetSparseData()[k].w[d];
    }

  /** Get the internal weight matrix */
  const WeightMatrix &GetWeightMatrix() const
    { return W; }

  /** Return partial derivative of whatever function F over the mesh */
  template<class T> T GetPartialDerivative(size_t d, size_t v, const T *F)  
    {
    size_t *ri = W.GetRowIndex();
    size_t *ci = W.GetColIndex();
    Weight *wgt = W.GetSparseData();

    T y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += wgt[i].w[d] * F[ci[i]];

    return y;
    }

  /** Compute the directional derivative of a function F over the mesh
   * in direction d (0,1) at vertex v */
  /*
  template<class T> T GetPartialDerivative(size_t d, size_t v, const T *F)
  {
    size_t *ri = level->nbr.GetRowIndex();
    size_t *ci = level->nbr.GetColIndex();

    double y = xVtxTangentWeights[d][v] * F[v];
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += xNbrTangentWeights[d][i] * F[ci[i]];

    return y;
  }*/

protected:

  // The matrix that stores all the Loop weights
  WeightMatrix W;

};


/**
 * This class computes the Riemannian gradient of function F
 * on a mesh with vertices X. It uses Loop's formula. It can
 * also compute the Jacobian of the gradient with respect to
 * X-s and F-s. In other words,
 *
 * gradF = gradF(X, F)
 * J_X = d(gradF_i)/d(X_j)
 * J_F = d(gradF_i)/d(F_j)
 */
class MeshGradientComputer 
{
public:
  typedef vnl_vector_fixed<double, 3> Vec;
  typedef vnl_matrix_fixed<double, 3, 3> Mat;
  typedef ImmutableSparseArray<Mat> SparseMatX;
  typedef ImmutableSparseArray<Vec> SparseMatF;

  /** Set mesh, initialize internal structures */
  void SetMesh(TriangleMesh *mesh);

  /** Compute gradient and (optionally) Jacobian for given X, F */
  void ComputeGradient(Vec *x, double *f, bool flagJacobian);

  /** A more generic method that allows arrays of structs to be passed in */
  void ComputeGradient(
    void *xVecPtr, size_t xVecStride, 
    void *xFunPtr, size_t xFunStrude,
    bool flagJacobian);

public:

  /** The jacobian with respect to X */
  SparseMatX jacX;

  /** The jacobian with respect to F */
  SparseMatF jacF;

  /** The gradient of F */
  std::vector<Vec> gradF;

private:
  TriangleMesh *xMesh;
  LoopTangentScheme xLoop;

};

/**
 * Use edge flipping to generate a Delaunay mesh from an input
 * triangle mesh. Based on the Matthew Fisher paper
 */
void MakeMeshDelaunay(TriangleMesh *src, TriangleMesh *dst);

inline ostream &operator << (ostream &out, const Triangle &t)
{
  out << "[V=(" << t.vertices[0] << "," << t.vertices[1] << "," << t.vertices[2] << ")";
  out << " N=(" << t.neighbors[0] << "," << t.neighbors[1] << "," << t.neighbors[2] << ")";
  out << " NE=(" << t.nedges[0] << "," << t.nedges[1] << "," << t.nedges[2] << ")]";
  return out;
}

#endif
