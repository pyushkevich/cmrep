#ifndef __BranchingSubdivisionSurface_h_
#define __BranchingSubdivisionSurface_h_

#include <vector>
#include <set>
#include <map>
#include "smlmath.h"
#include "Registry.h"
#include "SparseMatrix.h"
#include "SubdivisionSurface.h"

class vtkPolyData;

using namespace std;

// Useful global-level inline routines.
// TODO: are there faster ops for this?
// inline short ror(short x) { return (x + 1) % 3; }
// inline short rol(short x) { return (x + 2) % 3; }

// Constant used for unreasonable size_t objects
#define NOID 0xffffffff

/**
 * Subdivision surface representation
 */
class BranchingSubdivisionSurface
{
public:

  /**
   * This class represents a triangle and its relationship with its
   * neighbors. A triangle has three edges, and across each of the 
   * edges, it may be adjacent to 0, 1 or 2 other triangles. This
   * structure provides inline routines to access the neighboring
   * triangles.
   */
  struct Triangle
  {
    // Indices of the three vertices
    size_t vertices[3];
    
    // Number of neighbors across each edge of the triangle
    size_t n_nbr[3];

    // Starting index into the neighbor triangle array for each edge
    size_t i_nbr[4];

    // List of neighbors, indexed by i_nbr
    size_t *nbr;

    // The index of the shared edge in each of the neighbors
    short *nedges;

    // Get the neighboring triangle across edge iEdge. Since there may
    // be multiple neighboring triangles across this edge, this method
    // takes a second parameter jNbr that selects the right one
    size_t GetNeighbor(short iEdge, size_t jNbr) const
      { return nbr[i_nbr[iEdge] + jNbr]; }

    // Get the index of the vertex iEdge in the triangle returned by
    // the GetNeighbor(iEdge, jNbr) method
    short GetNeighborEdge(short iEdge, size_t jNbr) const
      { return nedges[i_nbr[iEdge] + jNbr]; }

    // Initializes to dummy values
    Triangle();

    // Destructor - cleans up memory
    ~Triangle();
  };

  /** 
   * A walk is a sequence of triangles that share common edges. All
   * edges crossed in a walk are shared by exactly two triangles, so
   * that the walk is well-defined by the starting and ending triangle
   */
  struct Walk
  {
    // First triangle in the walk and this vertex's index in it
    size_t t0; short v0;

    // Second triangle in the walk and this vertex's index in it
    size_t t1; short v1;

    // The size of the walk (number of triangles) - equal to valence for
    // internal vertices but not for boundary ones!!!
    size_t n;
  };

  /**
   * This structure describes a triangle on the side of an edge. It is 
   * called a 'wing' based on the 'winged edge' concept. However, I don't
   * really remember how winged edges work, so I just made it up on the fly.
   * You can say that I 'winged' it. :)
   */
  struct Wing
    {
    // The triangle index of this wing, or NOID to represent a wing that does
    // not exist.
    size_t t;

    // The indices of the three vertices in this triangle in relation to the 
    // edge in question. 
    short vCenter, vMoving, vOpposite;

    // Default constructor
    Wing() : t(NOID), vCenter(-1), vMoving(-1), vOpposite(-1) {}

    // Constructor
    Wing(size_t in_t, short in_vCenter, short in_vMoving, short in_vOpposite)
      : t(in_t), vCenter(in_vCenter), vMoving(in_vMoving), vOpposite(in_vOpposite) {}
    };

  /**
   * This structure describes a step in a walk around a vertex. It 
   * describes the edge from the fixed vertex to the moving vertex 
   * and the two faces that share it. One of the faces is 'front', 
   * i.e., will be visited next, and the other is 'back', i.e., was
   * just visited.
   */
  struct WalkStep
    {
    // The index of the moving vertex
    size_t iMoving;

    // A step in the walk consists of two wings on the sides of the edge.
    Wing front, back;

    // Constructor
    WalkStep(size_t in_moving, Wing in_front, Wing in_back) 
      : iMoving(in_moving), front(in_front), back(in_back) {}

    // Default constructor
    WalkStep() { iMoving = NOID; }
  };

  /**
   * A walk is simply a chain of walk steps
   */
  typedef std::list<WalkStep> VertexWalk;

  /**
   * Each vertex is associated with one or more walks.
   */
  struct Vertex
    {
    // First triangle in the walk and this vertex's index in it
    std::vector<VertexWalk> walks;
/*     
    // The valence of the the vertex
    size_t valence;

    // Whether the vertex is on the boundary.
    bool bnd;

    // Whether the vertex is on the seam
    bool seam;

    // Is this a boundary vertex
    bool IsBoundary() const { return bnd; }

    // Is this a seam vertex
    bool IsSeam() const { return seam; }

    // Is this an internal vertex
    bool IsInternal() const { return !bnd && !seam; }

    // Get the valence of the vertex. For internal vertices it is equal 
    // to the number of triangles that share the vertex, but for boundary 
    // vertices it is equal to 1 plus that (when the mesh has a 
    // disk topology)
    size_t Valence() const { return valence; }

    // Constructors
    Vertex() : bnd(false), seam(false), valence(0) {};  */
  };

  /** This structure is used to represent edges during VTK import */
  struct ImportEdge : public std::pair<size_t, size_t>
  {
    typedef std::pair<size_t, size_t> Superclass;

    ImportEdge(size_t v1, size_t v2) : 
      Superclass(std::min(v1,v2), std::max(v1,v2)) {};
  };

  struct MeshLevel
  {
    // List of triangles in this mesh level
    vector<Triangle> triangles;

    // List of vertices at the mesh level
    vector<Vertex> vertices;

    // Parent mesh level
    const MeshLevel *parent;

    // Sparse matrix mapping the mesh level to the parent level
    ImmutableSparseMatrix<double> weights;

    // Return the valence of a vertex
    /*
    size_t GetVertexValence(size_t ivtx)
      { return vertices[ivtx].Valence(); } */

    // Check if the vertex is on the boundary
    /*
    bool IsVertexInternal(size_t ivtx)
      { return vertices[ivtx].IsInternal(); } */

    // Set the mesh level as the 'root'. This means that its weight matrix
    // becomes identity and its parent is null
    void SetAsRoot()
      { parent = NULL; weights.Reset(); }

    // Constructor
    MeshLevel()
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
  static void ExportLevelToVTK(MeshLevel &src, vtkPolyData *mesh);

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

  /** Dump mesh structure */
  static void DumpTriangles(MeshLevel &mesh);

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

  // Recursive methods used for building up a walk around a vertex
  static void RecursiveBuildWalkForward(
    MeshLevel *mesh, VertexWalk &walk, vnl_matrix<int> &visited, Wing w);

  static void RecursiveBuildWalkBackward(
    MeshLevel *mesh, VertexWalk &walk, vnl_matrix<int> &visited, Wing w);

  // Compute the walks in a mesh level. This associates each vertex with a
  // triangle, which makes subsequent processing a lot easier
  static void ComputeWalks(MeshLevel *mesh);

};



class BranchingEdgeWalkAroundVertex
{
public:

  /** Constructor: takes a vertex in a fully initialized mesh level */
  BranchingEdgeWalkAroundVertex(
    const BranchingSubdivisionSurface::MeshLevel *mesh, 
    size_t iVertex, size_t iWalk)
    : walk(mesh->vertices[iVertex].walks[iWalk])
    {
    // Store the input
    this->mesh = mesh;
    this->iVertex = iVertex;
    this->iWalk = iWalk;

    // Set the position to zero
    it = walk.begin();
    }

  /** Check if this walk is closed or open (boundary vertex) */
  bool IsOpen()
    { return walk.back().front.t == NOID; }

  /**
   * Are we at the end of the walk? The end is reached when you have made a
   * full circle or moved past the end of the mesh. No operations should be
   * performed once you are at the end of the walk
   */
  bool IsAtEnd() const
    { return it == walk.end(); }

  /** Increment operator */
  BranchingEdgeWalkAroundVertex &operator++ ()
    { it++; return *this; }

  /** Go to the last edge in the walk (next step puts you at end) */
  void GoToLastEdge()
    { it--; }

  /** Go to the first edge in the walk (next step puts you at end) */
  void GoToFirstEdge()
    { it = walk.begin(); }

  /** Get the id of the 'fixed' vertex. This is always the same id */
  size_t FixedVertexId()
    { return iVertex; }

  /** Get the moving vertex id */
  size_t MovingVertexId()
    { return it->iMoving; }

  /** Get the triangle ahead */
  size_t TriangleAhead()
    { return it->front.t; }

  /** Get the triangle behind */
  size_t TriangleBehind()
    { return it->back.t; }

  /** Get the index of the fixed vertex in the triangle ahead */
  short FixedVertexIndexInTriangleAhead()
    { return it->front.vCenter; }

  /** Get the index of the fixed vertex in the triangle behind */
  short FixedVertexIndexInTriangleBehind()
    { return it->back.vCenter; }

  /** Get the index of the moving vertex in the triangle ahead */
  short MovingVertexIndexInTriangleAhead()
    { return it->front.vMoving; }

  /** Get the index of the moving vertex in the triangle behind */
  short MovingVertexIndexInTriangleBehind()
    { return it->back.vMoving; }

  /** Get the index of the face-connected vertex in triangle ahead */
  short OppositeVertexIndexInTriangleAhead()
    { return it->front.vOpposite; }

  /** Get the index of the face-connected vertex in triangle behind */
  short OppositeVertexIndexInTriangleBehind()
    { return it->back.vOpposite; }

  /** Get the id of the face-connected vertex ahead */
  size_t VertexIdAhead()
    {
    // Never call this if the walk is at its end
    assert(it != walk.end());

    // Get the pointer to the next triangle
    WalkIterator itNext = it; itNext++;

    // Get the vertex
    if(itNext == walk.end())
      {
      if(it->front.t == NOID) return NOID;
      else return walk.begin()->iMoving;
      }
    else return itNext->iMoving;
    }

  /** Get the id of the face-connected vertex behind */
  size_t VertexIdBehind()
    {
    if(it == walk.begin())
      {
      if(it->back.t == NOID) return NOID;
      else return walk.back().iMoving;
      }
    
    // Go back one step
    WalkIterator itPrev = it; itPrev--;
    return itPrev->iMoving;
    }

  /** 
   * Get the multiplicity of the egde between the fixed vertex and the moving
   * vertex. Multiplicity is the number of triangles that share the edge
   */
  size_t GetEdgeMultiplicity()
    {
    // One of the triangles must be valid
    assert((it->front.t != NOID) || (it->back.t != NOID));

    if(it->front.t != NOID)
      return 1 + mesh->triangles[it->front.t].n_nbr[it->front.vOpposite];
    else
      return 1 + mesh->triangles[it->back.t].n_nbr[it->back.vOpposite];
    }

private:
  // Some typedefs
  typedef BranchingSubdivisionSurface::VertexWalk WalkType;
  typedef WalkType::const_iterator WalkIterator;

  // The mesh being iterated around
  const BranchingSubdivisionSurface::MeshLevel *mesh;

  // The current state of the walk is defined by these variables
  const WalkType &walk;

  // An iterator into the walk (which is a list)
  WalkIterator it;

  // The index of the vertex walked around by this iterator
  size_t iVertex, iWalk;
};

#ifdef COMMENTOUT

/**
 * A scheme for computing tangent vectors from mesh levels
 */
class LoopTangentScheme
{
public:
  typedef BranchingSubdivisionSurface::MeshLevel MeshLevel;

  /** Constructor does nothing */
  LoopTangentScheme();
  ~LoopTangentScheme();

  /** Pass a mesh to the tangent scheme */
  void SetMeshLevel(MeshLevel *in_level);

  /** Direct access to the weights. This method returns the contribution
   * of a neighbor pointed by a walk iterator to the d-tangent at a vertex */
  double GetNeighborWeight(size_t d, BranchingEdgeWalkAroundVertex &walk)
    {
    return xNbrTangentWeights[d][walk.GetPositionInMeshSparseArray()];
    }

  /** Direct access to the weights. This method returns the contribution
   * of the vertex itself to a d-tangent. This is zero for vertices with 
   * regular valence */
  double GetOwnWeight(size_t d, size_t v)
    {
    return xVtxTangentWeights[d][v];
    }


  /** Compute the directional derivative of a function F over the mesh
   * in direction d (0,1) at vertex v */
  template<class T> T GetPartialDerivative(size_t d, size_t v, const T *F)
  {
    size_t *ri = level->nbr.GetRowIndex();
    size_t *ci = level->nbr.GetColIndex();

    double y = xVtxTangentWeights[d][v] * F[v];
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += xNbrTangentWeights[d][i] * F[ci[i]];

    return y;
  }

protected:

  // Mesh level for which tangents are computed
  MeshLevel *level;

  // Method to delete internally stored data
  void Reset();

  // An array of weights used to compute the tangent vectors at each vertex.
  // These correspond to the sparse matrix A, i.e., include the neighbors of the
  // vertex and the vertex itself
  double *xNbrTangentWeights[2];
  double *xVtxTangentWeights[2];
};

#endif //COMMENTOUT

inline ostream &operator << (ostream &out, const BranchingSubdivisionSurface::Triangle &t)
{
  // Vertices
  out << "[V=(" << t.vertices[0] << "," << t.vertices[1] 
    << "," << t.vertices[2] << ")";

  // Neighbors
  for(size_t i = 0; i < 3; i++)
    {
    out << " N[" << i << "]={";
    for(size_t j = 0; j < t.n_nbr[i]; j++)
      {
      out << "(" << t.nbr[t.i_nbr[i] + j];
      out << "," << t.nedges[t.i_nbr[i] + j];
      out << ")";
      }
    cout << "} ";
    }

  return out;
}

#endif
