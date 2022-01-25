#include "BranchingSubdivisionSurface.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vnl/vnl_vector_fixed.h"
#include <string>
#include <list>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

BranchingSubdivisionSurface::Triangle::Triangle()
{
  for(size_t i = 0; i < 3; i++)
  {
    vertices[i] = NOID;
    n_nbr[i] = 0;
    i_nbr[i] = 0;
  }
  i_nbr[3] = 0;
  nbr = NULL;
  nedges = NULL;
}

BranchingSubdivisionSurface::Triangle::~Triangle()
{
  if(i_nbr[3] > 0)
    {
    delete nbr;
    delete nedges;
    }
}

void
BranchingSubdivisionSurface
::SetOddVertexWeights(MutableSparseMatrix &W,
                      const MeshLevel *parent, MeshLevel *child, 
                      size_t t, size_t v)
{
  // Weight constants
  const static double W_INT_EDGE_CONN = 3.0 / 8.0;
  const static double W_INT_FACE_CONN = 1.0 / 8.0;
  const static double W_BND_EDGE_CONN = 1.0 / 2.0;

  // Get the parent and child triangles
  const Triangle &tp = parent->triangles[t / 4];
  Triangle &tc = child->triangles[t];

  // Get the child vertex index
  size_t ivc = tc.vertices[v];

  // Find out if this is a boundary/crease vertex. Since the child triangle 
  // is a center triangle, we must check if the parent triangle has a 
  // neighbor accross edge v
  if(tp.n_nbr[v] == 1)
  {
    // Get the neighbor of the parent triangle
    const Triangle &topp = parent->triangles[tp.GetNeighbor(v,0)];

    // Internal triangle. It's weights are 3/8 for the edge-connected parent
    // vertices and 1/8 for the face-connected vertices
    W(ivc,tp.vertices[ror(v)]) = W_INT_EDGE_CONN;
    W(ivc,tp.vertices[rol(v)]) = W_INT_EDGE_CONN;
    W(ivc,tp.vertices[v]) = W_INT_FACE_CONN;
    W(ivc,topp.vertices[tp.GetNeighborEdge(v,0)]) = W_INT_FACE_CONN;
  }
  else
  {
    // Only the edge-connected vertices are involved
    W(ivc,tp.vertices[ror(v)]) = W_BND_EDGE_CONN;
    W(ivc,tp.vertices[rol(v)]) = W_BND_EDGE_CONN;
  }
}

inline bool WalkStepPredicate (
  const BranchingSubdivisionSurface::WalkStep &p1,
  const BranchingSubdivisionSurface::WalkStep &p2)
{
  return p1.iMoving < p2.iMoving;
}


void BranchingSubdivisionSurface
::RecursiveBuildWalkForward(
  MeshLevel *mesh, VertexWalk &walk, vnl_matrix<int> &visited, Wing w)
{
  // Mark this triangle as visited
  visited[w.t][w.vCenter] = 1;

  // Get the reference to the current triangle
  Triangle &T = mesh->triangles[w.t];

  // Create a new step structure that will be appended to the walk. 
  // By default, the wings in the step are initialized to NOID
  WalkStep step;

  // Populate information about the back wing (simply by inversion of the 
  // moving and opposite vertex indices in the front wing of the last step)
  step.back.t = w.t;
  step.back.vCenter = w.vCenter;
  step.back.vMoving = w.vOpposite;
  step.back.vOpposite = w.vMoving;

  // Set the iMoving vertex for this step
  step.iMoving = T.vertices[step.back.vMoving];

  // Check if it is possible to continue walking forward (i.e., the edge ahead
  // is an interior edge, there is only one neighbor triangle)
  if(T.n_nbr[w.vMoving] == 1)
    {
    // Fill out the front step data structure
    step.front.t = T.GetNeighbor(w.vMoving, 0);
    step.front.vOpposite = T.GetNeighborEdge(w.vMoving, 0);

    // Get a handle to the triangle ahead
    Triangle &TF = mesh->triangles[step.front.t];

    // Fill out the other two indices by searching for the vertex number
    if(TF.vertices[ror(step.front.vOpposite)] == T.vertices[w.vCenter])
      {
      step.front.vCenter = ror(step.front.vOpposite);
      step.front.vMoving = rol(step.front.vOpposite);
      }
    else
      {
      step.front.vCenter = rol(step.front.vOpposite);
      step.front.vMoving = ror(step.front.vOpposite);
      }
    }

  // Add the step to the walk
  walk.push_back(step);

  // printf("FWK APPENDING iM=%d, front=(%d,%d,%d,%d), back=(%d,%d,%d,%d)\n",
  //   step.iMoving, 
  //   step.front.t, step.front.vCenter, step.front.vMoving, step.front.vOpposite,
  //   step.back.t, step.back.vCenter, step.back.vMoving, step.back.vOpposite);

  // Proceed to the next step if it's possible
  if(step.front.t != NOID && !visited[step.front.t][step.front.vCenter])
    RecursiveBuildWalkForward(mesh, walk, visited, step.front);
}


void BranchingSubdivisionSurface
::RecursiveBuildWalkBackward(
  MeshLevel *mesh, VertexWalk &walk, vnl_matrix<int> &visited, Wing w)
{
  // Mark this triangle as visited
  visited[w.t][w.vCenter] = 1;

  // Get the reference to the current triangle
  Triangle &T = mesh->triangles[w.t];

  // Create a new step structure that will be appended to the walk. 
  // By default, the wings in the step are initialized to NOID
  WalkStep step;

  // Populate information about the back wing (simply by inversion of the 
  // moving and opposite vertex indices in the front wing of the last step)
  step.front.t = w.t;
  step.front.vCenter = w.vCenter;
  step.front.vMoving = w.vOpposite;
  step.front.vOpposite = w.vMoving;

  // Set the iMoving vertex for this step
  step.iMoving = T.vertices[step.front.vMoving];

  // Check if it is possible to continue walking forward (i.e., the edge ahead
  // is an interior edge, there is only one neighbor triangle)
  if(T.n_nbr[w.vMoving] == 1)
    {
    // Fill out the front step data structure
    step.back.t = T.GetNeighbor(w.vMoving, 0);
    step.back.vOpposite = T.GetNeighborEdge(w.vMoving, 0);

    // Get a handle to the triangle ahead
    Triangle &TF = mesh->triangles[step.back.t];

    // Fill out the other two indices by searching for the vertex number
    if(TF.vertices[ror(step.back.vOpposite)] == T.vertices[w.vCenter])
      {
      step.back.vCenter = ror(step.back.vOpposite);
      step.back.vMoving = rol(step.back.vOpposite);
      }
    else
      {
      step.back.vCenter = rol(step.back.vOpposite);
      step.back.vMoving = ror(step.back.vOpposite);
      }
    }

  // Add the step to the walk
  walk.push_front(step);

  // printf("BWK APPENDING iM=%d, front=(%d,%d,%d,%d), back=(%d,%d,%d,%d)\n",
  //   step.iMoving, 
  //   step.front.t, step.front.vCenter, step.front.vMoving, step.front.vOpposite,
  //   step.back.t, step.back.vCenter, step.back.vMoving, step.back.vOpposite);

  // Proceed to the next step if it's possible
  if(step.back.t != NOID && !visited[step.back.t][step.back.vCenter])
    RecursiveBuildWalkBackward(mesh, walk, visited, step.back);
}


void BranchingSubdivisionSurface
::ComputeWalks(MeshLevel *mesh)
{
  size_t t = 0, v = 0, i = 0;

  // Empty the list of walks at each vertex
  for(i = 0; i < mesh->vertices.size(); i++)
    mesh->vertices[i].walks.clear();

  // Create a flag array indicating which vertices in which triangles have
  // already been visited.
  vnl_matrix<int> flagVisit(mesh->triangles.size(), 3);
  flagVisit.fill(0);

  // Loop over all triangles, all vertices
  for(t = 0; t < mesh->triangles.size(); t++) for(v = 0; v < 3; v++)
    {
    // If the triangle/vertex has been visited it can't be in a walk
    if(flagVisit[t][v] > 0) continue; 

    // Initialize a walk that we will fill in
    BranchingSubdivisionSurface::VertexWalk walk;

    // This is the index of the current vertex
    size_t ivtx = mesh->triangles[t].vertices[v];

    // Create wings for the forward and reverse walks
    Wing wForward(t, v, ror(v), rol(v)), wBackward(t, v, rol(v), ror(v));

    // printf("*** WALK AROUND VERTEX %d (t = %d, v = %d) ***\n",
    //   ivtx, t, v);

    // Build up a walk in the forward direction
    RecursiveBuildWalkForward(mesh, walk, flagVisit, wForward);

    // Check if the walk is a loop (last step's front.t is not NOID)
    if(walk.back().front.t != NOID)
      {
      // Sort the walk so it starts with the smallest vertex index
      rotate(walk.begin(),
        min_element(walk.begin(), walk.end(), &WalkStepPredicate),
        walk.end());
      }
    else
      {
      // Perform a backwards walk
      RecursiveBuildWalkBackward(mesh, walk, flagVisit, wBackward);
      }

    // Now, assign this walk to the vertex
    mesh->vertices[ivtx].walks.push_back(walk);

    /*
    // Report the walk
    printf("FINAL WALK AROUND VERTEX %d\n", ivtx);
    VertexWalk::iterator itw = walk.begin();
    for(; itw != walk.end(); ++itw)
      {
      printf("Moving: %d, TFront: %d, VFront %d, TBack: %d, VBack: %d\n",
        itw->iMoving, 
        itw->front.t, itw->front.vMoving, 
        itw->back.t, itw->back.vMoving);
      }
    cout << endl;
    */

    }
}

/* 
void BranchingSubdivisionSurface
::ComputeWalks(MeshLevel *mesh)
{
  size_t t = 0, v = 0, i = 0;

  // Empty the list of walks at each vertex
  for(i = 0; i < mesh->vertices.size(); i++)
    mesh->vertices[i].walks.clear();

  // Create a flag array indicating which vertices in which triangles have
  // already been visited.
  vnl_matrix<int> flagVisit(mesh->triangles.size(), 3, 0);

  // Loop over all triangles, all vertices
  for(t = 0; t < mesh->triangles.size(); t++) for(v = 0; v < 3; v++)
    {
    // If the triangle/vertex has been visited it can't be in a walk
    if(flagVisit[t][v] > 0) continue; 

    // Initialize a walk that we will fill in
    BranchingSubdivisionSurface::VertexWalk walk;

    // This is the index of the current vertex
    size_t ivtx = mesh->triangles[t].vertices[v];

    // The current position in the walk
    size_t tWalk = t; short vWalk = v;

    // The position at which the walk will loop around
    size_t tLoop = t;

    // Walk until reaching a non-internal edge
    printf("WALK AROUND vertex %d in T[%d], whose index is %d\n",
      v, t, ivtx);

    do
      {
      // Get a reference to the current triangle
      Triangle &T = mesh->triangles[tWalk];

      // Update the 'visited' flags
      flagVisit[tWalk][vWalk] = 1;

      // Check if the next triangle is walkable or not
      size_t tNext; short vNext;
      if(T.n_nbr[ror(vWalk)] != 1)
        {
        tNext = NOID;
        vNext = -1;
        }
      else
        {
        tNext = T.GetNeighbor(ror(vWalk), 0);
        vNext = ror(T.GetNeighborEdge(ror(vWalk), 0));
        }

      // Put the current edge in the triangle into the list
      cout << "Fwd walk: visiting vtx. " << T.vertices[rol(vWalk)] <<
          "; behind: (" << tWalk << "," << vWalk <<
          "); ahead: (" << tNext << "," << vNext << ")" <<  endl;
      walk.push_back(
        WalkStep(T.vertices[rol(vWalk)], tNext, vNext, tWalk, vWalk));


      // walks[ivtx].push_back( make_pair(
      //     T.vertices[rol(vWalk)],
      //     NeighborInfo(tNext, vNext, tWalk, vWalk)));

      // Update the position in the walk
      tWalk = tNext; vWalk = vNext;
      }
    while(tWalk != NOID && tWalk != tLoop);

    // Now, if we hit a NOID, we can need to walk in the opposite direction,
    // starting from the same triangle as before
    if(tWalk == NOID)
      {
      // Reset the position in the walk
      tWalk = t; vWalk = v;

      // Walk again, in a symmetrical loop
      do
        {
        // Get a reference to the current triangle
        Triangle &T = mesh->triangles[tWalk];

        // Update the 'visited' flags
        flagVisit[tWalk][vWalk] = 1;

        // Check if the next triangle is walkable or not
        size_t tNext; short vNext;
        if(T.n_nbr[rol(vWalk)] != 1)
          {
          tNext = NOID;
          vNext = -1;
          }
        else
          {
          tNext = T.GetNeighbor(rol(vWalk),0);
          vNext = rol(T.GetNeighborEdge(rol(vWalk),0));
          }

        // Put the current edge in the triangle into the list
        cout << "Rev walk: visiting vtx. " << T.vertices[ror(vWalk)] <<
           "; behind: (" << tNext << "," << vNext <<
           "); ahead: (" << tWalk << "," << vWalk << ")" << endl;
        walk.push_front(
          WalkStep(T.vertices[ror(vWalk)], tWalk, vWalk, tNext, vNext));

        // walks[ivtx].push_front( make_pair(
        //    T.vertices[ror(vWalk)],
        //    NeighborInfo(tWalk, vWalk, tNext, vNext)));

        // Update the position in the walk
        tWalk = tNext; vWalk = vNext;
        }
      while(tWalk != NOID);
      }
    else // i.e., the walk loops around
      {
      // Rotate the walk so that the smallest member is the first element
      // (for consistency with MMA code)
      rotate(walk.begin(),
        min_element(walk.begin(), walk.end(), &WalkStepPredicate),
        walk.end());
      }

    // Now, assign this walk to the vertex
    mesh->vertices[ivtx].walks.push_back(walk);

    // Report the walk
    printf("FINAL WALK AROUND VERTEX %d\n", ivtx);
    VertexWalk::iterator itw = walk.begin();
    for(; itw != walk.end(); ++itw)
      {
      printf("Moving: %d, TFront: %d, VFront %d, TBack: %d, VBack: %d\n",
        itw->iMoving, itw->tFront, itw->vFront, itw->tBack, itw->vBack);
      }
    cout << endl;
    }

  cout << "Walk Done\n";
}
*/

/**
 * This method should be called for child-level triangles whose index in the
 * parent triangle matches the vertex id (v).
 */
void BranchingSubdivisionSurface
::SetEvenVertexWeights(MutableSparseMatrix &W,
                       const MeshLevel *parent, MeshLevel *child, size_t t, size_t v)
{
  // Weight constants
  const static double W_BND_EDGE_CONN = 1.0 / 8.0;
  const static double W_BND_SELF = 3.0 / 4.0;

  // Get the child vertex index
  size_t ivc = child->triangles[t].vertices[v];

  // Get the corresponding parent triangle and vertex index
  size_t tp = t / 4;
  size_t ivp = parent->triangles[tp].vertices[v];

  // Find all edge-connected, crease-connected and sheet-connected vertices
  std::set<size_t> cEdge, cSeam, cSheet;
  typedef std::set<size_t>::iterator SetIt;

  // printf("*** SETTING WEIGHTS FOR EVEN VERTEX %d ***\n", ivc);

  // Go through all the walks for this vertex in the parent mesh, in order to
  // collect edge-connected, seam-connected and sheet-connected neighbors.
  for(size_t i = 0; i < parent->vertices[ivp].walks.size(); i++)
    {
    for(BranchingEdgeWalkAroundVertex it(parent, ivp, i); !it.IsAtEnd(); ++it)
      {
      size_t ivn = it.MovingVertexId();
      size_t mult = it.GetEdgeMultiplicity();
      if(mult == 1)
        {
        // printf("Adding boundary-edge (%d %d)\n", ivc, ivn);
        cEdge.insert(ivn);
        }
      else if(mult == 2)
        {
        // printf("Adding sheet-edge (%d %d)\n", ivc, ivn);
        cSheet.insert(ivn);
        }
      else if(mult == 3)
        {
        // printf("Adding seam-edge (%d %d)\n", ivc, ivn);
        cSeam.insert(ivn);
        }
      else
        throw "Wrong multiplicity for an edge";
      }    
    }

  // The vertex is internal
  if(cEdge.size() == 0 && cSeam.size() == 0)
    {
    // Get the number of neighbor vertices
    int n = cSheet.size();

    // Compute the beta constant
    double beta = (n > 3) ? 3.0 / (8.0 * n) : 3.0 / 16.0;

    // Assign beta to each of the edge-adjacent vertices
    for(SetIt it = cSheet.begin(); it != cSheet.end(); ++it)
      W(ivc, *it) = beta;

    // Assign the balance to the coincident vertex
    W(ivc, ivp) = 1.0 - n * beta;
    }

  // The vertex lies on the edge
  else if(cEdge.size() == 2 && cSeam.size() == 0)
    {
    // Add the point itself
    W(ivc, ivp) = W_BND_SELF;

    // Add the starting point
    for(SetIt it = cEdge.begin(); it != cEdge.end(); ++it)
      W(ivc, *it) = W_BND_EDGE_CONN;
    }

  // The vertex is on the seam
  else if(cEdge.size() == 0 && cSeam.size() == 2)
    {
    // Add the point itself
    W(ivc, ivp) = W_BND_SELF;

    // Add the starting point
    for(SetIt it = cSeam.begin(); it != cSeam.end(); ++it)
      W(ivc, *it) = W_BND_EDGE_CONN;
    }

  // The vertex is an edge-seam vertex
  else if(cSeam.size() == 1 && cEdge.size() == 1)
    {
    // Treat the seam vertex as an internal vertex, and ignore the edge vertex
    int n = cSheet.size() + cSeam.size();

    // Compute the beta constant
    double beta = (n > 3) ? 3.0 / (8.0 * n) : 3.0 / 16.0;

    // Assign beta to each of the edge-adjacent vertices
    for(SetIt it = cSheet.begin(); it != cSheet.end(); ++it)
      W(ivc, *it) = beta;

    // Assign beta to the seam vertex
    W(ivc, *cSeam.begin()) = beta;

    // Assign the balance to the coincident vertex
    W(ivc, ivp) = 1.0 - n * beta;
    }

  // The remaining possibilities are 'weird'. We keep the vertex in its place
  else
    {
    W(ivc, ivp) = 1.0;
    }
}

void
BranchingSubdivisionSurface
::RecursiveAssignVertexLabel(
  MeshLevel *mesh, size_t t, size_t v, size_t id)
{
  Triangle &T = mesh->triangles[t];

  // printf("RecursiveAssignVertexLabel: t = %d, v = %d, id = %d\n",t,v,id);

  // If the vertex has a label assigned, leave it alone
  if(T.vertices[v] != NOID)
    return;
    
  // Assign the label to the vertex
  mesh->triangles[t].vertices[v] = id;

  // printf("ASSIGNED!\n");

  // Visit the vertex in all neighboring triangles across both edges
  for(size_t q = 1; q <= 2; q++)
    {
    // The edge across which we want to go to 
    short w = (v+q) % 3;

    // Repeat for all adjacent triangles across that edge
    for(size_t i = 0; i < T.n_nbr[w]; i++)
      RecursiveAssignVertexLabel(mesh, 
        T.GetNeighbor(w, i), 
        (T.GetNeighborEdge(w, i) + q) % 3, id);
    }
}

void BranchingSubdivisionSurface::Subdivide(const MeshLevel *parent, MeshLevel *child)
{
  // Get the numbers of triangles before and after
  size_t ntParent = parent->triangles.size();
  size_t ntChild = 4 * ntParent;
  size_t i, j, k;

  // Connect the child to the parent
  child->parent = parent;

  // Reset the child mesh
  child->vertices.clear();
  child->triangles.clear();

  // Create a list of unique edges in the parent mesh. Since each edge is
  // subdivided, the edges can be used to assign indices to newly inserted
  // vertices.
  map<ImportEdge, size_t> edgemap;
  for(i = 0; i < ntParent; i++) for(j = 0; j < 3; j++)
    {
    // Create an edge representation
    ImportEdge ie(
      parent->triangles[i].vertices[j],
      parent->triangles[i].vertices[ror(j)]);

    // Check if the edge is already present and if it's not assign it an
    // index and store in the map
    if(edgemap.find(ie) == edgemap.end())
      {
      // printf("Odd vertex %d is on edge (%d %d)\n",
      //   edgemap.size() + parent->vertices.size(), 
      //   ie.first, ie.second);
      edgemap.insert(make_pair(ie, edgemap.size()));
      }
    }

  // Allocate the triangles in the child mesh
  child->triangles.resize(ntChild);

  // Allocate the vertices. The numbering scheme is such that a child's 
  // even vertex is the same as the number of the corresponding vertex in
  // the parent.
  child->vertices.resize(parent->vertices.size() + edgemap.size());

  // Subdivide each triangle into four
  for (i = 0; i < ntParent; i++)
    {
    // Get pointers to the four children
    const Triangle &pt = parent->triangles[i];
    Triangle *ct[4] = {
                        &child->triangles[4*i], &child->triangles[4*i+1],
                        &child->triangles[4*i+2], &child->triangles[4*i+3]};

    // Get the indices in the child mesh for the three odd vetices that will
    // be inserted opposite to the three vertices in the parent triangle
    size_t iodd[3];
    for (j = 0; j < 3; j++)
      {
      ImportEdge ie(pt.vertices[ror(j)], pt.vertices[rol(j)]);
      iodd[j] = parent->vertices.size() + edgemap[ie];
      }

    // Set the neighbor information for the three outer children of this
    // parent triangle
    for (j = 0; j < 3; j++)
      {
      // Assign a vertex number to the even vertex in the child
      ct[j]->vertices[j] = pt.vertices[j];

      // Assign vertex numbers to the odd vertices in the child
      ct[j]->vertices[ror(j)] = iodd[rol(j)];
      ct[j]->vertices[rol(j)] = iodd[ror(j)];

      // Count the number of neighbors and allocate data structures
      size_t n_total = pt.n_nbr[ror(j)] + pt.n_nbr[rol(j)] + 1;
      ct[j]->nbr = new size_t[n_total];
      ct[j]->nedges = new short[n_total];

      // Compute the neighbors across each edge
      size_t ipos = 0;
      for (k = 0; k < 3; k++)
        {
        if(k == j)
          {
          ct[j]->n_nbr[k] = 1;
          ct[j]->nbr[ipos] = 4*i + 3;
          ct[j]->nedges[ipos] = j;
          ipos++;
          }
        else
          {
          ct[j]->n_nbr[k] = pt.n_nbr[k];
          for (size_t a = 0; a < ct[j]->n_nbr[k]; a++)
            {
            // This is the parent's neighbor across the edge 'k'
            size_t nbr_p = pt.GetNeighbor(k, a); 

            // This is the index of the vertex opposite edge 'k' in nbr_p
            short q = pt.GetNeighborEdge(k, a);

            // Get a handle to the neighboring triangle
            const Triangle &ptn = parent->triangles[nbr_p];

            // Figure out what is the index of pt's vertex j in triangle ptn
            short idx_c = (ptn.vertices[ror(q)] == pt.vertices[j]) ? ror(q) : rol(q);

            // Set the actual neigbor index
            ct[j]->nbr[ipos] = nbr_p * 4 + idx_c;
            ct[j]->nedges[ipos] = q;
            ipos++;
            }
          }
        ct[j]->i_nbr[k+1] = ipos;
        }
      }

    // Now set the vertex / neighbor info for the internal triangle
    ct[3]->nbr = new size_t[3];
    ct[3]->nedges = new short[3];
    for (k = 0; k < 3; k++)
      {
      // Set neighbor info
      ct[3]->n_nbr[k] = 1;
      ct[3]->i_nbr[k+1] = k+1;
      ct[3]->nbr[k] = 4 * i + k;
      ct[3]->nedges[k] = k;

      // Set vertex info
      ct[3]->vertices[k] = iodd[k];
      }
    }

  // Compute the weights of the child vertices wrt parent vertices
  MutableSparseMatrix W(child->vertices.size(), parent->vertices.size());

  // Assign weights to the even and odd vertices
  for(i = 0; i < ntParent; i++) for(j = 0; j < 3; j++)
    {
    if(W.empty_row(child->triangles[4 * i + j].vertices[j]))
      SetEvenVertexWeights(W, parent, child, 4*i+j, j);
    if(W.empty_row(child->triangles[4 * i + 3].vertices[j]))
      SetOddVertexWeights(W, parent, child, 4*i+3, j);
    }

  // Copy the sparse matrix into immutable form
  child->weights.SetFromVNL(W);

  // If the parent's parent is not NULL, we need to multiply the sparse
  // matrices of the parent and child
  if(parent->parent)
    ImmutableSparseMatrix<double>::Multiply(
      child->weights, child->weights, parent->weights);

  // Compute the walks in the child mesh
  ComputeWalks(child);
}

void
BranchingSubdivisionSurface::
RecursiveSubdivide(const MeshLevel *src, MeshLevel *dst, size_t n)
{
  // N has to be greater than 1
  if(n == 0)
  { *dst = *src; return; }
  else if(n == 1)
  { Subdivide(src, dst); return; }

  // Create an array of intermedial mesh levels
  MeshLevel *temp = new MeshLevel[n - 1];

  // Subdivide the intermediate levels
  const MeshLevel *parent = src;
  for(size_t i = 0; i < n-1; i++)
  {
    Subdivide(parent, temp + i);
    parent = temp + i;
  }

  // Subdivide the last level
  Subdivide(parent, dst);

  // Delete the intermediates
  delete[] temp;

  // Set the parent pointer in dst to src (bypass intermediates)
  dst->parent = src;
}

void BranchingSubdivisionSurface::ExportLevelToVTK(MeshLevel &src, vtkPolyData *mesh)
{
  // The important thing with importing and exporting mesh levels is that the
  // order of triangles and vertices must be maintained.

  // Allocate an array of cells
  mesh->Allocate(src.triangles.size());

  // Add all of the triangles in the mesh
  for(size_t i = 0; i < src.triangles.size(); i++)
  {
    vtkIdType pts[3];
    pts[0] = src.triangles[i].vertices[0];
    pts[1] = src.triangles[i].vertices[1];
    pts[2] = src.triangles[i].vertices[2];
    mesh->InsertNextCell(VTK_TRIANGLE, 3, pts);
  }
}

void BranchingSubdivisionSurface::ImportLevelFromVTK(vtkPolyData *mesh, MeshLevel &dest)
{
  size_t i;
  short j;

  // Prepare the mesh
  mesh->BuildCells();

  // A representation of a triangle at an edge (triangle index, edge 
  // index in triangle)
  typedef std::pair<size_t, short> TRep;

  // A container associating edges with triangles
  typedef std::multimap<ImportEdge, TRep> EdgeMap;

  // An edge map for accumulating neighborhood info
  EdgeMap edges;

  // Get the number of triangles in the mesh
  size_t nTriangles = mesh->GetNumberOfCells();
  dest.triangles.resize(nTriangles);

  // Set the number of vertices in the mesh
  size_t nVertices = mesh->GetNumberOfPoints();
  dest.vertices.resize(nVertices);

  // For each triangle, compute the neighbor. This can be done by first enumerating
  // all the edges in the mesh. For each edge there will be one or two triangles
  for(i = 0; i < nTriangles; i++)
  {
    // Get the points from the current triangle
    vtkIdType npts;
    const vtkIdType *pts;
    mesh->GetCellPoints(i, npts, pts);

    // If the number of points is not 3, return with an exception
    if(npts != 3) throw string("Mesh contains cells other than triangles");

    // Associate each edge with a triangle
    for(j = 0; j < 3; j++)
    {
      // Set the vertices in each triangle
      dest.triangles[i].vertices[j] = pts[j];

      // Create a structure representing the edge opposite to vertex j
      ImportEdge edge(pts[ror(j)], pts[rol(j)]);

      // Associate the edge with a t-rep
      edges.insert(make_pair(edge, TRep(i, j)));
    }
  }

  // Take a second pass through the triangles, this time assign all the neighbor info
  for(i = 0; i < nTriangles; i++)
    {
    // Get a shorthand for the triangle
    Triangle &T = dest.triangles[i];

    // First, how many total neighbors are there
    size_t ntotal = 0;

    // The first value in neighbor array index is always 0
    T.i_nbr[0] = 0;

    // First pass sets up array indices 
    for(j = 0; j < 3; j++)
      {
      // Get the current edge
      ImportEdge edge(T.vertices[ror(j)], T.vertices[rol(j)]);

      // Set the number of neighbors
      T.n_nbr[j] = edges.count(edge) - 1; 

      // Increment total number of neighbors
      ntotal += T.n_nbr[j];

      // Set the neighbor array index
      T.i_nbr[j+1] = ntotal;
      }

    // Allocate the nedge, nbr structures
    T.nbr = new size_t[ntotal];
    T.nedges = new short[ntotal];

    // Now, populate the arrays
    for(j = 0; j < 3; j++)
      {
      // Get the current edge
      ImportEdge edge(T.vertices[ror(j)], T.vertices[rol(j)]);

      // Loop over all triangles associated with that edge
      size_t a = 0;
      EdgeMap::iterator it;
      for(it = edges.lower_bound(edge); it != edges.upper_bound(edge); ++it)
        {
        // Check if it's the same triangle
        if(it->second.first == i) continue;
        
        // Set the neighbor and index values
        T.nbr[ T.i_nbr[j] + a ] = it->second.first;
        T.nedges[ T.i_nbr[j] + a ] = it->second.second;
        a++;
        }
      }
    }
  
  // Set the mesh's parent to NULL
  dest.parent = NULL;

  // Finally, compute the walks in this mesh
  ComputeWalks(&dest);
}

void BranchingSubdivisionSurface
::DumpTriangles(MeshLevel &m)
{
  size_t i,a;
  short j;

  for(i = 0; i < m.triangles.size(); i++)
    {
    Triangle &t = m.triangles[i];
    for(j = 0; j < 3; j++)
      {
      printf("T[%d]\'s vertex %d is %d\n", (int) i, j, (int) t.vertices[j]);
      printf("T[%d] has %d neighbors across edge %d\n", (int) i, (int) t.n_nbr[j], j);
      for(a = 0; a < t.n_nbr[j]; a++)
        {
        printf("T[%d]'s neighbor %d across edge %d is T[%d] (reciprocal edge %d)\n",
          (int) i, (int) a, j, (int) t.GetNeighbor(j, a), t.GetNeighborEdge(j, a));
        }
      }
    }
}

void BranchingSubdivisionSurface
::ApplySubdivision(vtkPolyData *src, vtkPolyData *target, MeshLevel &m)
{
  size_t i;

  // Initialize the target mesh
  vtkPoints *tpoints = vtkPoints::New();
  vtkCellArray *tcells = vtkCellArray::New();

  // Add the points to the output mesh
  for(i = 0; i < m.vertices.size(); i++)
  {
    // Output point
    vnl_vector_fixed<double, 3> p(0.0);

    // Iterate over columns
    typedef ImmutableSparseMatrix<double>::RowIterator IteratorType;
    for(IteratorType it = m.weights.Row(i); !it.IsAtEnd(); ++it)
    {
      vnl_vector_fixed<double, 3> x(src->GetPoints()->GetPoint(it.Column()));
      p += x * (double) it.Value();
    }

    // Add the point to the output mesh
    tpoints->InsertNextPoint(p.data_block());
  }

  // Add the cells to the mesh
  for(i = 0; i < m.triangles.size(); i++)
  {
    Triangle &t = m.triangles[i];
    vtkIdType ids[3];
    ids[0] = t.vertices[0];
    ids[1] = t.vertices[1];
    ids[2] = t.vertices[2];
    tcells->InsertNextCell(3, ids);
  }

  // Set the points in the mesh
  target->SetPoints(tpoints);
  target->SetPolys(tcells);
  tpoints->Delete();
  tcells->Delete();
}

void BranchingSubdivisionSurface
::ApplySubdivision(SMLVec3d *xSrc, double *rhoSrc, SMLVec3d *xTrg, double *rhoTrg, MeshLevel &m)
{
  size_t i;

  // Add the points to the output mesh
  for(i = 0; i < m.vertices.size(); i++)
  {
    // Output point
    xTrg[i].fill(0.0); rhoTrg[i] = 0;

    // Iterate over columns
    typedef ImmutableSparseMatrix<double>::RowIterator IteratorType;
    for(IteratorType it = m.weights.Row(i); !it.IsAtEnd(); ++it)
    {
      xTrg[i] += it.Value() * xSrc[it.Column()];
      rhoTrg[i] += it.Value() * rhoSrc[it.Column()];
    }
  }
}

void BranchingSubdivisionSurface
::ApplySubdivision(const double *xsrc, double *xdst, size_t nComp, MeshLevel &m)
{
  typedef ImmutableSparseMatrix<double>::RowIterator IteratorType;

  // Get the total number of values
  size_t n = nComp * m.vertices.size();

  // Set the target array to zero
  std::fill_n(xdst, n, 0.0);

  // Iterate over the output vertices
  for(size_t iv = 0; iv < m.vertices.size(); iv++)
  {
    size_t i = iv * nComp;

    // Iterate over contributing input vertices
    for(IteratorType it = m.weights.Row(iv); !it.IsAtEnd(); ++it)
    {
      double w = it.Value();
      size_t j = it.Column() * nComp;

      // Iterate over the component
      for(size_t k = 0; k < nComp; k++)
        xdst[i+k] += w * xsrc[j+k];
    }
  }
}

bool BranchingSubdivisionSurface::CheckMeshLevel (MeshLevel &mesh)
{
  // Check the following rules for all triangles in the mesh
  // 1. I am one of the neighbors of each of my neighbor triangles
  // 2. The indexes of an edge shared by two neighboring triangles are correct
  size_t nerr = 0;
  for(size_t i = 0; i < mesh.triangles.size(); i++)
    {
    // This is the triangle we want to test
    Triangle &t = mesh.triangles[i];

    // Loop over its vertices
    for(size_t j = 0; j < 3; j++)
      {
      // Repeat for each of the neighbor triangles
      for(size_t a = 0; a < t.n_nbr[j]; a++)
        {
        // Get the neighbor triangle index
        size_t in = t.GetNeighbor(j, a);

        // Get the index of shared edge with t in the neighbor triangle
        short jn = t.GetNeighborEdge(j, a);

        // Get the shorthand to the neighbor triangle
        Triangle &tn = mesh.triangles[in];

        // Check the number of edges across the shared edge is correct
        if(tn.n_nbr[jn] != t.n_nbr[j])
          {
          printf("Error %d: T[%d].n_nbr[%d] != T[%d].n_nbr[%d]\n",
            (int) nerr++, (int) i, (int) j, (int) in, (int) jn); 
          }

        // Check that one of the neighbors of tn across jn is t
        size_t tCount = 0;
        for(size_t an = 0; an < tn.n_nbr[jn]; an++)
          {
          if(tn.GetNeighbor(jn, an) == i)
            {
            // Record that t has been found
            tCount++;

            // Check that the index matches
            if((int) j != tn.GetNeighborEdge(jn, an))
              {
              printf("Error %d: T[%d].nedge[%d,%d] != T[%d].nedge[%d,%d]\n",
                (int) nerr++, (int) i, (int) j, (int) a, (int) in, (int) jn, (int) an); 
              }
            }
          }

        // Check that the count is 1
        if(tCount == 0)
          {
          printf("Error %d: T[%d] is not a neighbor of its neighbor T[%d]\n",
            (int) nerr++, (int) i, (int) in); 
          }
        else if(tCount > 1)
          {
          printf("Error %d: T[%d] appears %d times as a neighbor of its neighbor T[%d]\n",
            (int) nerr++, (int) i, (int) tCount, (int) in); 
          }
        }
      }
    }

  // Check 
         
/* (old rule 3)
      for(size_t k = 1; k < 3; k++)
      {
        if(t.neighbors[(j+k) % 3] != NOID)
        {
          Triangle &tk = mesh.triangles[t.neighbors[(j+k) % 3]];
          if(t.vertices[j] != tk.vertices[(t.nedges[(j+k) % 3] + k ) % 3])
          {
            cout << "Error " << nerr++ <<
            " Rule 3 violated for i = " << i << ", j = " << j << " and k = " << k << endl;
          }
        }
      } */

  return (nerr == 0);
}


/******************************************************************
 CODE FOR THE LOOP TANGENT SCHEME
 *****************************************************************/

#ifdef COMMENTOUT

LoopTangentScheme::LoopTangentScheme()
{
  level = NULL;
}

LoopTangentScheme::~LoopTangentScheme()
{
  Reset();
}

void LoopTangentScheme::Reset()
{
  if(level)
  {
    level = NULL;
    delete xNbrTangentWeights[0];
    delete xNbrTangentWeights[1];
    delete xVtxTangentWeights[0];
    delete xVtxTangentWeights[1];
  }
}

void LoopTangentScheme::SetMeshLevel(MeshLevel *in_level)
{
  // Clean up if necessary
  Reset();

  // Store the pointer to the mesh level
  level = in_level;

  // Initialize the tangent weights. These arrays assign a weight to all the
  // neighboring vertices. But sometimes the vertex's own weight also needs
  // to be included. How to handle this?
  size_t nSparseEntries = level->nbr.GetNumberOfSparseValues();

  for(size_t d = 0; d < 2; d++)
  {
    xNbrTangentWeights[d] = new double[nSparseEntries];
    xVtxTangentWeights[d] = new double[level->nVertices];
    fill_n(xNbrTangentWeights[d], nSparseEntries, 0.0);
    fill_n(xVtxTangentWeights[d], level->nVertices, 0.0);
  }

  // Compute tangent weight at each vertex
  for(size_t i = 0; i < level->nVertices; i++)
  {
    // Create an iterator around the vertex
    BranchingEdgeWalkAroundVertex it(level, i);

    // Get the valence of the vertex
    size_t n = it.Valence();

    // If the vertex is internal, weights are easy to compute
    if(!it.IsOpen())
    {
      for(size_t j = 0; !it.IsAtEnd(); ++it, ++j)
      {
        // size_t idxA = xMapVertexNbrToA[it.GetPositionInMeshSparseArray()];
        size_t k = it.GetPositionInMeshSparseArray();
        xNbrTangentWeights[0][k] = cos(j * vnl_math::pi * 2.0 / n);
        xNbrTangentWeights[1][k] = sin(j * vnl_math::pi * 2.0 / n);
      }
    }
    else
    {
      // Get the index of the first neighbor
      size_t kFirst = it.GetPositionInMeshSparseArray();

      // Move the index to the last neighbor
      it.GoToLastEdge();
      size_t kLast = it.GetPositionInMeshSparseArray();

      // We can now set the along-the-edge tangent vector
      xNbrTangentWeights[1][kFirst] = 1.0;
      xNbrTangentWeights[1][kLast] = -1.0;

      // Now we branch according to the valence
      if(n == 2)
      {
        xNbrTangentWeights[0][kFirst] = -1.0;
        xNbrTangentWeights[0][kLast] = -1.0;
        xVtxTangentWeights[0][i] = 2.0;
      }
      else if (n == 3)
      {
        // Move iterator to the previous (middle) vertex
        it.GoToFirstEdge(); ++it;

        // Get the corresponding A index
        size_t kMiddle = it.GetPositionInMeshSparseArray();

        // Set the weights
        xNbrTangentWeights[0][kMiddle] = -1.0;
        xVtxTangentWeights[0][i] = 1.0;
      }
      else
      {
        // Compute the angle theta
        double theta = vnl_math::pi / (n - 1);
        xNbrTangentWeights[0][kFirst] = xNbrTangentWeights[0][kLast] = sin(theta);

        // Assign the weights to intermediate vertices
        it.GoToFirstEdge(); ++it;
        for(size_t j = 1; j < n - 1; ++j, ++it)
        {
          size_t kInt = it.GetPositionInMeshSparseArray();
          xNbrTangentWeights[0][kInt] = (2.0 * cos(theta) - 2.0) * sin(theta * j);
        }
      }
    }
  }
}

#endif // COMMENTOUT

