#include "SubdivisionSurface.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vnl/vnl_vector_fixed.h"
#include <string>
#include <list>
#include <map>
#include <algorithm>

using namespace std;

void
SubdivisionSurface
::SetOddVertexWeights(MutableSparseMatrix &W, size_t ivc,
                      const MeshLevel *parent, size_t t, size_t v,
                      bool flat_mode)
{
  // Weight constants
  const static double W_INT_EDGE_CONN = 3.0 / 8.0;
  const static double W_INT_FACE_CONN = 1.0 / 8.0;
  const static double W_BND_EDGE_CONN = 1.0 / 2.0;

  // Get the parent and child triangles
  const Triangle &tp = parent->triangles[t];

  // Find out if this is a boundary vertex. Since the child triangle is a center
  // triangle, we must check if the parent triangle has a neighbor accross edge v
  if(tp.neighbors[v] != NOID && !flat_mode)
  {
    // Get the neighbor of the parent triangle
    const Triangle &topp = parent->triangles[tp.neighbors[v]];

    // Internal triangle. It's weights are 3/8 for the edge-connected parent
    // vertices and 1/8 for the face-connected vertices
    W(ivc,tp.vertices[(v+1)%3]) = W_INT_EDGE_CONN;
    W(ivc,tp.vertices[(v+2)%3]) = W_INT_EDGE_CONN;
    W(ivc,tp.vertices[v]) = W_INT_FACE_CONN;
    W(ivc,topp.vertices[tp.nedges[v]]) = W_INT_FACE_CONN;
  }
  else
  {
    // Only the edge-connected vertices are involved
    W(ivc,tp.vertices[(v+1)%3]) = W_BND_EDGE_CONN;
    W(ivc,tp.vertices[(v+2)%3]) = W_BND_EDGE_CONN;
  }
}


/**
 * This method should be called for child-level triangles whose index in the
 * parent triangle matches the vertex id (v).
 */
void SubdivisionSurface
::SetEvenVertexWeights(MutableSparseMatrix &W, size_t ivc, 
                       const MeshLevel *parent, size_t ivp,
                       bool flat_mode)
{
  // Weight constants
  const static double W_BND_EDGE_CONN = 1.0 / 8.0;
  const static double W_BND_SELF = 3.0 / 4.0;

  // Handle flat mode first
  if(flat_mode)
    {
    W(ivc, ivp) = 1.0;
    return;
    }

  // Get the vertex object for this vertex
  EdgeWalkAroundVertex it(parent, ivp);
  if(it.IsOpen())
  {
    // Add the point itself
    W(ivc, ivp) = W_BND_SELF;

    // Add the starting point
    W(ivc, it.MovingVertexId()) = W_BND_EDGE_CONN;

    // Get the ending point
    it.GoToLastEdge();
    W(ivc, it.MovingVertexId()) = W_BND_EDGE_CONN;
  }
  else
  {
    // Compute the beta constant
    int n = it.Valence();
    double beta = (n > 3) ? 3.0 / (8.0 * n) : 3.0 / 16.0;

    // Assign beta to each of the edge-adjacent vertices
    for( ; !it.IsAtEnd(); ++it)
      W(ivc, it.MovingVertexId()) = beta;

    // Assign the balance to the coincident vertex
    W(ivc, ivp) = 1.0 - n * beta;
  }
}

void
SubdivisionSurface
::RecursiveAssignVertexLabel(
  MeshLevel *mesh, size_t t, size_t v, size_t id)
{
  if(mesh->triangles[t].vertices[v] != NOID)
    return;
  else
    mesh->triangles[t].vertices[v] = id;
  if(mesh->triangles[t].neighbors[(v+1) % 3] != NOID)
    RecursiveAssignVertexLabel(mesh, mesh->triangles[t].neighbors[(v+1) % 3],
                               (mesh->triangles[t].nedges[(v+1) % 3] + 1) % 3, id);
  if(mesh->triangles[t].neighbors[(v+2) % 3] != NOID)
    RecursiveAssignVertexLabel(mesh, mesh->triangles[t].neighbors[(v+2) % 3],
                               (mesh->triangles[t].nedges[(v+2) % 3] + 2) % 3, id);
}


void SubdivisionSurface::SubdivideSelected(
  const MeshLevel *parent, MeshLevel *child, std::set<size_t> tsel, bool flat_mode)
{
  cout << "Sub with " << tsel.size() << " selected triangles" << endl;

  // Get the numbers of triangles before and after
  size_t ntParent = parent->triangles.size();
  size_t nvParent = parent->nVertices;
  size_t nvChild = nvParent;
  size_t i, j;

  // Parent triangles pointer
  const std::vector<Triangle> &ptl = parent->triangles; 

  // Array that keeps track of new vertices inserted into existing edges
  // (or NOID to indicate the edge is not being split)
  std::vector< std::vector<size_t> > split(ntParent, std::vector<size_t>(3, NOID));

  // Create a queue of triangles, put them all on the queue
  std::list<size_t> tqueue;
  for(i = 0; i < ntParent; i++) tqueue.push_back(i);

  // Split edges if (1) triangle is in tsel; (2) there are already 2 split edges
  while(tqueue.size() > 0)
    {
    size_t t = tqueue.front(); tqueue.pop_front();
    bool mark = (tsel.find(t) != tsel.end());
    for(j = 0; j < 3; j++)
      {
      if(split[t][j] == NOID)
        {
        if(mark || (split[t][(j+1) % 3] != NOID && split[t][(j+2) % 3] != NOID))
          {
          split[t][j] = nvChild;
          if(ptl[t].neighbors[j] != NOID)
            {
            split[ptl[t].neighbors[j]][ptl[t].nedges[j]] = nvChild;
            tqueue.push_back(ptl[t].neighbors[j]);
            }
          nvChild++;
          }
        }
      }
    }

  // Create a mesh generator, to which we will add all child triangles
  TriangleMeshGenerator tmg(child, nvChild);

  // At this point every triangle will have 0, 1 or 3 split edges. This means the triangle
  // will also have 1, 2 or 4 children. Add these triangles to the generator
  for(i = 0; i < ntParent; i++)
    {
    const Triangle &pt = ptl[i];
    size_t nsplit = 
      (split[i][0] == NOID ? 0 : 1) + 
      (split[i][1] == NOID ? 0 : 1) + 
      (split[i][2] == NOID ? 0 : 1);
    assert(nsplit != 2);

    if(nsplit == 0)
      {
      tmg.AddTriangle(pt.vertices[0], pt.vertices[1], pt.vertices[2], pt.label);
      }
    else if(nsplit == 1)
      {
      j = (split[i][0] != NOID) ? 0 : ( (split[i][1] != NOID) ? 1 : 2 );
      tmg.AddTriangle(pt.vertices[j], pt.vertices[(j+1) % 3], split[i][j], pt.label);
      tmg.AddTriangle(pt.vertices[j], split[i][j], pt.vertices[(j+2) % 3], pt.label);
      }
    else if(nsplit == 3)
      {
      tmg.AddTriangle(pt.vertices[0], split[i][2], split[i][1], pt.label);
      tmg.AddTriangle(split[i][2], pt.vertices[1], split[i][0], pt.label);
      tmg.AddTriangle(split[i][1], split[i][0], pt.vertices[2], pt.label);
      tmg.AddTriangle(split[i][0], split[i][1], split[i][2], pt.label);
      }
    }

  // Populate the new mesh
  tmg.GenerateMesh(); 

  // Now, assign Loop weights to the new vertices
  // Create a mutable sparse matrix representation for the weights
  MutableSparseMatrix W(nvChild, nvParent);

  // Visit each of the even vertices, assigning it an id and weights
  for(i = 0; i < nvParent; i++)
    SetEvenVertexWeights(W, i, parent, i, flat_mode);

  // Visit the odd vertices 
  for(i = 0; i < ntParent; i++) for(j = 0; j < 3; j++)
    if(split[i][j] != NOID && W.empty_row(split[i][j]))
      SetOddVertexWeights(W, split[i][j], parent, i, j, flat_mode);

  // Copy the sparse matrix into immutable form
  child->weights.SetFromVNL(W);

  // If the parent's parent is not NULL, we need to multiply the sparse
  // matrices of the parent and child
  if(parent->parent)
    ImmutableSparseMatrix<double>::Multiply(
      child->weights, child->weights, parent->weights);

  // Compute the walks in the child mesh
  child->ComputeWalks();
}

void SubdivisionSurface::Subdivide(const MeshLevel *parent, MeshLevel *child, bool flat_mode)
{
  // Get the numbers of triangles before and after
  size_t ntParent = parent->triangles.size();
  size_t ntChild = 4 * ntParent;
  size_t i, j;

  // Connect the child to the parent
  child->parent = parent;

  // Initialize the number of vertices in the new mesh
  child->nVertices = 0;

  // Subdivide each triangle into four
  child->triangles.resize(ntChild);
  std::fill(child->triangles.begin(), child->triangles.end(), Triangle());

  for (i = 0; i < ntParent; i++)
  {
    // Get pointers to the four ctren
    const Triangle &pt = parent->triangles[i];
    Triangle *ct[4] = {
                        &child->triangles[4*i], &child->triangles[4*i+1],
                        &child->triangles[4*i+2], &child->triangles[4*i+3]};

    // Copy the label
    for (j = 0; j < 4; j++)
      ct[j]->label = pt.label;

    // Set the neighbors within this triangle
    for (j = 0; j < 3; j++)
    {
      // Assign the neighborhoods within the pt triangle
      ct[j]->neighbors[j] = 4*i + 3;
      ct[3]->neighbors[j] = 4*i + j;
      ct[j]->nedges[j] = j;
      ct[3]->nedges[j] = j;

      // Assign neighborhoods outside the pt triangle
      if (pt.neighbors[(j+1) % 3] != NOID)
      {
        ct[j]->neighbors[(j+1) % 3] =
          pt.neighbors[(j+1) % 3] * 4 + ((pt.nedges[(j+1) % 3] + 1) % 3);
        ct[j]->nedges[(j+1) % 3] = pt.nedges[(j+1) % 3];
      }
      if (pt.neighbors[(j+2) % 3] != NOID)
      {
        ct[j]->neighbors[(j+2) % 3] =
          pt.neighbors[(j+2) % 3] * 4 + ((pt.nedges[(j+2) % 3] + 2) % 3);
        ct[j]->nedges[(j+2) % 3] = pt.nedges[(j+2) % 3];
      }
    }
  }

  // Compute the number of edges in the pt graph
  size_t neParent = 0;
  for(i = 0; i < ntParent; i++) for(j = 0; j < 3; j++)
      neParent += (parent->triangles[i].neighbors[j] == NOID) ? 2 : 1;
  neParent >>= 1;

  // Assign vertex ids and weights. Under this scheme, the vertex ids
  // of the pt and ct should match (for even vertices) and odd
  // vertices appear at the end of the list

  // Create a mutable sparse matrix representation for the weights
  MutableSparseMatrix W(parent->nVertices + neParent, parent->nVertices);

  // Visit each of the even vertices, assigning it an id and weights
  for(i = 0; i < ntParent; i++)
    {
    for(j = 0; j < 3; j++)
      {
      if (child->triangles[4 * i + j].vertices[j] == NOID)
        {
        RecursiveAssignVertexLabel(child, 4*i+j, j, parent->triangles[i].vertices[j]);
        size_t ivc = child->triangles[4*i+j].vertices[j];
        size_t ivp = parent->triangles[i].vertices[j];
        assert(ivc == ivp); 
        SetEvenVertexWeights(W, ivc, parent, ivp, flat_mode);
        }
      }
    }

  // Visit each of the odd vertices, assigning it an id, and weights
  child->nVertices = parent->nVertices;
  for(i = 0; i < ntParent; i++) 
    {
    for(j = 0; j < 3; j++)
      {
      if (child->triangles[4 * i + 3].vertices[j] == NOID)
        {
        RecursiveAssignVertexLabel(child, 4*i+3, j, child->nVertices++);
        size_t ivc = child->triangles[4*i+3].vertices[j];
        SetOddVertexWeights(W, ivc, parent, i, j, flat_mode);
        }
      }
    }

  // Copy the sparse matrix into immutable form
  child->weights.SetFromVNL(W);

  // If the parent's parent is not NULL, we need to multiply the sparse
  // matrices of the parent and child
  if(parent->parent)
    ImmutableSparseMatrix<double>::Multiply(
      child->weights, child->weights, parent->weights);

  // Compute the walks in the child mesh
  child->ComputeWalks();
}

void
SubdivisionSurface::
RecursiveSubdivide(const MeshLevel *src, MeshLevel *dst, size_t n, bool flat_mode)
{
  // N has to be greater than 1
  if(n == 0)
  { 
    // Set the source to have the same structure as the destination
    *dst = *src; 

    // Make the destination be its own root
    dst->SetAsRoot();

    // Done
    return; 
  }
  
  else if(n == 1)
  { Subdivide(src, dst, flat_mode); return; }

  // Create an array of intermedial mesh levels
  MeshLevel *temp = new MeshLevel[n - 1];

  // Subdivide the intermediate levels
  const MeshLevel *parent = src;
  for(size_t i = 0; i < n-1; i++)
  {
    Subdivide(parent, temp + i, flat_mode);
    parent = temp + i;
  }

  // Subdivide the last level
  Subdivide(parent, dst, flat_mode);

  // Delete the intermediates
  delete[] temp;

  // Set the parent pointer in dst to src (bypass intermediates)
  dst->parent = src;
}

void
SubdivisionSurface::
PickTrianglesWithLabel(const MeshLevel &src, size_t label, MeshLevel &dst, std::vector<size_t> &vtx_full_to_sub)
{
  // We need to maintain a mapping from new vertices to parent vertices
  std::vector<size_t> vtx_sub_to_full;
  vtx_full_to_sub.resize(src.nVertices, NOID);
  std::vector<size_t> tri_sub_vtx;

  // Figure out the new triangles and the vertex mapping
  for(size_t t = 0; t < src.triangles.size(); t++)
    {
    if(src.triangles[t].label == label)
      {
      // This triangle gets copied. Remap its vertices
      for(size_t k = 0; k < 3; k++)
        {
        // Figure out the new index of the vertex
        size_t v_full = src.triangles[t].vertices[k];
        if(vtx_full_to_sub[v_full] == NOID)
          {
          vtx_full_to_sub[v_full] = vtx_sub_to_full.size();
          vtx_sub_to_full.push_back(v_full);
          }

        // Store the vertex
        tri_sub_vtx.push_back(vtx_full_to_sub[v_full]);
        }
      }
    }

  // Create a mesh builder for the destination mesh
  TriangleMeshGenerator tmg(&dst, vtx_sub_to_full.size());

  // Add all the new vertices
  for(size_t j = 0; j < tri_sub_vtx.size(); j += 3)
    tmg.AddTriangle(tri_sub_vtx[j], tri_sub_vtx[j+1], tri_sub_vtx[j+2], label);

  // Generate the mesh
  tmg.GenerateMesh();

  // We need a matrix from new vertices to parent vertices. The matrix is just ones and zeros
  ImmutableSparseMatrix<double>::VNLSourceType w_vnl(vtx_sub_to_full.size(), src.nVertices);
  for(size_t j = 0; j < vtx_sub_to_full.size(); j++)
    w_vnl(j, vtx_sub_to_full[j]) = 1.0;

  // Set the parent
  dst.parent = &src;

  // Set the weights
  dst.weights.SetFromVNL(w_vnl);

  // If the parent's parent is not NULL, we need to multiply the sparse
  // matrices of the parent and child
  if(src.parent)
    ImmutableSparseMatrix<double>::Multiply(
      dst.weights, dst.weights, src.weights);
}

std::set<size_t> 
SubdivisionSurface::
GetBoundaryTriangles(const MeshLevel *src)
{
  std::set<size_t> bag;
  for(size_t t = 0; t < src->triangles.size(); t++)
    for(size_t j = 0; j < 3; j++)
      {
      EdgeWalkAroundVertex it(src, src->triangles[t].vertices[j]);
      if(it.IsOpen())
        {
        bag.insert(t); 
        break;
        }
      }
  return bag;
}

void
SubdivisionSurface::
RecursiveSubdivideBoundary(const MeshLevel *src, MeshLevel *dst, size_t n, bool flat_mode)
{
  // N has to be greater than 1
  if(n == 0)
  { 
    // Set the source to have the same structure as the destination
    *dst = *src; 

    // Make the destination be its own root
    dst->SetAsRoot();

    // Done
    return; 
  }
  
  else if(n == 1)
    { 
    SubdivideSelected(src, dst, GetBoundaryTriangles(src), flat_mode); 
    return; 
    }

  // Create an array of intermedial mesh levels
  MeshLevel *temp = new MeshLevel[n - 1];

  // Subdivide the intermediate levels
  const MeshLevel *parent = src;
  for(size_t i = 0; i < n-1; i++)
  {
    SubdivideSelected(parent, temp + i, GetBoundaryTriangles(parent), flat_mode);
    parent = temp + i;
  }

  // Subdivide the last level
  SubdivideSelected(parent, dst, GetBoundaryTriangles(parent), flat_mode);

  // Delete the intermediates
  delete[] temp;

  // Set the parent pointer in dst to src (bypass intermediates)
  dst->parent = src;
}

void SubdivisionSurface::ExportLevelToVTK(const MeshLevel &src, vtkPolyData *mesh)
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

void SubdivisionSurface::ImportLevelFromVTK(vtkPolyData *mesh, MeshLevel &dest)
{
  // Create the generator
  TriangleMeshGenerator tmgen(&dest, mesh->GetNumberOfPoints());

  // Prepare the VTK mesh
  mesh->BuildCells();

  // For each triangle, compute the neighbor. This can be done by first enumerating
  // all the edges in the mesh. For each edge there will be one or two triangles
  for(size_t i = 0; i < (size_t) mesh->GetNumberOfCells(); i++)
  {
    // Get the points from the current triangle
    vtkIdType npts;
    const vtkIdType *pts;
    mesh->GetCellPoints(i, npts, pts);

    // If the number of points is not 3, return with an exception
    if(npts != 3) throw string("Mesh contains cells other than triangles");

    // Add the vertices
    tmgen.AddTriangle(pts[0], pts[1], pts[2]);
  }

  // Generate the mesh after all triangles added
  tmgen.GenerateMesh();

  // Set the mesh's parent to NULL
  dest.parent = NULL;
}

void SubdivisionSurface
::ApplySubdivision(vtkPolyData *src, vtkPolyData *target, MeshLevel &m)
{
  size_t i;

  // Initialize the target mesh
  vtkPoints *tpoints = vtkPoints::New();
  vtkCellArray *tcells = vtkCellArray::New();

  // Add the points to the output mesh
  for(i = 0; i < m.nVertices; i++)
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

void SubdivisionSurface
::ApplySubdivision(SMLVec3d *xSrc, double *rhoSrc, SMLVec3d *xTrg, double *rhoTrg, MeshLevel &m)
{
  size_t i;

  // Add the points to the output mesh
  for(i = 0; i < m.nVertices; i++)
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

void SubdivisionSurface
::ApplySubdivision(const double *xsrc, double *xdst, size_t nComp, MeshLevel &m)
{
  typedef ImmutableSparseMatrix<double>::RowIterator IteratorType;

  // Get the total number of values
  size_t n = nComp * m.nVertices;

  // Set the target array to zero
  std::fill_n(xdst, n, 0.0);

  // Iterate over the output vertices
  for(size_t iv = 0; iv < m.nVertices; iv++)
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

bool SubdivisionSurface::CheckMeshLevel (MeshLevel *mesh)
{
  // Check the following rules for all triangles in the mesh
  // 1. T[T[i].nbr[j]].nbr[T[i].ne[j]] == i for all i, j    (i am my neighbors neighbor)
  // 2. T[T[i].nbr[j]].ne[T[i].ne[j]] == j for all i, j    (i am my neighbors neighbor)
  // 3. T[i].v[j] == T[T[i].nbr[j+k]].v[T[i].ne[j+k]+k] for all i, j, k=(1, 2),
  // with modulo 3 addition
  size_t nerr = 0;
  for(size_t i = 0; i < mesh->triangles.size(); i++)
  {
    Triangle &t = mesh->triangles[i];
    for(size_t j = 0; j < 3; j++)
    {
      if(t.neighbors[j] != NOID)
      {
        Triangle &tn = mesh->triangles[t.neighbors[j]];
        if(tn.neighbors[t.nedges[j]] != i)
        {
          cout << "Error " << nerr++ <<
          " Rule 1 violated for i = " << i << " and j = " << j << endl;
        }
        if(tn.nedges[t.nedges[j]] != (int) j)
        {
          cout << "Error " << nerr++ <<
          " Rule 2 violated for i = " << i << " and j = " << j << endl;
        }
      }
      for(size_t k = 1; k < 3; k++)
      {
        if(t.neighbors[(j+k) % 3] != NOID)
        {
          Triangle &tk = mesh->triangles[t.neighbors[(j+k) % 3]];
          if(t.vertices[j] != tk.vertices[(t.nedges[(j+k) % 3] + k ) % 3])
          {
            cout << "Error " << nerr++ <<
            " Rule 3 violated for i = " << i << ", j = " << j << " and k = " << k << endl;
          }
        }
      }
    }
  }

  return (nerr == 0);
}







