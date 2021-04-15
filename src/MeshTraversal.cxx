#include "MeshTraversal.h"
#include "MedialException.h"
#include <algorithm>

using namespace std;

inline bool NeighborPairPredicate (
  const pair<size_t, NeighborInfo> &p1,
  const pair<size_t, NeighborInfo> &p2)
{
  return p1.first < p2.first;
}

void TriangleMesh::ComputeWalks()
{
  // Create a temporary mutable structure that represents walk info
  typedef std::pair<size_t, NeighborInfo> Entry;
  typedef std::list<Entry> Row;
  std::vector<Row> walks(this->nVertices);

  // Loop over all triangles, all vertices
  for(size_t t = 0; t < this->triangles.size(); t++) for(size_t v = 0; v < 3; v++)
    {
      // This is the index of the current vertex
      size_t ivtx = this->triangles[t].vertices[v];

      // Check if the walk around the vertex has already been generated
      if(walks[ivtx].size() == 0)
      {
        // The current position in the walk
        size_t tWalk = t; short vWalk = v;

        // The position at which the walk will loop around
        size_t tLoop = t;

        // Walk until reaching a NOID triangle or looping around
        // cout << "Walk around vtx. " << ivtx << endl;
        do
        {
          // Get a reference to the current triangle
          Triangle &T = this->triangles[tWalk];

          // Represent the next triangle in the walk
          size_t tNext = T.neighbors[ror(vWalk)];
          short vNext = ror(T.nedges[ror(vWalk)]);

          // Put the current edge in the triangle into the list
          // cout << "Fwd walk: visiting vtx. " << T.vertices[rol(vWalk)] <<
          //    "; behind: (" << tWalk << "," << vWalk <<
          //    "); ahead: (" << tNext << "," << vNext << ")" <<  endl;
          walks[ivtx].push_back( make_pair(
                                   T.vertices[rol(vWalk)],
                                   NeighborInfo(tNext, vNext, tWalk, vWalk)));

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
            Triangle &T = this->triangles[tWalk];

            // Represent the next triangle in the walk
            size_t tNext = T.neighbors[rol(vWalk)];
            short vNext = rol(T.nedges[rol(vWalk)]);

            // Put the current edge in the triangle into the list
            // cout << "Rev walk: visiting vtx. " << T.vertices[ror(vWalk)] <<
            //   "; behind: (" << tNext << "," << vNext <<
            //   "); ahead: (" << tWalk << "," << vWalk << ")" << endl;
            walks[ivtx].push_front( make_pair(
                                      T.vertices[ror(vWalk)],
                                      NeighborInfo(tWalk, vWalk, tNext, vNext)));

            // Update the position in the walk
            tWalk = tNext; vWalk = vNext;
          }
          while(tWalk != NOID);
        }
        else
        {
          // Rotate the walk so that the smallest member is the first element
          // (for consistency with MMA code)
          rotate(walks[ivtx].begin(),
                 min_element(walks[ivtx].begin(), walks[ivtx].end(), &NeighborPairPredicate),
                 walks[ivtx].end());
        }
      }
    }

  // Now we have visited all the vertices and computed a walk around each one.
  // All that is left is to transfer this result into a sparse matrix
  this->nbr.SetFromSTL(walks, this->triangles.size());
}

/******************************************************************
 CODE FOR THE LOOP TANGENT SCHEME
 *****************************************************************/
void LoopTangentScheme::SetMesh(TriangleMesh *mesh)
{
  const TriangleMesh::NeighborMatrix &nmat = mesh->GetNeighborMatrix();

  // Create new sparse array structure based on the mesh, but expanded to
  // include the diagonal of the matrix
  size_t i;
  size_t nSparse = nmat.GetNumberOfSparseValues() + nmat.GetNumberOfRows();
  size_t *xRowIndex = new size_t[nmat.GetNumberOfRows() + 1];
  size_t *xColIndex = new size_t[nSparse];
  Weight *xSparse = new Weight[nSparse];

  // Populate the rows and colums
  xRowIndex[0] = 0;
  for(i = 0; i < nmat.GetNumberOfRows(); i++)
    {
    // Get the number of neighbors
    size_t nnbr = nmat.GetRowIndex()[i+1] - nmat.GetRowIndex()[i];

    // One more column for every row
    xRowIndex[i+1] = xRowIndex[i] + nnbr + 1;

    // Fill the first value in the column: diagonal element
    xColIndex[xRowIndex[i]] = i;

    // Fill the rest of the column values
    for(size_t j = 0; j < nnbr; j++)
      xColIndex[xRowIndex[i] + j + 1] = nmat.GetColIndex()[nmat.GetRowIndex()[i] + j];
    }

  // Clear out the sparse values
  for(i = 0; i < nSparse; i++)
    xSparse[i].w[0] = xSparse[i].w[1] = xSparse[i].w[2] = 0.0;

  // Now, compute the weights for the loop scheme
  for(i = 0; i < nmat.GetNumberOfRows(); i++)
    {
    // Create an iterator around the vertex
    EdgeWalkAroundVertex it(mesh, i);

    // Get the valence of the vertex
    size_t n = it.Valence();

    // If the vertex is internal, weights are easy to compute
    if(!it.IsOpen())
      {
      // Terms for the limit surface calculation
      double beta = (n > 3) ? 3.0 / (8.0 * n) : 3.0 / 16.0;
      double xi = 1.0 / (3. / (8. * beta) + n);

      for(size_t j = 0; !it.IsAtEnd(); ++it, ++j)
        {
        size_t k = it.GetPositionInMeshSparseArray() + i + 1;

        // Compute the weights for the tangents
        xSparse[k].w[0] = cos(j * vnl_math::pi * 2.0 / n);
        xSparse[k].w[1] = sin(j * vnl_math::pi * 2.0 / n);

        // Compute the weight for the limit surface
        xSparse[k].w[2] = xi;
        }

      size_t kSelf = xRowIndex[i];
      xSparse[kSelf].w[2] = 1.0 - n * xi;
      }
    else
      {
      // Get the index of the first neighbor
      size_t kFirst = it.GetPositionInMeshSparseArray() + i + 1;

      // Move the index to the last neighbor
      it.GoToLastEdge();
      size_t kLast = it.GetPositionInMeshSparseArray() + i + 1;

      // Get the index of the vertex itself
      size_t kSelf = xRowIndex[i];

      // We can now set the along-the-edge tangent vector
      xSparse[kFirst].w[1] = 1.0;
      xSparse[kLast].w[1] = -1.0;

      // Now we branch according to the valence
      if(n == 2)
        {
        xSparse[kFirst].w[0] = -1.0;
        xSparse[kLast].w[0] = -1.0;
        xSparse[kSelf].w[0] = 2.0;
        }
      else if (n == 3)
        {
        // Move iterator to the previous (middle) vertex
        it.GoToFirstEdge(); ++it;

        // Get the corresponding A index
        size_t kMiddle = it.GetPositionInMeshSparseArray() + i + 1;

        // Set the weights
        xSparse[kMiddle].w[0] = -1.0;
        xSparse[kSelf].w[0] = 1.0;
        }
      else
        {
        // Compute the angle theta
        double theta = vnl_math::pi / (n - 1);
        xSparse[kFirst].w[0] = xSparse[kLast].w[0] = sin(theta);

        // Assign the weights to intermediate vertices
        it.GoToFirstEdge(); ++it;
        for(size_t j = 1; j < n - 1; ++j, ++it)
          {
          size_t kInt = it.GetPositionInMeshSparseArray() + i + 1;
          xSparse[kInt].w[0] = (2.0 * cos(theta) - 2.0) * sin(theta * j);
          }
        }

      // Limit surface is quite simple
      xSparse[kSelf].w[2] = 0.6;
      xSparse[kFirst].w[2] = 0.2; 
      xSparse[kLast].w[2] = 0.2; 
      }
    }

  // Create the sparse array
  W.SetArrays(
    nmat.GetNumberOfRows(), nmat.GetNumberOfColumns(),
    xRowIndex, xColIndex, xSparse);
}

/******************************************************************
 CODE FOR THE MESH GENERATOR
 *****************************************************************/
TriangleMeshGenerator::TriangleMeshGenerator(
  TriangleMesh *target, size_t nVertices)
{
  xTargetMesh = target;  

  xTargetMesh->nVertices = nVertices;
  xTargetMesh->triangles.clear();
}

void TriangleMeshGenerator::AddTriangle(size_t v0, size_t v1, size_t v2, size_t label)
{
  size_t v[3] = {v0, v1, v2};
  size_t i = xTargetMesh->triangles.size();

  // Associate each half-edge with a triangle
  Triangle tri;
  tri.label = label;
  for(size_t j = 0; j < 3; j++)
  {
    // Set the vertices in each triangle
    tri.vertices[j] = v[j];

    // Create the key and value
    HalfEdge he(v[(j+1) % 3], v[(j+2) % 3]);
    TriangleRep trep(i, (short) j);

    // Insert the half-edge and check for uniqueness
    pair<TriangleMap::iterator, bool> rc = tmap.insert(make_pair(he, trep));
    if(!rc.second)
      {
      std::ostringstream oss;
      oss << "Half-edge [" << v[(j+1) % 3] << "," << v[(j+2) % 3] << "] " <<
        "appears twice in the mesh, making the mesh irregular. This is often" <<
        "caused by meshes where an entire triangle lies on the medial edge" <<
        "curve. The solution is to remove the bad triangle from the mesh" << std::endl;
      throw MedialModelException(oss.str().c_str());
      }
  }

  // Store the triangle
  xTargetMesh->triangles.push_back(tri);
}

void TriangleMeshGenerator::GenerateMesh()
{
  // Take a pass through all the half-edges. For each, set the 
  // corresponding triangle's neighbor and neighbor index
  for(TriangleMap::iterator it = tmap.begin(); it != tmap.end(); ++it)
  {
    TriangleRep &trep = it->second;
    HalfEdge opposite(it->first.second, it->first.first);
    TriangleMap::iterator itopp = tmap.find(opposite);
    if(itopp != tmap.end())
    {
      xTargetMesh->triangles[trep.first].neighbors[trep.second] = itopp->second.first;
      xTargetMesh->triangles[trep.first].nedges[trep.second] = itopp->second.second;
    }
  }

  // Finally, compute the walks in this mesh
  xTargetMesh->ComputeWalks();  
}


/******************************************************************
 CODE FOR THE MESH GRADIENT COMPUTER
 *****************************************************************/
void MeshGradientComputer::SetMesh(TriangleMesh *mesh)
{
  // Store the mesh
  this->xMesh = mesh;
  xLoop.SetMesh(mesh);

  // Initialize the sparse matrices
  jacX.SetFromReference(xLoop.GetWeightMatrix(), Mat(0.0));
  jacF.SetFromReference(xLoop.GetWeightMatrix(), Vec(0.0));
  
  // Initialize the gradient vector
  gradF.resize(mesh->nVertices);
}

void MeshGradientComputer::ComputeGradient(Vec *x, double *f, bool flagJacobian)
{
  this->ComputeGradient(x, sizeof(Vec), f, sizeof(double), flagJacobian);
}

void MeshGradientComputer::ComputeGradient(
    void *xVecPtr, size_t xVecStride, 
    void *xFunPtr, size_t xFunStride,
    bool flagJacobian)
{
  // Cast the arrays to char arrays so we can add strides
  char *caVecPtr = static_cast<char *>(xVecPtr);
  char *caFunPtr = static_cast<char *>(xFunPtr);

  // Get the row and column index of the sparse matrix
  size_t *xRowIdx = jacX.GetRowIndex();
  size_t *xColIdx = jacX.GetColIndex();

  // Get the weights as an array
  const LoopTangentScheme::Weight *xWeights 
    = xLoop.GetWeightMatrix().GetSparseData();

  // Loop over all vertices
  for(size_t i = 0; i < xMesh->nVertices; i++)
    {
    double Fi[2];
    Vec Xi[2];

    // Compute the partial derivatives
    for(size_t d = 0; d < 2; d++)
      {
      // Initialize the accumulators
      Fi[d] = 0; Xi[d].fill(0.0);

      // Loop over the adjacents
      for(size_t ic = xRowIdx[i]; ic < xRowIdx[i+1]; ic++)
        {
        size_t j = xColIdx[ic];
        Vec x = *reinterpret_cast<Vec *>(caVecPtr + xVecStride * j);
        double f = *reinterpret_cast<double *>(caFunPtr + xFunStride * j);
        
        double w = xWeights[ic].w[d];
        Fi[d] += w * f;
        Xi[d] += w * x;
        }
      }

    // Compute the covariant metric tensor
    double gcov[2][2], gcon[2][2];
    gcov[0][0] = dot_product(Xi[0], Xi[0]);    
    gcov[0][1] = gcov[1][0] = dot_product(Xi[0], Xi[1]);
    gcov[1][1] = dot_product(Xi[1], Xi[1]);

    // Compute one over the determinant
    double ginvdet = 1.0 / (gcov[0][0] * gcov[1][1] - gcov[1][0] * gcov[1][0]);

    // Compute the contravariant tensor
    gcon[0][0] = ginvdet * gcov[1][1];
    gcon[0][1] = gcon[1][0] = -ginvdet * gcov[0][1];
    gcon[1][1] = ginvdet * gcov[0][0];

    // Compute the gradient
    gradF[i] = 
      (gcon[0][0] * Fi[0] + gcon[1][0] * Fi[1]) * Xi[0] + 
      (gcon[0][1] * Fi[0] + gcon[1][1] * Fi[1]) * Xi[1];

    // Jacobian computation has to be requested
    if(flagJacobian)
      {
      // Now, compute the Jacobian terms. This will require a loop
      // We have:
      // \frac{d gradF}{d X_{,q}} = 
      //   = -F_{,b} (g^{ac}g^{qb}+g^{aq}g^{cb}) X_{,c} X_{,a}^T + F_{,b} g^{qb} I_{3 \cross 3}
      // \frac{d gradF}{d F_{,q}} = 
      //   = -g^{bq} X_{,b}
      // There is a Mathematica notebook that verifies this (TestMeshGradientComputation.nb)
      Mat coeffXi[2]; 
      Vec coeffFi[2];
      size_t a, b, c, q;

      // First compute the simple terms
      for(q=0;q<2;q++) 
        {
        coeffXi[q].fill(0.0); coeffFi[q].fill(0.0);
        for(b=0;b<2;b++) 
          {
          // The F derivative is simple
          coeffFi[q] += gcon[b][q] * Xi[b];

          // So is the second part of the X derivative (add to diagonal)
          double rhs = gcon[b][q] * Fi[b];
          coeffXi[q](0,0) += rhs; coeffXi[q](1,1) += rhs; coeffXi[q](2,2) += rhs;
          }
        }

      // To compute the second part, we loop over a/c first
      for(a=0;a<2;a++) for(c=0;c<2;c++)
        {
        // Compute the matrix once
        Mat M = outer_product(Xi[a], Xi[c]);

        for(q=0;q<2;q++)
          {
          // We will accumulate the weight on the left hand side of X_{,c} X_{,a}^T
          double w = 0;
          for(b=0;b<2;b++)
            w += Fi[b] * (gcon[a][c] * gcon[q][b] + gcon[a][q] * gcon[c][b]);
          coeffXi[q] -= M * w;
          }
        }

      // Now that we have the derivative of the gradient w.r.t. Xu, Xv, Fu, Fv, we
      // can compute the full Jacobian
      for(size_t ic = xRowIdx[i]; ic < xRowIdx[i+1]; ic++)
        {
        // Get the weight 
        double w0 = xWeights[ic].w[0];
        double w1 = xWeights[ic].w[1];
        
        // Set the Jacobians
        jacX.GetSparseData()[ic] = w0 * coeffXi[0] + w1 * coeffXi[1];
        jacF.GetSparseData()[ic] = w0 * coeffFi[0] + w1 * coeffFi[1];
        }
      }
    }  
}

/******************************************************************
 DELAUNAY REMESHING CODE
 *****************************************************************/
#include <deque>


struct FullEdge
{
  size_t t1, t2;
  short q1, q2;
  double length;
  bool mark;
  size_t id;

  void Replace(size_t told, size_t tnew, short qnew)
    {
    if(t1 == told)
      { t1 = tnew; q1 = qnew; }
    else if(t2 == told)
      { t2 = tnew; q2 = qnew; }
    else assert(0);
    }

};

void SanityCheck(
                 TriangleMesh *mesh,
                 std::vector< vnl_vector_fixed<unsigned int, 3> > &tdidx,
                 std::vector<FullEdge> &edges)
  {
  for(size_t i = 0; i < edges.size(); i++)
    {
    FullEdge &edge = edges[i];
    if(edge.t1 != NOID)
      {
      if(tdidx[edge.t1][edge.q1] != edge.id)
        cout << "Error IDX" << endl;
      if(mesh->triangles[edge.t1].neighbors[edge.q1] != edge.t2)
        cout << "Error Mesh" << endl;
      }
    if(edge.t2 != NOID)
      {
      if(tdidx[edge.t2][edge.q2] != edge.id)
        cout << "Error IDX" << endl;
      if(mesh->triangles[edge.t2].neighbors[edge.q2] != edge.t1)
        cout << "Error Mesh" << endl;
      }
    }
  }



void TriangleMesh::MakeDelaunay(vnl_vector_fixed<double,3> *V)
{
  // Structures to store edges
  std::vector<FullEdge> edges;
  typedef vnl_vector_fixed<unsigned int, 3> ShVec;
  std::vector<ShVec> teidx(triangles.size(), ShVec((unsigned int) 0));

  // Create a queue to store edges to process
  std::deque<size_t> edge_que;

  // Generate a list of all edges in the mesh
  for(size_t i = 0; i < triangles.size(); i++)
    {
    for(short q = 0; q < 3; q++)
      {
      // Add edges so that t1 < t2
      if(triangles[i].neighbors[q] == NOID || triangles[i].neighbors[q] > i)
        {
        // Create an edge descriptor
        FullEdge fedge;
        fedge.t1 = i;
        fedge.t2 = triangles[i].neighbors[q];
        fedge.q1 = q;
        fedge.q2 = triangles[i].nedges[q];
        fedge.length = (
          V[this->triangles[i].vertices[ror(q)]] - 
          V[this->triangles[i].vertices[rol(q)]]).magnitude();
        fedge.mark = 1;
        fedge.id = edges.size();
        edges.push_back(fedge);
        edge_que.push_back(fedge.id);

        // Point to this edge
        teidx[fedge.t1][fedge.q1] = fedge.id;
        if(fedge.t2 != NOID)
          teidx[fedge.t2][fedge.q2] = fedge.id;
        }
      }
    }

  int n_flips = 0;


  // Proceed while there are edges on the queue
  while(edge_que.size())
    {
    SanityCheck(this, teidx, edges);

    // Get the current edge of interest
    size_t eij = edge_que.front();
    edge_que.pop_front();

    // Get the triangle and the vertex
    size_t t = edges[eij].t1; short q = edges[eij].q1;

    // If the edge is marked zero, continue (processed opp edge)
    if(edges[eij].mark == 0) 
      {
      continue;
      }

    // Unmark the edge
    edges[eij].mark = 0;

    // Get the opposite triangle    
    size_t topp = edges[eij].t2;
    short qopp = edges[eij].q2;

    // Branch based on boundary or not? Also cannot flip when triangles have different labels
    if(topp != NOID && this->triangles[t].label == this->triangles[topp].label)
      {
      // Get the four edges
      size_t eik = teidx[t][rol(q)];
      size_t ejk = teidx[t][ror(q)];
      size_t eil = teidx[topp][ror(qopp)];
      size_t ejl = teidx[topp][rol(qopp)];

      // Get the four lengths 
      double lij = edges[eij].length;
      double lik = edges[eik].length;
      double ljk = edges[ejk].length;
      double lil = edges[eil].length;
      double ljl = edges[ejl].length;

      // Compute the two cotan weights
      double tan_half_k = sqrt( 
        ((lij - lik + ljk) * (lij + lik - ljk)) / 
        ((lij + lik + ljk) * (-lij + lik + ljk)));
      double tan_half_l = sqrt( 
        ((lij - lil + ljl) * (lij + lil - ljl)) / 
        ((lij + lil + ljl) * (-lij + lil + ljl)));

      // Compute the weight
      double cotan_weight = 
        (1 - tan_half_k * tan_half_k) / tan_half_k + 
        (1 - tan_half_l * tan_half_l) / tan_half_l;

      // If the cotan weight is less than zero, we have to flip
      if(cotan_weight < 0)
        {
        // cout << "Flip " << t << " edge " << q << endl; 
        n_flips++;
        if(n_flips > 0 && n_flips % 1000 == 0)
          {
          cout << "Flipped " << n_flips << "\t Que Size: " << edge_que.size() << endl;
          }

        // Create a pair of triangles and perform all the flip operations
        Triangle tnew, toppnew;
        Triangle told = this->triangles[t], toppold = this->triangles[topp];
        tnew.vertices[q] = told.vertices[ror(q)];
        tnew.neighbors[q] = topp;
        tnew.nedges[q] = qopp;
        tnew.label = told.label;
        
        tnew.vertices[rol(q)] = told.vertices[q];
        tnew.neighbors[rol(q)] = toppold.neighbors[ror(qopp)];
        tnew.nedges[rol(q)] = toppold.nedges[ror(qopp)];
        if(tnew.neighbors[rol(q)] != NOID)
          {
          this->triangles[tnew.neighbors[rol(q)]].neighbors[tnew.nedges[rol(q)]] = t;
          this->triangles[tnew.neighbors[rol(q)]].nedges[tnew.nedges[rol(q)]] = rol(q);
          }

        tnew.vertices[ror(q)] = toppold.vertices[qopp];
        tnew.neighbors[ror(q)] = told.neighbors[rol(q)];
        tnew.nedges[ror(q)] = told.nedges[rol(q)];
        if(tnew.neighbors[ror(q)] != NOID)
          {
          this->triangles[tnew.neighbors[ror(q)]].neighbors[tnew.nedges[ror(q)]] = t;
          this->triangles[tnew.neighbors[ror(q)]].nedges[tnew.nedges[ror(q)]] = ror(q);
          }

        toppnew.vertices[qopp] = toppold.vertices[ror(qopp)];
        toppnew.neighbors[qopp] = t;
        toppnew.nedges[qopp] = q;
        toppnew.label = told.label;
        
        toppnew.vertices[rol(qopp)] = toppold.vertices[qopp];
        toppnew.neighbors[rol(qopp)] = told.neighbors[ror(q)];
        toppnew.nedges[rol(qopp)] = told.nedges[ror(q)];
        if(toppnew.neighbors[rol(qopp)] != NOID)
          {
          this->triangles[toppnew.neighbors[rol(qopp)]].neighbors[toppnew.nedges[rol(qopp)]] = topp;
          this->triangles[toppnew.neighbors[rol(qopp)]].nedges[toppnew.nedges[rol(qopp)]] = rol(qopp);
          }

        toppnew.vertices[ror(qopp)] = told.vertices[q];
        toppnew.neighbors[ror(qopp)] = toppold.neighbors[rol(qopp)];
        toppnew.nedges[ror(qopp)] = toppold.nedges[rol(qopp)];
        if(toppnew.neighbors[ror(qopp)] != NOID)
          {
          this->triangles[toppnew.neighbors[ror(qopp)]].neighbors[toppnew.nedges[ror(qopp)]] = topp;
          this->triangles[toppnew.neighbors[ror(qopp)]].nedges[toppnew.nedges[ror(qopp)]] = ror(qopp);
          }

        // Stick the new triangles in there
        this->triangles[t] = tnew;
        this->triangles[topp] = toppnew;

        // Now, we've got to update the five edges 
        edges[eij].length = (
           V[tnew.vertices[ror(q)]] - 
           V[tnew.vertices[rol(q)]]).magnitude();

        // Change mapping of edges to triangles
        edges[eik].Replace(t, t, ror(q));
        edges[eil].Replace(topp, t, rol(q));
        edges[ejk].Replace(t, topp, rol(qopp));
        edges[ejl].Replace(topp, topp, ror(qopp));

        // Change the triangle-to-edge index
        teidx[t][rol(q)] = eil;
        teidx[t][ror(q)] = eik;
        teidx[topp][rol(qopp)] = ejk;
        teidx[topp][ror(qopp)] = ejl;
        
        // Mark each of the neighbor edges
        size_t enbr[] = {eik, eil, ejk, ejl};
        for(size_t p = 0; p < 4; p++)
          if(edges[enbr[p]].mark == 0)
            {
            edges[enbr[p]].mark = 1;
            edge_que.push_back(enbr[p]);
            }
        }
      }
    }
  cout << endl;
}


// We need to instantiate sparse matrix of NeighborInfo objects
#include "SparseMatrix.txx"
template class ImmutableSparseArray<NeighborInfo>;
template class ImmutableSparseArray<LoopTangentScheme::Weight>;
template class ImmutableSparseArray<MeshGradientComputer::Vec>;
template class ImmutableSparseArray<MeshGradientComputer::Mat>;

