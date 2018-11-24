#ifndef __PointSetOctree_h_
#define __PointSetOctree_h_

#include <vnl_matrix.h>
#include <vnl_vector.h>
#include <vnl_vector_fixed.h>
#include <vector>

template <class TFloat, unsigned int VDim>
class PointSetOctree
{
public:

  // Vectors and matrices
  typedef vnl_vector_fixed<TFloat, VDim> VecN;
  typedef vnl_matrix<TFloat> Matrix;

  // Number of children
  static const unsigned int NC = 1 << VDim;
  static const unsigned int NIL = (unsigned int)(-1);

  // Actual node
  struct Node 
    {
    // Indices of child nodes
    unsigned int children[NC];

    // Whether the node is a leaf
    bool isleaf;

    // Leaf nodes contain vertices
    unsigned int vtx;

    // Center of the node
    VecN x_ctr, x_rad;

    // Number of children in the node
    unsigned int nv;

    // Initialize a node as leaf
    Node(unsigned int vertex_id)
      {
      isleaf = true;
      for(unsigned int i = 0; i < NC; i++)
        children[i] = NIL;
      vtx = vertex_id;
      nv = 1;
      }
    };

  // Data associated with a node - stuff that depends on q/p
  struct NodeData
    {
    // The accumulated data for computing the center and total
    // potential of the node
    VecN q_ctr, p_total;

    // The actual bounds of the octree
    VecN q_lb, q_ub;
    };

  // Place a list of points into an octree
  void Build(const Matrix &q);

  // Update the node data with new positions and weights
  void Update(const Matrix &q, const Matrix &p, unsigned int i_node = 0);

  // Approximate Gaussian transform at a point
  void Approximate(const VecN &x, VecN &p, unsigned int &count, unsigned int &found_node,
    TFloat d2_rel_thresh, TFloat f, unsigned int i_node = 0);

protected:

  // Insert a vertex into the octree
  void InsertVertex(unsigned int i_node, unsigned int depth, const Matrix &q, unsigned int i_vtx);

  // Create a child node inside an existing one
  void CreateChildNode(Node &node, unsigned int k, unsigned int i_vtx);



  std::vector<Node> nodes;
  std::vector<NodeData> node_data;

};



#endif // __PointSetOctree_h_
