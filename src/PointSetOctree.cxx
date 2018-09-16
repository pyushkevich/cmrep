#include "PointSetOctree.h"

template <class TFloat, unsigned int VDim>
void
PointSetOctree<TFloat, VDim>
::CreateChildNode(Node &node, unsigned int k, unsigned int i_vtx)
{
  // Assign the child node
  node.children[k] = nodes.size();

  // Create a new node
  nodes.push_back(Node(i_vtx));
  Node &child = nodes.back();

  // TODO: this is redundant - we can just have a global list of radii
  child.x_rad = node.x_rad * (TFloat) 0.5;

  // Compute the center of child node
  for(unsigned int a = 0; a < VDim; a++)
    {
    if(k & (1 << a))
      child.x_ctr[a] = node.x_ctr[a] + node.x_rad[a] * 0.5;
    else
      child.x_ctr[a] = node.x_ctr[a] - node.x_rad[a] * 0.5;
    }
}

template <class TFloat, unsigned int VDim>
void
PointSetOctree<TFloat, VDim>
::InsertVertex(unsigned int i_node, unsigned int depth, const Matrix &q, unsigned int i_vtx)
{
  // Get the current node
  Node &node = nodes[i_node];

  // If the node is a leaf, it will require splitting
  if(node.isleaf)
    {
    // Find what octant the existing vertex should go to
    unsigned int k = 0;
    for(unsigned int a = 0; a < VDim; a++)
      if(q(node.vtx, a) > node.x_ctr[a])
        k += (1 << a);

    // Create the child node for that octant
    CreateChildNode(node, k, node.vtx);

    // Update the current node's status
    node.isleaf = false;
    node.vtx = NIL;

    // Repeat the same operation -- now the node is no longer leaf
    InsertVertex(i_node, depth, q, i_vtx);
    }

  else
    {
    // The node is not a leaf. Find which child corresponds to the node being inserted
    unsigned int k = 0;
    for(unsigned int a = 0; a < VDim; a++)
      if(q(i_vtx, a) > node.x_ctr[a])
        k += (1 << a);

    // Does a child for this location already exist?
    if(node.children[k] == NIL)
      {
      // If not, create a new leaf node
      CreateChildNode(node, k, i_vtx);
      }
    else
      {
      // Does a child already exist there? Then perform recursion
      InsertVertex(node.children[k], depth+1, q, i_vtx);
      }

    // Increment number of children
    node.nv++;
    }
}

/**
 * Build an octree from scratch
 */
template <class TFloat, unsigned int VDim>
void
PointSetOctree<TFloat, VDim>
::Build(const Matrix &q)
{
  // Initialize the data storage
  unsigned int n = q.rows(); unsigned int n_reserve = (unsigned int)(n * log(n) / log(2.0));
  nodes.clear(); nodes.reserve(n_reserve);

  // Compute the bounding box of the octree
  VecN lb(0.0), ub(0.0);
  for(unsigned int i = 0; i < q.rows(); i++)
    {
    for(unsigned int a = 0; a < VDim; a++)
      {
      lb[a] = (i == 0) ? q(i,a) : std::min(q(i,a), lb[a]);
      ub[a] = (i == 0) ? q(i,a) : std::max(q(i,a), ub[a]);
      }
    }

  // Create the root node and place the first vertex in it
  nodes.push_back(Node(0));
  Node &root = nodes.front();

  // Set the bounding box of the center
  root.x_ctr = (ub + lb) * (TFloat) 0.5;
  root.x_rad = (ub - lb) * (TFloat) 0.5;

  // Insert into the octree
  for(unsigned int i = 1; i < q.rows(); i++)
    InsertVertex(0, 0, q, i);

  // Allocate the node data
  node_data.resize(nodes.size(), NodeData());
}

template <class TFloat, unsigned int VDim>
void
PointSetOctree<TFloat, VDim>
::Update(const Matrix &q, const Matrix &p, unsigned int i_node)
{
  const Node &node = nodes[i_node];
  NodeData &nd = node_data[i_node];

  if(node.isleaf)
    {
    // If the node is a leaf, it's center of mass and total weight/momentum is 
    // just that of its vertex, and it's radius is zero
    nd.q_ctr = q.get_row(node.vtx);
    nd.q_lb = nd.q_ub = nd.q_ctr;
    nd.p_total = p.get_row(node.vtx);
    }
  else
    {
    // Go over all the children of the node and incorporate their information
    nd.q_ctr.fill(0.0);
    nd.p_total.fill(0.0);

    bool found_child = false;
    for(unsigned int j = 0; j < NC; j++)
      {
      if(node.children[j] != NIL)
        {
        // Descend into the child node
        Update(q, p, node.children[j]);

        const Node &child = nodes[node.children[j]];
        const NodeData &cnd = node_data[node.children[j]];

        // Update the center and total momentum
        nd.q_ctr += cnd.q_ctr * (TFloat) child.nv;
        nd.p_total += cnd.p_total;

        // Update the bounding boxes
        if(!found_child)
          {
          nd.q_lb = cnd.q_lb;
          nd.q_ub = cnd.q_ub;
          found_child = true;
          }
        else
          {
          for(unsigned int a = 0; a < VDim; a++)
            {
            nd.q_lb[a] = std::min(nd.q_lb[a], cnd.q_lb[a]);
            nd.q_ub[a] = std::max(nd.q_ub[a], cnd.q_ub[a]);
            }
          }
        }
      }

    // Normalize the center
    nd.q_ctr *= (1.0 / node.nv);
    }
}

template <class TFloat, unsigned int VDim>
void
PointSetOctree<TFloat, VDim>
::Approximate(const VecN &x, VecN &px, unsigned int &count, 
  TFloat d2_rel_thresh, TFloat f, unsigned int i_node)
{
  // Get the node info
  const Node &node = nodes[i_node];
  NodeData &nd = node_data[i_node];

  // Compute the delta from the point to the node center (we need it anyway)
  VecN dx = x - nd.q_ctr;

  // Do we descend down this node or use it's center and bulk potential?
  bool descend;
  if(node.isleaf)
    {
    descend = false;
    }
  else if(i_node == 0)
    {
    descend = true;
    }
  else
    {
    // Use the relative distance threshold test
    TFloat d2_rel = 0.0;
    for(unsigned int a = 0; a < VDim; a++)
      {
      TFloat del = dx[a] / (nd.q_ub[a] - nd.q_lb[a]);
      d2_rel += del * del;
      }
    descend = d2_rel < d2_rel_thresh;
    }

  if(descend)
    {
    for(unsigned int j = 0; j < NC; j++)
      if(node.children[j] != NIL)
        Approximate(x, px, count, d2_rel_thresh, f, node.children[j]);
    }
  else
    {
    px += nd.p_total * exp(f * dx.squared_magnitude());
    count++;
    }
}

template class PointSetOctree<float, 2>;
template class PointSetOctree<float, 3>;
template class PointSetOctree<double, 2>;
template class PointSetOctree<double, 3>;
