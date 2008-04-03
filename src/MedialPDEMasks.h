#ifndef __MedialPDEMasks_h_
#define __MedialPDEMasks_h_

#include "smlmath.h"
#include <vector>

/******************************************************************************
 * MASK: This class represents a mask on a finite difference grid
 *****************************************************************************/
class FiniteDifferenceMask
{
public:
  // Matrix and vector typedefs
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;
  
  // Constructor: initializes the data to the number of nodes
  FiniteDifferenceMask(size_t nNodes);
  virtual ~FiniteDifferenceMask();

  // Initialize the location of the mask, also specifying rotation
  // Weights must be recomputed after this operation, and flips and
  // transposition should be done afterwards.
  virtual void SetLocation(size_t iu, size_t iv)
    { this->iu = iu; this->iv = iv; }

  // Flip the mask around a cardinal axis by changing the u and/or v offsets
  // of the nodes in the mask. Pass in +1 or -1 for each parameter
  // Weights must be recomputed after this operation
  virtual void FlipMask(bool uFlip, bool vFlip);

  // Transpose the mask, i.e, replace u offsets by v offsets and vice versa.
  // Weights must be recomputed after this operation
  virtual void TransposeMask();

  // Sort the entries in the mask in the row-major order. 
  virtual void SortNodes();

  // Compute the weight matrix associated with the mask
  virtual void ComputeWeights(const double *uGrid, const double *vGrid);

  // Call right after computing weights to allow faster jet computation
  void OptimizeWeights(size_t uNodes, size_t vNodes);
  
  // Compute the first derivatives of a function
  double ComputeOneJet(const double *F, double &Fu, double &Fv);

  // Compute the second derivatives of the function
  double ComputeTwoJet(const Mat &F,
    double &Fu, double &Fv, double &Fuu, double &Fuv, double &Fvv);

  // Get the contribution of a grid node to the first derivatives
  virtual double GetNodeComponentInOneJet( size_t iNode, double &Wu, double &Wv);

  // Get the contribution of a grid node to the second derivatives
  virtual double GetNodeComponentInTwoJet( size_t iNode, 
    double &Wu, double &Wv, double &Wuu, double &Wuv, double &Wvv);

  // Get the location of the site
  size_t GetLocationU() { return iu; }
  size_t GetLocationV() { return iv; }

  // Get the offset of a site covered by the mask
  int GetOffsetU(size_t iNode) { return qu[iNode]; }
  int GetOffsetV(size_t iNode) { return qv[iNode]; }

  // Get the number of nodes in this mask
  size_t Size() { return n; }
  
  // Print a report describing the site
  void PrintReport();

protected:
  // The meaning of the columns in W, i.e., ordering of the weights
  enum WeightMeaning
    { F00 = 0, F10, F01, F20, F02, F11, F30, F03, F21, F12 };

  // Optimized array for weight storage
  static const size_t NUM_WEIGHTS;
  size_t *fwIndex[6];
  double *fwWeight[6];
  size_t fwCount[6];
  
  // A list of weight vectors for the first six partial derivatives.
  Mat W;
  
  // The number of nodes involved
  size_t n;

  // The position of the center node in the grid
  size_t iu, iv, iraw;

  // The u and v offsets of the nodes. Conventionally they are encoded in
  // counterclockwise rings: central point is assigned 0, then starting from
  // the right of the central point we have 1, 2, ..., and so on.
  std::vector<int> qu, qv;

  // Whether the mask has been transposed
  bool flagTransposed;

};

/**
 * This is a parent class for masks that consist of a row and a column of
 * nodes that meet at the central point. These are separable into a finite
 * difference computation in U and another in V.
 */
class CrossFDMask : public FiniteDifferenceMask
{
public:
  CrossFDMask(size_t n) : FiniteDifferenceMask(n) {}

protected:
  /** 
   * This method uses the fact that the samples fall onto two coordinate axis
   * to introduce some heuristics into the weight computation. Basically, the
   * sample points are separated into two arrays, (vertical and horizontal)
   * and for each array, we solve a linear system to compute the weights
   */
  void ComputeWeights(const double *uGrid, const double *vGrid);
};

/**
 * This is a uniformly spaced finite difference mask. It assumes that the
 * grid is evenly spaced around the center node. It uses a 3x3 neighborhood
 * to compute the finite differences 
 */
class UniformFDMask : public FiniteDifferenceMask
{
public:
  // Constructor
  UniformFDMask();
  
  // Compute the weights using standard formula
  void ComputeWeights(const double *uGrid, const double *vGrid);

  // Disable transposing of this mask - it's already symmetric
  void TransposeMask() {}

  // Disable flipping of this mask
  void FlipMask(bool, bool) {}

  // The mask is already sorted
  void SortNodes() {}

private:
  // Grid spacing (it's uniform)
  double du, dv;
};

/**
 * Mask of the following layout. It is used when the grid is non-uniform. It
 * may be used adjacent to corner nodes.
 * 
 *     9  
 *  3  2
 *  4  0  1  8
 *  5  6  7
 *
 */  
class DoublyNonuniformFDMask : public FiniteDifferenceMask
{
public:
  // Constructor: 10 nodes
  DoublyNonuniformFDMask();
};

/** 
 * Mask of the following layout. It's used when the grid is isotropic in one
 * direction but not in another. By default, the U axis is assumed to be
 * anisotropic, but you can transpose the grid easily
 * 
 *
 * 4  3  2
 * 5  0  1  9
 * 6  7  8
 *
 * 
 *
class SinglyNonuniformFDMask : public FiniteDifferenceMask
{
public:
  // Constructor: 10 nodes
  SinglyNonuniformFDMask() : FiniteDifferenceMask(10) {}
  
  // Set the location for the mask. Flips are very important. With 1,1 the
  // mask is located in the upper left corner
  void SetLocation(size_t u, size_t v);
  
  // Compute the weights
  virtual void ComputeWeights(const double *uGrid, const double *vGrid);
};
*/

/**
 * Mask of the following layout. Used at the corners of grids for computing
 * forward/backward difference gradients. The approximation of the gradient
 * will be of second order
 * 
 *  4
 *  2    
 *  0  1  3 
 */  
class CornerFDMask : public CrossFDMask
{
public:
  // Constructor: 5 nodes
  CornerFDMask();
};


/**
 * Mask of the following layout. Used at the borders of grids for computing
 * forward/backward difference gradients. The approximation of the gradient
 * will be of second order
 * 
 *  2
 *  0  1  4
 *  3     
 */  
class BorderFDMask : public CrossFDMask                    
{
public:
  // Constructor: 5 nodes
  BorderFDMask();
};

/**
 * 2
 * 0  1
 */
class SimpleCornerFDMask : public CrossFDMask
{
public:
  SimpleCornerFDMask();
};

/**
 * 2
 * 0  1
 * 3
 */
class SimpleBorderFDMask : public CrossFDMask
{
public:
  SimpleBorderFDMask();
};

#endif
