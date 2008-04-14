#ifndef __CartesianGridMedialIterationContext_h_
#define __CartesianGridMedialIterationContext_h_

#include "MedialIterationContext.h"
#include "MeshTraversal.h"

class CartesianGridMedialIterationContext
: public MedialIterationContext
{
public:

  /** 
   * Constructor, takes the size of the cartesian grid. The problem with the
   * cartesian grid is that its faces are quads and not triangles. So in a
   * way, picking one direction in which to compute triangles will somehow
   * bias the computations. How do we get around that?
   * 
   * For now we define it like this:
   * *---*---*
   * |1/2|3/4|
   * *---*---*
   * *5/6|7/8|
   * *---*---*
   *
   * Another question in the order in which we iterate: first along the u axis
   * or first along the v axis? Here we will iterate along the u axis first,
   * i.e., ------------->
   *       ------------->
   *       ------------->
   * where u in the horizontal axis.      
   */
  CartesianGridMedialIterationContext(size_t nu, size_t nv)
    {
    // Save nu and nv
    this->nu = nu; this->nv = nv;

    // Set the sizes of the arrays
    nAtoms = nu * nv;
    nTriangles = 2 * (nu - 1) * (nv - 1);

    // Allocate the arrays
    xTriMap = new VTuple[nTriangles];
    xBndMap = new size_t[nAtoms + 1];
    
    // Generate the mapping from atom index / boundary side into 
    // the boundary index
    size_t iu, iv, i = 0; xBndMap[0] = 0;
    for(iv = 0; iv < nv; iv++) for(iu = 0; iu < nu; iu++)
      {
      bool isbnd = (iv == 0 || iu == 0 || iu == nu - 1 || iv == nv - 1);
      xBndMap[i+1] = xBndMap[i] + (isbnd ? 1 : 2);
      ++i;
      }

    // Generate inverse boundary mapping
    xInvBndMap = new AtomSidePair[xBndMap[nAtoms]];
    for(size_t i = 0; i < nAtoms; i++)
      {
      xInvBndMap[xBndMap[i + 0] - 0].atom = i;
      xInvBndMap[xBndMap[i + 0] - 0].side = 0;
      xInvBndMap[xBndMap[i + 1] - 1].atom = i;
      xInvBndMap[xBndMap[i + 1] - 1].side = 1;
      }

    // Generate the triangle mapping
    i = 0;
    for(iv = 0; iv < nv - 1; iv++) for(iu = 0; iu < nu - 1; iu++)
      {
      size_t i00 = iv * nu + iu;
      xTriMap[i].v[0] = i00;
      xTriMap[i].v[1] = xTriMap[i+1].v[2] = i00 + 1;
      xTriMap[i].v[2] = xTriMap[i+1].v[1] = i00 + nu;
      xTriMap[i+1].v[0] = i00 + nu + 1;
      i+=2; 
      }

    // The passed in mesh can be stored as the medial mesh
    xMedialMesh = new TriangleMesh();
    TriangleMeshGenerator tgmed(xMedialMesh, this->GetNumberOfAtoms());
    for(MedialTriangleIterator it(this); !it.IsAtEnd(); ++it)
      {
      tgmed.AddTriangle(        
        it.GetAtomIndex(0), it.GetAtomIndex(1), it.GetAtomIndex(2));
      }
    tgmed.GenerateMesh();
    
    // Create a new boundary mesh
    xBoundaryMesh = new TriangleMesh();
    TriangleMeshGenerator tgbnd(xBoundaryMesh, this->GetNumberOfBoundaryPoints());
    for(MedialBoundaryTriangleIterator it(this); !it.IsAtEnd(); ++it)
      {
      tgbnd.AddTriangle(
        it.GetBoundaryIndex(0), it.GetBoundaryIndex(1), it.GetBoundaryIndex(2));
      }
    tgbnd.GenerateMesh();
    };

  ~CartesianGridMedialIterationContext()
    { 
    delete xBndMap; 
    delete xTriMap; 
    delete xBoundaryMesh;
    delete xMedialMesh;
    delete xInvBndMap;
    }

  /**
   * Return the index of at atom at position iu, iv
   */
  size_t GetAtomIndex(size_t iu, size_t iv) const
    { return iv * nu + iu; }

  /** 
   * Get a data structure for the medial triagle mesh. This allows various
   * neighborhood inquiries to be performed, and is good for computing
   * discrete first and second order derivative properties
   */
  TriangleMesh *GetMedialMesh()
    { return this->xMedialMesh; }

  /**
   * Get a data structure for the boundary triangle mesh 
   */
  TriangleMesh *GetBoundaryMesh()
    { return this->xBoundaryMesh; }


private:
  // Dimension of the grid
  size_t nu, nv;

  // Iterable triangle meshes
  TriangleMesh *xMedialMesh, *xBoundaryMesh;

};


#endif
