#ifndef __SubdivisionSurfaceMedialIterationContext_h_
#define __SubdivisionSurfaceMedialIterationContext_h_

#include "SubdivisionSurface.h"
#include "MedialIterationContext.h"
#include "MeshTraversal.h"

class SubdivisionSurfaceMedialIterationContext 
: public MedialIterationContext
{
public:
  typedef SubdivisionSurface::MeshLevel MeshLevel;

  SubdivisionSurfaceMedialIterationContext(MeshLevel *inMesh)
    {
    // Set the array sizes
    nAtoms = inMesh->nVertices;
    nTriangles = inMesh->triangles.size();
    
    // Allocate the arrays
    xTriMap = new VTuple[nTriangles];
    xBndMap = new size_t[nAtoms + 1];
    
    // Generate the mapping from atom index / boundary side into 
    // the boundary index
    xBndMap[0] = 0;
    for(size_t i = 0; i < nAtoms; i++)
      xBndMap[i+1] = xBndMap[i] + (inMesh->IsVertexInternal(i) ? 2 : 1);   

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
    for(size_t j = 0; j < nTriangles; j++)
      for(size_t k = 0; k < 3; k++)
        xTriMap[j].v[k] = inMesh->triangles[j].vertices[k];

    // The passed in mesh can be stored as the medial mesh
    xMedialMesh = inMesh;
    
    // Create a new boundary mesh
    xBoundaryMesh = new TriangleMesh();

    // Generate the mesh from the triangle information
    TriangleMeshGenerator tgbnd(xBoundaryMesh, this->GetNumberOfBoundaryPoints());

    // Iterate over the boundary triangles, adding points
    for(MedialBoundaryTriangleIterator it(this); !it.IsAtEnd(); ++it)
      {
      tgbnd.AddTriangle(
        it.GetBoundaryIndex(0), it.GetBoundaryIndex(1), it.GetBoundaryIndex(2));
      }

    // Generate it!
    tgbnd.GenerateMesh();
    }

  ~SubdivisionSurfaceMedialIterationContext()
    {
    delete xTriMap;
    delete xBndMap;
    delete xBoundaryMesh;
    delete xInvBndMap; 
    }

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
  TriangleMesh *xMedialMesh, *xBoundaryMesh;
};

#endif
