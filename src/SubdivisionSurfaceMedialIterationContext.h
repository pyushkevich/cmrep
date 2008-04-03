#ifndef __SubdivisionSurfaceMedialIterationContext_h_
#define __SubdivisionSurfaceMedialIterationContext_h_

#include "SubdivisionSurface.h"
#include "MedialIterationContext.h"

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

    // Generate the triangle mapping
    for(size_t j = 0; j < nTriangles; j++)
      for(size_t k = 0; k < 3; k++)
        xTriMap[j].v[k] = inMesh->triangles[j].vertices[k];
    }

  ~SubdivisionSurfaceMedialIterationContext()
    {
    delete xTriMap;
    delete xBndMap;
    }
};

#endif
