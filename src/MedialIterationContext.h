#ifndef __MedialIterationContext_h_
#define __MedialIterationContext_h_

#include <cstdlib>

class TriangleMesh;

class MedialIterationContext
{
public:
  size_t GetNumberOfAtoms() const
    { return nAtoms; }

  size_t GetNumberOfBoundaryPoints() const
    { return xBndMap[nAtoms]; }

  size_t GetNumberOfInternalPoints(size_t ncuts) const
    { return (ncuts + 1) * GetNumberOfBoundaryPoints() + GetNumberOfAtoms(); }

  size_t GetNumberOfTriangles() const
    { return nTriangles; }

  size_t GetNumberOfBoundaryTriangles() const
    { return nTriangles + nTriangles; }

  size_t GetNumberOfCells(size_t ncuts) const
    { return (ncuts + 1) * 2 * GetNumberOfBoundaryTriangles(); }

  size_t GetNumberOfProfileIntervals(size_t ncuts) const
    { return (ncuts + 1) * GetNumberOfBoundaryPoints(); }
   
  bool IsEdgeAtom(size_t i)
    { return xBndMap[i+1] == 1 + xBndMap[i]; }
  
  size_t GetAtomIndexInTriangle(size_t iTri, size_t jVert)
    { return xTriMap[iTri].v[jVert]; }

  size_t GetBoundaryPointIndex(size_t iAtom, size_t iSide)
    { return xBndMap[iAtom + iSide] - iSide; }

  size_t GetInternalPointIndex(size_t iAtom, size_t iSide, size_t iDepth)
    { 
    return iDepth == 0 ? iAtom : 
      nAtoms + (iDepth - 1) * xBndMap[nAtoms] 
      + GetBoundaryPointIndex(iAtom, iSide);
    }
      
  size_t GetProfileIntervalIndex(size_t iBoundary, size_t iDepth)
    { return iDepth * xBndMap[nAtoms] + iBoundary; }

  size_t GetBoundaryPointAtomIndex(size_t iBnd)
    { return xInvBndMap[iBnd].atom; }

  size_t GetBoundaryPointSide(size_t iBnd)
    { return xInvBndMap[iBnd].side; }

  /** 
   * Get a data structure for the medial triagle mesh. This allows various
   * neighborhood inquiries to be performed, and is good for computing
   * discrete first and second order derivative properties
   */
  virtual TriangleMesh *GetMedialMesh() = 0;

  /**
   * Get a data structure for the boundary triangle mesh 
   */
  virtual TriangleMesh *GetBoundaryMesh() = 0;

  virtual ~MedialIterationContext() {}

protected:

  // Can't create an instance of this class
  MedialIterationContext()
    { nAtoms = nTriangles = 0; xBndMap = NULL; xTriMap = NULL; }


  // This is an array of form 0 2 4 5 6 8 ... 
  // difference of 2 between x[i+1] and x[i] indicates i is internal.
  // Then the boundary index for atom i, side j (0,1) is x[i+j]-j
  size_t *xBndMap;

  // This is the inverse of xBndMap
  struct AtomSidePair { size_t atom, side; };
  AtomSidePair *xInvBndMap;

  // This is a mapping from triangles to vertices. It is probably redundant
  // but it makes the implementation a lot easier if we have a copy of this
  // mapping stored here
  struct VTuple { size_t v[3]; } *xTriMap;  

  // Number of atoms and triangles
  size_t nAtoms, nTriangles;
};

#endif
