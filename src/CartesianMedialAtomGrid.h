#ifndef __CartesianMedialAtomGrid_h_
#define __CartesianMedialAtomGrid_h_

#include "MedialAtomGrid.h"

/**
 * An atom grid that lives on a cartesian unit square
 */
class CartesianMedialAtomGrid : public MedialAtomGrid
{
public:
  // Initialize, specifying the number of atoms in x and y directions
  CartesianMedialAtomGrid(size_t m, size_t n);

  size_t GetAtomIndex(size_t i, size_t j)
    { return xIndex[j] + i; }

  size_t GetBoundaryPointIndex(size_t k, size_t iSide)
    { return xBndIndex[(k << 1) + iSide]; }
  
  size_t GetGridSize(size_t d) 
    { return d==0 ? m : n; }

  size_t GetInternalPointIndex(size_t iAtom, size_t iSide, size_t iDepth)
    {
    if(iDepth == 0) return iAtom;
    return numAtoms + (iDepth - 1) * numBndPts + GetBoundaryPointIndex(iAtom, iSide);
    }

  size_t GetProfileIntervalIndex(size_t iBoundary, size_t iDepth)
    {
    return iDepth * numBndPts + iBoundary;
    }

  // Create various iterators
  MedialAtomIterator *NewAtomIterator();
  MedialQuadIterator *NewQuadIterator();
  MedialBoundaryPointIterator *NewBoundaryPointIterator();
  MedialBoundaryQuadIterator *NewBoundaryQuadIterator();
  MedialInternalCellIterator *NewInternalCellIterator(size_t nCuts);
  MedialInternalPointIterator *NewInternalPointIterator(size_t nCuts);
  MedialProfileIntervalIterator *NewProfileIntervalIterator(size_t nCuts);

  // Get the number of iterated objects 
  size_t GetNumberOfAtoms()
    { return numAtoms; }

  size_t GetNumberOfQuads()
    { return numQuads; }

  size_t GetNumberOfBoundaryPoints()
    { return numBndPts; }

  size_t GetNumberOfBoundaryQuads()
    { return numBndQuads; }

  size_t GetNumberOfCells(size_t nCuts)
    { return (nCuts + 1) * 2 * numQuads; }

  size_t GetNumberOfInternalPoints(size_t nCuts)
    { return (nCuts + 1) * GetNumberOfBoundaryPoints() + GetNumberOfAtoms(); }

  size_t GetNumberOfProfileIntervals(size_t nCuts)
    { return (nCuts + 1) * GetNumberOfBoundaryPoints(); }

private:
  // A way to map a pair of indices to a linear index
  vector<size_t> xIndex, xBndIndex;

  // Dimensions of the grid and various sizes
  size_t m, n, nEdge, nInner;
  size_t numAtoms, numQuads, numBndQuads, numBndPts;
};


void TestCartesianGrid();

#endif
