#ifndef __MedialAtomIterators_h_
#define __MedialAtomIterators_h_

#include "MedialIterationContext.h"

/**
 * An abstract iterator for parsing points in the medial atom grid
 */
class MedialAtomIterator
{
public:
  MedialAtomIterator(MedialIterationContext *context)
    { this->context = context; i = 0; }
  
  size_t GetIndex() const
    { return i; }
  
  bool IsEdgeAtom() const
    { return context->IsEdgeAtom(i); }

  MedialAtomIterator &operator++()
    { ++i; return *this; }
  
  bool IsAtEnd() const
    { return i >= context->GetNumberOfAtoms(); }
  
  void GoToBegin()
    { i = 0; }

private:
  // Context object defining our behavior
  MedialIterationContext *context;

  // Atom index
  size_t i;  
};

/**
 * An iterator for parsing quad cells in a medial atom grid
 */
class MedialTriangleIterator
{
public:
  MedialTriangleIterator(MedialIterationContext *context)
    { this->context = context; i = 0; }
  
  /** Get the index of the atom at the j-th vertex of the triangle. */
  size_t GetAtomIndex(size_t j) const
    { return context->GetAtomIndexInTriangle(i, j); }

  MedialTriangleIterator &operator++()
    { i++; return *this; }

  bool IsAtEnd() const
    { return i >= context->GetNumberOfTriangles(); }
  
  void GoToBegin()
    { i = 0; }

  size_t GetIndex() const
    { return i; }

private:
  // Context object defining our behavior
  MedialIterationContext *context;

  // Triangle index in the context
  size_t i;
};

/**
 * An iterator for parsing boundary points associated with a medial grid. This
 * iterator will visit every boundary point exactly once.
 */
class MedialBoundaryPointIterator
{
public:
  MedialBoundaryPointIterator(MedialIterationContext *context)
    : itAtom(context)
    { this->context = context; this->GoToBegin(); }

  size_t GetIndex() const
    { return context->GetBoundaryPointIndex(itAtom.GetIndex(), iSide); }

  size_t GetOppositeIndex() const
    { return context->GetBoundaryPointIndex(itAtom.GetIndex(), 1 - iSide); }
  
  size_t GetAtomIndex() const
    { return itAtom.GetIndex(); }

  bool IsEdgeAtom() const
    { return itAtom.IsEdgeAtom(); }

  size_t GetBoundarySide() const
    { return iSide; }

  MedialBoundaryPointIterator &operator++() 
    {
    if(itAtom.IsEdgeAtom() || iSide == 1)
      { ++itAtom; iSide = 0; }
    else
      iSide = 1;
    return *this;
    }
  
  bool IsAtEnd() const
    { return itAtom.IsAtEnd(); }
  
  void GoToBegin() 
    { itAtom.GoToBegin(); iSide = 0; }

private:
  // Context object defining our behavior
  MedialIterationContext *context;

  // Atom iterator
  MedialAtomIterator itAtom;

  // Triangle index in the context
  size_t iSide;
};

/**
 * An iterator for parsing quad cells in a medial atom grid
 */
class MedialBoundaryTriangleIterator
{
public:
  MedialBoundaryTriangleIterator(MedialIterationContext *context)
    : itTriangle(context)
    { this->context = context; this->GoToBegin(); }
  
  /** Get the boundary index corresponding to the j-th vertex in this
   * triangle. The loop 0, 1, 2 is guaranteed to be clockwise */
  size_t GetBoundaryIndex(size_t j) const
    { return context->GetBoundaryPointIndex(GetAtomIndex(j), iSide); }

  /** Get the atom index corresponding to the j-th vertex in the triangle */
  size_t GetAtomIndex(size_t j) const
    { return itTriangle.GetAtomIndex( iSide == 1 ? j : 2 - j); }

  /** Get the index of the corresponding triangle in the medial mesh */
  size_t GetMedialTriangleIndex() const
    { return itTriangle.GetIndex(); }

  size_t GetBoundarySide() const
    { return iSide; }

  size_t GetIndex() const
    { return iSide + (itTriangle.GetIndex() << 1); }

  MedialBoundaryTriangleIterator &operator++()
    {
    if(iSide == 0) 
      { iSide = 1; }
    else 
      { iSide = 0; ++itTriangle; }
    return *this;
    }
  
  bool IsAtEnd() const
    { return itTriangle.IsAtEnd(); }
  
  void GoToBegin()
    { itTriangle.GoToBegin(); iSide = 0; }
  
private:
  MedialIterationContext *context;
  MedialTriangleIterator itTriangle;
  size_t iSide;
};

/** 
 * An iterator that goes through the internal points in an m-rep. This
 * iterator first steps over all the atoms on the medial surface (points at
 * depth 0) and then for each depth level it steps over all the boundary
 * points 
 */
class MedialInternalPointIterator
{
public:
  MedialInternalPointIterator(MedialIterationContext *context, size_t nCuts)
    : itAtom(context), itBnd(context)
    {
    this->context = context;
    // Max depth is nCuts + 1, i.e. point index counting from 0 (medial atom)
    // to -maxdepth or +maxdepth
    iMaxDepth = nCuts + 1;
    GoToBegin();
    }
    
  // Get the index of this internal point in the parent context
  size_t GetIndex() const
    { 
    // assert(iIndex == context->GetInternalPointIndex(
    //   GetAtomIndex(), GetBoundarySide(), iDepth))
    return iIndex;
    }
  
  // Get the index of the corresponding atom on the medial surface
  size_t GetAtomIndex() const
    { return iDepth == 0 ? itAtom.GetIndex() : itBnd.GetAtomIndex(); }

  // Get the corresponding boundary index (-1 returned for depth=0)
  size_t GetBoundaryIndex() const
    { return iDepth == 0 ? (size_t) -1 : itBnd.GetIndex(); }

  // Is the point on the profile that corresponds to a medial edge?
  bool IsEdgeAtom() const
    { return iDepth == 0 ? itAtom.IsEdgeAtom() : itBnd.IsEdgeAtom(); }

  // Get a non-negative value t that is zero on medial axis, 1 on boundary
  double GetRelativeDistanceToMedialAxis() const
    { return xRelativeDepth; }

  // This is only defined for points that are not on the medial axis and 
  // not on the crest (edge of the medial surface)
  size_t GetBoundarySide() const
    { return iDepth == 0 ? 0 : itBnd.GetBoundarySide(); }
  
  // Get the maximum depth (corresponding to the boundary)
  size_t GetMaxDepth()  const
    { return iMaxDepth; }
  
  // Get the unsigned depth of the point; 0 means on the medial axis and
  // MaxDepth means we are at the level of the boundary
  size_t GetDepth()  const
    { return iDepth; } 

  // We first iterate over the medial atoms, then over the atoms closest to
  // the medial axis, and so on, until we get to the boundary
  MedialInternalPointIterator &operator++()
    {
    // Go to the next internal point
    if(iDepth == 0)
      {
      // At depth 0, iterate over the medial atoms
      ++itAtom;
      if(itAtom.IsAtEnd())
        SetDepth(1);
      }
    else
      {
      // At other depth, iterate over the given layer
      ++itBnd;
      if(itBnd.IsAtEnd())
        { SetDepth(iDepth+1); itBnd.GoToBegin(); }
      }
    iIndex++;
    return *this;
    }
  
  bool IsAtEnd() const
    { return iDepth > iMaxDepth; } 
  
  void GoToBegin()
    {
    itAtom.GoToBegin(); itBnd.GoToBegin();
    SetDepth(0);
    iIndex = 0; 
    }
  
private:
  // Internal function to update the depth
  void SetDepth(size_t iNewDepth) 
    {
    iDepth = iNewDepth;
    if(iDepth == 0) xRelativeDepth = 0.0;
    else if(iDepth == iMaxDepth) xRelativeDepth = 1.0;
    else xRelativeDepth = iDepth * 1.0 / iMaxDepth;
    }
  
  MedialIterationContext *context;
  MedialAtomIterator itAtom;
  MedialBoundaryPointIterator itBnd;
  size_t iDepth, iMaxDepth;
  double xRelativeDepth;
  size_t iIndex;
};

/** 
 * An iterator that goes through the internal cells in an m-rep. These cells
 * are shaped tike Toblerone bars: triangle in one cross-section and a
 * rectangle in the other cross-section. The triangles are 'parallel' to the
 * medial axis and boundary. 
 */
class MedialInternalCellIterator
{
public:
  MedialInternalCellIterator(MedialIterationContext *context, size_t nCuts) 
    : itTriangle(context)
    {
    this->context = context;
    iMaxDepth = nCuts;
    GoToBegin();
    }
    
  // Get the index corresponding to the j-th vertex of the triangle defining
  // the cross-section of the cell (j = 0, 1, 2)
  size_t GetAtomIndex(size_t j) const
    { return itTriangle.GetAtomIndex(j); }

  // Return the internal point index of one of the six vertices of the cell.
  // The first index is the corner of the triangle face (see GetAtomIndex) and
  // the second is the index of the face relative to the medial axis (0 =
  // nearer to the medial axis, 1 = farther from it) 
  size_t GetInternalPointIndex(size_t jVertex, size_t jDepth) const
    { 
    return context->GetInternalPointIndex(
      itTriangle.GetAtomIndex(jVertex), 
      itTriangle.GetBoundarySide(), iDepth + jDepth); 
    }

  // Get the index of the profile interval at one of three toblerone edges
  size_t GetProfileIntervalIndex(size_t j) const
    {
    return context->GetProfileIntervalIndex(
      itTriangle.GetBoundaryIndex(j), iDepth);
    }
    
  /** Get the depth of the side k (0 inner, 1 outer) */
  size_t GetDepth(size_t k) const
    { return iDepth + k; }
  
  MedialInternalCellIterator &operator++()
    {
    ++itTriangle;
    if(itTriangle.IsAtEnd())
      { ++iDepth; itTriangle.GoToBegin(); }
    return *this;
    }

  bool IsAtEnd() const
    { return iDepth > iMaxDepth; }
  
  void GoToBegin()
    {
    iDepth = 0;
    itTriangle.GoToBegin();
    }

private:
  MedialIterationContext *context;
  MedialBoundaryTriangleIterator itTriangle;
  size_t iDepth, iMaxDepth;
};

/** An iterator that goes through point pairs on the profiles that connect the
 * medial axis to the boundary. There are four such intervals for every
 * internal cell */
class MedialProfileIntervalIterator
{
public:
  MedialProfileIntervalIterator(MedialIterationContext *context, size_t nCuts) 
    : itBnd(context)
    {
    iMaxDepth = nCuts;
    this->context = context;
    GoToBegin();
    }
  
  size_t GetInnerPointIndex() const
    { 
    return context->GetInternalPointIndex(
      itBnd.GetAtomIndex(), itBnd.GetBoundarySide(), iDepth);
    }
  
  size_t GetOuterPointIndex() const
    { 
    return context->GetInternalPointIndex(
      itBnd.GetAtomIndex(), itBnd.GetBoundarySide(), iDepth + 1);
    }
  
  size_t GetIndex() const
    { return context->GetProfileIntervalIndex(itBnd.GetIndex(), iDepth); }

  MedialProfileIntervalIterator &operator++()
    {
    ++itBnd;
    if(itBnd.IsAtEnd())
      {
      iDepth++;
      itBnd.GoToBegin();
      }
    return *this;
    }
  
  bool IsAtEnd() const
    { return iDepth > iMaxDepth; }

  void GoToBegin()
    {
    itBnd.GoToBegin();
    iDepth = 0;
    }

private:
  MedialBoundaryPointIterator itBnd;
  MedialIterationContext *context;
  size_t iMaxDepth, iDepth;
};

#endif
