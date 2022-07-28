#ifndef _VTKMeshHalfEdgeWrapper_h_
#define _VTKMeshHalfEdgeWrapper_h_

#include "vtkPolyData.h"

class VTKMeshHalfEdgeWrapper
{
public:
  VTKMeshHalfEdgeWrapper(vtkPolyData *mesh)
    {
    // Prepare the mesh
    xMesh = mesh;
    xMesh->BuildCells();
    xMesh->BuildLinks();
    
    // Initialize the index data structure
    nVertices = mesh->GetNumberOfPoints();
    nFaces = mesh->GetNumberOfCells();
    xAdjacencyIndex = new unsigned int[nVertices+1];
    xAdjacencyIndex[0] = 0;

    // Count the number of edges leaving each vertex
    for(unsigned int iVtx=0; iVtx < nVertices; iVtx++)
      {
      // Get the number of cells for this point
      vtkIdType nCells;
      vtkIdType *xDummy = NULL;
      xMesh->GetPointCells(iVtx, nCells, xDummy);

      // Add to the adjacency index
      xAdjacencyIndex[iVtx+1] = xAdjacencyIndex[iVtx] + nCells;
      }

    // Set the number of half-edges
    nHalfEdges = xAdjacencyIndex[nVertices];

    // Allocate the edge-related structures
    xAdjacency = new unsigned int[nHalfEdges];
    xFace = new unsigned int[nHalfEdges];
    xFlipEdge = new unsigned int[nHalfEdges];
    xNextEdge = new unsigned int[nHalfEdges];
    xFaceEdges = new unsigned int[nFaces];

    // Allocate an additional array that keeps track of how many half-edges
    // have been added for each vertex
    unsigned int *xTemp = new unsigned int[nVertices];
    memset(xTemp, 0, sizeof(unsigned int) * nVertices);

    // Traverse all the cells in the VTK mesh. The assumption here is that each cell
    // is traversed in the consistant, counter-clockwise order.
    for(unsigned int iCell = 0; iCell < (unsigned int) xMesh->GetNumberOfCells(); iCell++)
      {
      // Get the points for this cell
      vtkIdType nPoints;
      const vtkIdType *xPoints;
      xMesh->GetCellPoints(iCell, nPoints, xPoints);

      // Walk around the list of points
      for(unsigned int j = 0; j < (unsigned int) nPoints; j++)
        {
        // Get the head and the tail of the current half-edge
        unsigned int iTail = xPoints[j], iHead = xPoints[(j+1) % nPoints];

        // Get the index of the current half-edge
        unsigned int index = xAdjacencyIndex[iTail] + xTemp[iTail];

        // Set the head of the current half-edge
        xAdjacency[index] = iHead;

        // Set the face corresponding to the current half-edge
        xFace[index] = iCell;
        xFaceEdges[iCell] = index;

        // Set the next half edge
        xNextEdge[index] = xAdjacencyIndex[iHead] + xTemp[iHead];
        }

      // Increment all the edge counters
      for(unsigned int k = 0;  k < (unsigned int) nPoints; k++)
        xTemp[xPoints[k]]++;
      }

    // Check the consistency of the temp array
    for(unsigned int iTest=0; iTest < nVertices; iTest++)
      if(xTemp[iTest] != xAdjacencyIndex[iTest+1] - xAdjacencyIndex[iTest])
        throw "Consistency check failed in VTKMeshHalfEdgeWrapper::VTKMeshHalfEdgeWrapper";

    // Clean up the temp array
    delete[] xTemp;

    // Now we should have created all the half-edges. What's left is to establish
    // the symmetry links between them
    unsigned int iTail, iHead, iEdge, jEdge;
    for(iTail = 0; iTail < nVertices; iTail++)
      {
      for(iEdge = xAdjacencyIndex[iTail]; iEdge < xAdjacencyIndex[iTail+1]; iEdge++)
        {
        iHead = xAdjacency[iEdge];
        for(jEdge = xAdjacencyIndex[iHead]; jEdge < xAdjacencyIndex[iHead+1]; jEdge++)
          {
          if(xAdjacency[jEdge] == iTail)
            xFlipEdge[iEdge] = jEdge;
          }
        }
      }
    }

  ~VTKMeshHalfEdgeWrapper()
    {
    delete xAdjacency;
    delete xAdjacencyIndex;
    delete xFlipEdge;
    delete xFace;
    delete xNextEdge;
    delete xFaceEdges;
    }

  vtkPolyData *GetPolyData() const
    { return xMesh; }

  /** Get the number of vertices in the graph */
  unsigned int GetNumberOfVertices() const 
    { return nVertices; }

  /** Get the number of edges in the graph. This returns the number of
   * bidirectional edges, which is twice the number of directional edges.
   * This way it is possible to have different weights for different 
   * directions through the edge */
  unsigned int GetNumberOfHalfEdges() const
    { return nHalfEdges; }

  /** Get number of vertices adjacent to a given vertex */
  unsigned int GetVertexNumberOfEdges(unsigned int iVertex) const
    { return xAdjacencyIndex[iVertex+1] - xAdjacencyIndex[iVertex]; }

  /** Get the given edge at a vertex */
  unsigned int GetVertexHalfEdge(unsigned int iVertex, unsigned int iEdge) const
    { return xAdjacencyIndex[iVertex] + iEdge; }

  /** Get the vertex at the head of the n-th halfedge */
  unsigned int GetHalfEdgeVertex(unsigned int iHalf) const
    { return xAdjacency[iHalf]; }

  /** Get the vertex at the head of the n-th halfedge */
  unsigned int GetHalfEdgeTailVertex(unsigned int iHalf) const
    { return xAdjacency[xFlipEdge[iHalf]]; }

  /** Get the flip edge of a given edge */
  unsigned int GetHalfEdgeOpposite(unsigned int iHalf) const
    { return xFlipEdge[iHalf]; }

  /** Get the edge following a given edge */
  unsigned int GetHalfEdgeNext(unsigned int iHalf) const
    { return xNextEdge[iHalf]; }

  /** Get the face to the left of a given edge */
  unsigned int GetHalfEdgeFace(unsigned int iHalf) const
    { return xFace[iHalf]; }

  /** Get one of the half-edges assigned to a face */
  unsigned int GetFaceHalfEdge(unsigned int iFace) const
    { return xFaceEdges[iFace]; }

  /** Get the edge between two vertices (if any) */
  bool GetHalfEdgeBetweenVertices(unsigned int iTail, unsigned int iHead, unsigned int &iEdge)
    {
    for(iEdge = xAdjacencyIndex[iTail]; iEdge < xAdjacencyIndex[iTail+1]; iEdge++)
      if(xAdjacency[iEdge] == iHead)
        return true;
    return false;
    }

  /** Get the adjacency index (for METIS, etc) */
  unsigned int *GetAdjacencyIndex() const
    { return xAdjacencyIndex; }

  /** Get the adjacency array (for METIS, etc) */
  unsigned int *GetAdjacency() const
    { return xAdjacency; }

protected:
  // Reference to the wrapped mesh
  vtkPolyData *xMesh;

  // Number of vertices and half-edges in the mesh
  unsigned int nVertices, nHalfEdges, nFaces;

  // Position in the half-edge array of the half-edges exiting each vertex
  unsigned int *xAdjacencyIndex;

  // List of vertices pointed at by each half-edge
  unsigned int *xAdjacency;

  // List of faces to the left of each half-edge
  unsigned int *xFace;

  // List of complement half-edges
  unsigned int *xFlipEdge;

  // List of subsequent half-edges around the face
  unsigned int *xNextEdge;

  // An arbitrary half-edge for each face
  unsigned int *xFaceEdges;

};

#endif // _VTKMeshHalfEdgeWrapper_h_
