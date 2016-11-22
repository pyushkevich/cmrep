#ifndef __VTKMeshShortestDistance_h_
#define __VTKMeshShortestDistance_h_

#include <vnl/vnl_vector_fixed.h>

#include <vtkCellLocator.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkFeatureEdges.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkLoopSubdivisionFilter.h>
#include "VTKMeshHalfEdgeWrapper.h"

#include <vector>
#include <list>
#include <utility>

#include "ShortestPath.h"

// Code for dealing with the VTK float->double transition

using namespace std;

template<class T>
inline vnl_vector_fixed<double, 3> CastVTK(T* ptr)
{
  return vnl_vector_fixed<double, 3>(ptr[0], ptr[1], ptr[2]);
}

/** 
 * Hierarchy of function objects used to compute edge weights
 * in the shortest path application
 */
class MeshEdgeWeightFunction 
{
public:
  virtual double GetEdgeWeight(
    vtkPolyData *mesh, vtkIdType x1, vtkIdType x2) const = 0;
  virtual ~MeshEdgeWeightFunction() {}
};

class UnitLengthMeshEdgeWeightFunction : public MeshEdgeWeightFunction
{
public:
  virtual double GetEdgeWeight(
    vtkPolyData *mesh, vtkIdType x1, vtkIdType x2) const
    { return 1.0; }
  virtual ~UnitLengthMeshEdgeWeightFunction() {}
};


class EuclideanDistanceMeshEdgeWeightFunction : public MeshEdgeWeightFunction
{
public:
  virtual double 
    GetEdgeWeight(vtkPolyData *mesh, vtkIdType x1, vtkIdType x2) const
    {
    // Get the two points
    vnl_vector_fixed<double,3> p1 = CastVTK(mesh->GetPoint(x1));
    vnl_vector_fixed<double,3> p2 = CastVTK(mesh->GetPoint(x2));

    // Compute the associated weight
    return (p1 - p2).two_norm();
    }
  virtual ~EuclideanDistanceMeshEdgeWeightFunction() {}
};
  
class PitchBasedMeshEdgeWeightFunction : public MeshEdgeWeightFunction
{
public:
  PitchBasedMeshEdgeWeightFunction() : m_PitchFactor(1.0) {}

  void SetPitchFactor(double c)
    { this->m_PitchFactor = c; };

  double GetPitchFactor() const
    { return m_PitchFactor; }
  
  virtual double 
    GetEdgeWeight(vtkPolyData *mesh, vtkIdType x1, vtkIdType x2) const
    {
    // Get the two points 
    vnl_vector_fixed<double,3> p1 = CastVTK(mesh->GetPoint(x1));
    vnl_vector_fixed<double,3> p2 = CastVTK(mesh->GetPoint(x2));

    // Get the normals associated with the points
    vnl_vector_fixed<double,3> n1 = CastVTK( 
      mesh->GetPointData()->GetNormals()->GetTuple3(x1));
    vnl_vector_fixed<double,3> n2 = CastVTK( 
      mesh->GetPointData()->GetNormals()->GetTuple3(x2));

    // Compute the 'pitch', i.e., the change in the normal w.r.t. step
    // along the edge, as projected on the edge
    double xPitch = fabs(dot_product(n2 - n1,p2 - p1));
    double xDistance = (p2 - p1).two_norm();

    // Compute the associated weight of the edge
    return (m_PitchFactor * xPitch + (1 - m_PitchFactor) * xDistance);
    }
  virtual ~PitchBasedMeshEdgeWeightFunction() {}
private:
  double m_PitchFactor;
};
  
/***************************************************************************
 * This class takes as an input a mesh and computes a graph that can
 * be used to calculate Dijkstra's shortest path on the surface. 
 *
 * The newer version of this class does not preprocess the mesh. Do all the
 * preprocessing on the mesh externally!
 *
 * In addition to providing graph shortest path, this class can also be 
 * used to find the closest vertex to a given point in space and to intersect
 * the graph by a ray. Later this functionality oughtta be moved to another 
 * class
 **************************************************************************/
class VTKMeshShortestDistance 
{
public:
  // type definitions
  typedef vnl_vector_fixed<double,3> Vec;

  /** A callback interface used in conjunction with point checking */
  class ICellChecher { 
  public:
    virtual bool CheckCell(vtkIdType iCell) = 0; 
    virtual ~ICellChecher(){}
  };
  
  /** Constructor */
  VTKMeshShortestDistance();
  ~VTKMeshShortestDistance();

  /** Specify the mesh to use for computing distances */
  void SetInputMesh(VTKMeshHalfEdgeWrapper *mesh);

  /** Specify the edge weight function object to use in order to 
   * compute the costs of edge traversal */
  void SetEdgeWeightFunction(MeshEdgeWeightFunction *function)
    { m_WeightFunctionPtr = function; }

  /** Compute the edge graph */
  void ComputeGraph();

  /** Compute shortest distances from a vertex on a mesh to other vertices */
  void ComputeDistances(
    vtkIdType iStartNode, 
    double xMaxDistance = DijkstraShortestPath<float>::INFINITE_WEIGHT);

  /** Compute the shortest distance from a list of start nodes */
  void ComputeDistances(const list<vtkIdType> &iStartNodes);

  /** Get the distance between start node and given node */
  float GetVertexDistance(vtkIdType iNode) const 
    { return m_ShortestPath->GetDistanceArray()[iNode]; }

  /** Use this to get the path between start node and given node */
  vtkIdType GetVertexPredecessor(vtkIdType iNode) const 
    { return m_ShortestPath->GetPredecessorArray()[iNode]; }

  /** Check if the given node is connected to the source node */
  bool IsVertexConnected(vtkIdType iNode) const
    { 
    return (m_ShortestPath->GetPredecessorArray()[iNode] 
      != DijkstraAlgorithm::NO_PATH);
    }
    
  /** This is a helper method: find vertex whose Euclidean distance to a
    given point is minimal */
  vtkIdType FindClosestVertexInSpace(Vec vec)
    // { return fltLocator->FindClosestPoint(vec[0], vec[1], vec[2]); }
    { return fltLocator->FindClosestPoint(vec.data_block()); }

  /** Find the cell closest to the specified point */
  vtkIdType FindClosestCellInSpace(Vec vec)
    {
    Vec xClosestPoint;
    vtkIdType iCell;
    int subid;

    // VTK Garbage
    vnl_vector_fixed<double, 3> v1(vec[0], vec[1], vec[2]);
    vnl_vector_fixed<double, 3> v2(
      xClosestPoint[0], xClosestPoint[1], xClosestPoint[2]);
    double dist2;
    
    fltCellLocator->FindClosestPoint(
      v1.data_block(), v2.data_block(),iCell,subid,dist2);
    return iCell;
    }

  /** Get the edge mesh to which the indices map */
  vtkPolyData *GetInputMesh() const 
    { return m_SourceMesh; }

  /** Get the weight associated with an edge */
  double GetEdgeWeight(vtkIdType x1, vtkIdType x2) const
    { return m_WeightFunctionPtr->GetEdgeWeight(m_SourceMesh,x1,x2); }
  
  /** Given a ray, find a point on the mesh that is closest to that 
    * and (optinally) satisfies some condition specified by 
    * making a callback to cbCellChecker object */
  bool PickPoint(Vec xStart, Vec xEnd, vtkIdType &point, 
    ICellChecher *cbCellChecher = NULL) const;

  /** Given a ray, find a cell clostest to that ray */
  bool PickCell(Vec xStart, Vec xEnd, vtkIdType &cell) const;

  /** Get the edge weight for a vertex */
  float GetEdgeWeight(unsigned int iEdge) const
    { return m_EdgeWeights[iEdge]; }

private:

  // Clean up graph structures
  void DeleteGraphData();

  // Mesh dimensions
  unsigned int m_NumberOfEdges, m_NumberOfVertices;

  // The array of weights of the graph edges
  float *m_EdgeWeights;

  // The structure used to compute the shortest paths on the mesh
  typedef GraphVoronoiDiagram<float> DijkstraAlgorithm;
  DijkstraAlgorithm *m_ShortestPath;

  // VTK filters
  vtkPointLocator *fltLocator;
  vtkCellLocator *fltCellLocator;

  // VTK source poly-data
  vtkPolyData *m_SourceMesh;
  VTKMeshHalfEdgeWrapper *m_HalfEdge;

  // Function object used to compute edge distances
  MeshEdgeWeightFunction *m_WeightFunctionPtr;

  // Default function object used to compute edge distances
  EuclideanDistanceMeshEdgeWeightFunction m_DefaultWeightFunction;
};

#endif
