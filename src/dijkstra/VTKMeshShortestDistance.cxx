#include "VTKMeshShortestDistance.h"
#include "vtkIdList.h"

VTKMeshShortestDistance
::VTKMeshShortestDistance()
{
  // Set the distance function
  m_WeightFunctionPtr = &m_DefaultWeightFunction;

  // Initialize the filters
  fltLocator = vtkPointLocator::New();
  fltCellLocator = vtkCellLocator::New();

  // Initialize the graph to NULL
  m_ShortestPath = NULL;
  m_EdgeWeights = NULL;
};

void
VTKMeshShortestDistance
::SetInputMesh(VTKMeshHalfEdgeWrapper *wrapper)
{
  // Store the input mesh
  m_HalfEdge = wrapper;
  m_SourceMesh = wrapper->GetPolyData();

  // Construct a locator
  fltLocator->SetDataSet(m_SourceMesh);
  fltLocator->BuildLocator();

  // Construct the cell locator
  fltCellLocator->SetDataSet(m_SourceMesh);
  fltCellLocator->BuildLocator();

  // Set the number of vertices
  m_NumberOfVertices = m_SourceMesh->GetNumberOfPoints();
  m_NumberOfEdges = m_HalfEdge->GetNumberOfHalfEdges();
}

void
VTKMeshShortestDistance
::ComputeGraph()
{
  // Clean up the old graph info
  DeleteGraphData();

  // Allocate the graph data
  m_EdgeWeights = new float[m_NumberOfEdges];

  // Compute the length of each half-edge in the graph
  for(unsigned int i=0;i<m_NumberOfEdges;i++)
    {
    vtkIdType iHead = m_HalfEdge->GetHalfEdgeVertex(i);
    vtkIdType iTail = m_HalfEdge->GetHalfEdgeTailVertex(i);

    // Compute the edge weight
    m_EdgeWeights[i] = 
      (float) m_WeightFunctionPtr->GetEdgeWeight(m_SourceMesh, iHead, iTail);
    }

  // Create the shortest path object
  m_ShortestPath = new DijkstraAlgorithm( 
    m_NumberOfVertices, 
    m_HalfEdge->GetAdjacencyIndex(), 
    m_HalfEdge->GetAdjacency(), 
    m_EdgeWeights);
}

VTKMeshShortestDistance
::~VTKMeshShortestDistance()
{
  DeleteGraphData();
  fltLocator->Delete();
  fltCellLocator->Delete();
}

bool
VTKMeshShortestDistance
::PickCell(Vec xStart, Vec xEnd, vtkIdType &point) const
{
  // Ugly VTK
  double v1[3], v2[3], ptLine[3], pCoords[3];
  int subId;
  vtkIdType cellid;
  double t; 

  v1[0] = xStart[0]; v1[1] = xStart[1]; v1[2] = xStart[2];
  v2[0] = xEnd[0]; v2[1] = xEnd[1]; v2[2] = xEnd[2];

  // Compute the intersection with the line
  fltCellLocator->IntersectWithLine( 
    v1,v2, 0.001, t, ptLine, pCoords, subId, cellid);

  return subId == 0;
}

bool 
VTKMeshShortestDistance
::PickPoint(Vec xStart, Vec xEnd, vtkIdType &point, ICellChecher *cbCell) const
{
  double v1[3], v2[3], ptLine[3], pCoords[3];
  double t; 
  int subId; vtkIdType iCell;
  
  v1[0] = xStart[0]; v1[1] = xStart[1]; v1[2] = xStart[2];
  v2[0] = xEnd[0]; v2[1] = xEnd[1]; v2[2] = xEnd[2];

  do
    {
    // cout << "Searching ray " << xStart << "  to  " << xEnd << endl;
    
    // Compute the intersection with the line
    int rc = fltCellLocator->IntersectWithLine(
      v1, v2, 0.001, t, ptLine, pCoords, subId, iCell);

    // cout << "   RC: " << rc << " \t T = " << t << " \t iCell " << iCell << endl;

    // If no intersection found, return false
    if(!rc) return false;

    // Increase the starting vector by t + epsilon
    xStart += (t + 0.001) * (xEnd - xStart);
    }
  while(cbCell && !cbCell->CheckCell(iCell));

  // cout << "Tracing from " << xStart << " to " << xEnd << endl;
  // cout << "Intersection at t = " << t << ", ptline " << ptLine << endl;
  // cout << "Vertex ID = " << subId << endl;

  // Find the vertex closest to the intersection
  point = fltLocator->FindClosestPoint(ptLine);
  return subId == 0;
}

void 
VTKMeshShortestDistance
::ComputeDistances(const list<vtkIdType> &lSources)
{
  // Create an array from the list
  unsigned int *lSourceArray = new unsigned int[lSources.size()];
  unsigned int iSource = 0;

  list<vtkIdType>::const_iterator it = lSources.begin();
  while(it != lSources.end())
    lSourceArray[iSource++] = *it++;

  m_ShortestPath->ComputePathsFromManySources(iSource, lSourceArray);
  
  delete[] lSourceArray;
}

void 
VTKMeshShortestDistance
::ComputeDistances(vtkIdType iStartNode, double xMaxDistance)
{
  m_ShortestPath->ComputePathsFromSource(iStartNode, xMaxDistance);
}

void 
VTKMeshShortestDistance
::DeleteGraphData()
{
  if(m_ShortestPath)
    {
    delete m_ShortestPath; m_ShortestPath = NULL;
    delete m_EdgeWeights;
    }
}
