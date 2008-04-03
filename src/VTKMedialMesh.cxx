#include "CodeTimer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnstructuredGrid.h"
#include "MedialAtom.h"
#include "MedialAtomGrid.h"
#include "GenericMedialModel.h"
#include "ITKImageWrapper.h"
#include "itkImage.h"
#include "OptimizationTerms.h"

vtkFloatArray *AddMedialScalarField(vtkPolyData *target, GenericMedialModel *model, char *name)
{
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetName(name);
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(model->GetNumberOfAtoms());
  target->GetPointData()->AddArray(array);
  return array;
}

vtkFloatArray *AddBndScalarField(vtkPolyData *target, GenericMedialModel *model, char *name)
{
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetName(name);
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(model->GetNumberOfBoundaryPoints());
  target->GetPointData()->AddArray(array);
  return array;
}

vtkFloatArray *AddBndTriScalarField(vtkPolyData *target, GenericMedialModel *model, char *name)
{
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetName(name);
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(model->GetNumberOfBoundaryTriangles());
  target->GetCellData()->AddArray(array);
  return array;
}

vtkFloatArray *AddMedialVectorField(vtkPolyData *target, GenericMedialModel *model, char *name)
{
  vtkFloatArray *array = vtkFloatArray::New();
  array->SetName(name);
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(model->GetNumberOfAtoms());
  target->GetPointData()->AddArray(array);
  return array;
}

void ExportMedialMeshToVTK(
  GenericMedialModel *xModel,
  ITKImageWrapper<float> *xImage,
  const char *file)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xModel->GetNumberOfAtoms());
  
  // Allocate the polydata
  vtkPolyData *pMedial = vtkPolyData::New();
  pMedial->Allocate(xModel->GetNumberOfTriangles());
  pMedial->SetPoints(lPoints);

  // Allocate the array of normals
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(xModel->GetNumberOfAtoms());
  pMedial->GetPointData()->SetNormals(lNormals);

  // Allocate the scalar arrays
  vtkFloatArray *lMetric = AddMedialScalarField(pMedial, xModel, "Covariant Tensor Determinant");
  vtkFloatArray *lRho = AddMedialScalarField(pMedial, xModel, "Rho Function");
  vtkFloatArray *lRadius = AddMedialScalarField(pMedial, xModel, "Radius Function");
  vtkFloatArray *lDummy1 = AddMedialScalarField(pMedial, xModel, "Dummy1");
  vtkFloatArray *lBending = AddMedialScalarField(pMedial, xModel, "Bending Energy");
  vtkFloatArray *lRegularity = AddMedialScalarField(pMedial, xModel, "Regularity Penalty");
  vtkFloatArray *lAngle = AddMedialScalarField(pMedial, xModel, "Metric Angle");
  vtkFloatArray *lCoordU = AddMedialScalarField(pMedial, xModel, "U Coordinate");
  vtkFloatArray *lCoordV = AddMedialScalarField(pMedial, xModel, "V Coordinate");
  vtkFloatArray *lMeanCurv = AddMedialScalarField(pMedial, xModel, "Mean Curvature");
  vtkFloatArray *lGaussCurv = AddMedialScalarField(pMedial, xModel, "Gauss Curvature");
  vtkFloatArray *lKappa1 = AddMedialScalarField(pMedial, xModel, "Kappa1");
  vtkFloatArray *lKappa2 = AddMedialScalarField(pMedial, xModel, "Kappa2");
  vtkFloatArray *lNormal = AddMedialVectorField(pMedial, xModel, "Atom Normal");
  vtkFloatArray *lStretch = AddMedialScalarField(pMedial, xModel, "Stretch");
  vtkFloatArray *lCurvPen = AddMedialScalarField(pMedial, xModel, "Curvature Penalty Feature");
  vtkFloatArray *lAreaElement = AddMedialScalarField(pMedial, xModel, "Area Element");

  vtkFloatArray *lSpoke1 = AddMedialVectorField(pMedial, xModel, "Spoke1");
  vtkFloatArray *lSpoke2 = AddMedialVectorField(pMedial, xModel, "Spoke2");

  vtkFloatArray *lContraOffDiag = 
    AddMedialScalarField(pMedial, xModel, "Off Diagonal Term of Contravariant MT");

  vtkFloatArray *lXu = AddMedialVectorField(pMedial, xModel, "Xu");
  vtkFloatArray *lXv = AddMedialVectorField(pMedial, xModel, "Xv");

  // Allocate and add the image intensity array
  bool flagImage = xImage && xImage->IsImageLoaded();
  vtkFloatArray *lImage = vtkFloatArray::New();
  if(flagImage)
    {
    lImage->SetName("Image");
    lImage->SetNumberOfComponents(1);
    lImage->SetNumberOfTuples(xModel->GetNumberOfAtoms());
    pMedial->GetPointData()->AddArray(lImage);
    }
  
  // Get the internals of the medial model
  MedialIterationContext *xGrid = xModel->GetIterationContext();
  MedialAtom *xAtoms = xModel->GetAtomArray();

  // Add all the points
  for(MedialAtomIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    MedialAtom &a = xAtoms[i];

    SMLVec3d X = a.X; SMLVec3d N = a.N;
    lPoints->InsertNextPoint(X[0], X[1], X[2]);
    lNormals->SetTuple3(i, N[0], N[1], N[2]);
    lXu->SetTuple3(i, a.Xu[0], a.Xu[1], a.Xu[2]);
    lXv->SetTuple3(i, a.Xv[0], a.Xv[1], a.Xv[2]);
    lMetric->SetTuple1(i, a.G.g);
    lRadius->SetTuple1(i, a.R);
    lRho->SetTuple1(i, a.xLapR);

    lCoordU->SetTuple1(i, a.u);
    lCoordV->SetTuple1(i, a.v);

    lMeanCurv->SetTuple1(i, a.xMeanCurv);
    lGaussCurv->SetTuple1(i, a.xGaussCurv);
    lKappa1->SetTuple1(i, a.xMeanCurv - sqrt( a.xMeanCurv *  a.xMeanCurv -  a.xGaussCurv));
    lKappa2->SetTuple1(i, a.xMeanCurv + sqrt( a.xMeanCurv *  a.xMeanCurv -  a.xGaussCurv));

    lAreaElement->SetTuple1(i, a.aelt);
    
    lNormal->SetTuple3(i, a.N(0), a.N(1), a.N(2));

    // Compute the stretch ???
    lStretch->SetTuple1(i, a.G.xChristoffelSecond[0][0][0]);

    // Add the spoke vectors
    lSpoke1->SetTuple3(i, a.R * a.xBnd[0].N[0], a.R * a.xBnd[0].N[1], a.R * a.xBnd[0].N[2]);
    lSpoke2->SetTuple3(i, a.R * a.xBnd[1].N[0], a.R * a.xBnd[1].N[1], a.R * a.xBnd[1].N[2]);

    // Compute sum of the squares of principal curvatures
    double k2 = 4 * a.xMeanCurv * a.xMeanCurv - 2 * a.xGaussCurv;
    lCurvPen->SetTuple1(i, k2);

    // Set the bending energy
    lBending->SetTuple1(i,
      dot_product(a.Xuu, a.Xuu) + dot_product(a.Xvv,a.Xvv) + 2.0 * dot_product(a.Xuv, a.Xuv));

    // Set the regularity energy
    double reg1 = a.G.xChristoffelSecond[0][0][0] + a.G.xChristoffelSecond[1][0][1];
    double reg2 = a.G.xChristoffelSecond[0][1][0] + a.G.xChristoffelSecond[1][1][1];
    double reg = reg1 * reg1 + reg2 * reg2;
    lRegularity->SetTuple1(i, reg);

    // Set the angle between Xu and Xv
    double dp = dot_product(a.Xu, a.Xv);

    lAngle->SetTuple1(i,
      (dp * dp) / (a.Xu.squared_magnitude() * a.Xv.squared_magnitude()));
    lContraOffDiag->SetTuple1(i, 
      a.G.xCovariantTensor[0][1] * a.G.xCovariantTensor[0][1] / a.G.g);

    // Sample the image along the middle
    if(flagImage)
      {
      itk::Index<3> idx;
      idx[0] = xAtoms[i].uIndex;
      idx[1] = xAtoms[i].vIndex;
      idx[2] = xImage->GetImageSize(2) >> 1;
      lImage->SetTuple1(i, xImage->GetInternalImage()->GetPixel(idx));
      }

    double du = xAtoms[i].u - 0.25 * floor(4 * xAtoms[i].u + 0.5);
    double dv = xAtoms[i].v - 0.25 * floor(4 * xAtoms[i].v + 0.5);
    double del = std::min(du * du, dv * dv);
    double q = exp(-0.01 * del * del);
    lDummy1->SetTuple1(i, q);
    }

  // Add all the quads
  for(MedialTriangleIterator itt(xGrid); !itt.IsAtEnd(); ++itt)
    {
    vtkIdType xTri[4];
    xTri[0] = itt.GetAtomIndex(0);
    xTri[1] = itt.GetAtomIndex(1);
    xTri[2] = itt.GetAtomIndex(2);
    pMedial->InsertNextCell(VTK_TRIANGLE, 3, xTri);
    }

  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();

  lCoordU->Delete();
  lCoordV->Delete();
  lRho->Delete();
  lRadius->Delete();
  lMetric->Delete();
  lNormals->Delete();
  lPoints->Delete();
  lDummy1->Delete();
  lImage->Delete();
  lBending->Delete();
  lRegularity->Delete();
  lContraOffDiag->Delete();
  lAngle->Delete();
  lMeanCurv->Delete();
  lGaussCurv->Delete();
  lAreaElement->Delete();
  pMedial->Delete();
  lSpoke2->Delete();
  lSpoke1->Delete();
}

vtkUnstructuredGrid *
ExportVolumeMeshToVTK(GenericMedialModel *xModel, size_t nSamples)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xModel->GetNumberOfInternalPoints(nSamples));
  
  // Allocate the polydata
  vtkUnstructuredGrid *pMedial = vtkUnstructuredGrid::New();
  pMedial->Allocate(xModel->GetNumberOfCells(nSamples));
  pMedial->SetPoints(lPoints);

  // Get the internals of the medial model
  MedialIterationContext *xGrid = xModel->GetIterationContext();
  MedialAtom *xAtoms = xModel->GetAtomArray();

  // Add all the points
  for(MedialInternalPointIterator it(xGrid,nSamples); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    SMLVec3d X = GetInternalPoint(it, xAtoms);
    lPoints->InsertPoint(i, X[0], X[1], X[2]);
    }

  // Add all the quads
  for(MedialInternalCellIterator itc(xGrid,nSamples); !itc.IsAtEnd(); ++itc)
    {
    vtkIdType xWedge[6];
    xWedge[0] = itc.GetInternalPointIndex(0, 0);
    xWedge[1] = itc.GetInternalPointIndex(1, 0);
    xWedge[2] = itc.GetInternalPointIndex(2, 0);
    xWedge[3] = itc.GetInternalPointIndex(0, 1);
    xWedge[4] = itc.GetInternalPointIndex(1, 1);
    xWedge[5] = itc.GetInternalPointIndex(2, 1);
    pMedial->InsertNextCell(VTK_WEDGE, 6, xWedge);
    }

  /*
  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();
  */

  return pMedial;
}


void ExportBoundaryMeshToVTK(
  GenericMedialModel *xModel,
  ITKImageWrapper<float> *xImage,
  const char *file)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xModel->GetNumberOfBoundaryPoints());
  
  // Allocate the array of normals
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(xModel->GetNumberOfBoundaryPoints());
  

  // Allocate the polydata
  vtkPolyData *pMedial = vtkPolyData::New();
  pMedial->Allocate(xModel->GetNumberOfBoundaryTriangles());
  pMedial->SetPoints(lPoints);
  
  // Add the arrays to the poly data
  pMedial->GetPointData()->SetNormals(lNormals);

  // Some more data arrays
  vtkFloatArray *lMeanCurv =  AddBndScalarField(pMedial, xModel, "Mean Curvature");
  vtkFloatArray *lGaussCurv = AddBndScalarField(pMedial, xModel, "Gauss Curvature");
  vtkFloatArray *lCurvPen = AddBndScalarField(pMedial, xModel, "Curvature Penalty");
  vtkFloatArray *lSqrMeanCrv = AddBndScalarField(pMedial, xModel, "SqrMeanCrv");
  vtkFloatArray *lCellJac = AddBndTriScalarField(pMedial, xModel, "Jacobian");

  vtkFloatArray *lDiffPenalty = 
    AddBndScalarField(pMedial, xModel, "Diffeomorphic Penalty");

  // Allocate the image intensity array
  bool flagImage = xImage && xImage->IsImageLoaded();
  vtkFloatArray *lImage = vtkFloatArray::New();
  if(flagImage)
    {
    lImage->SetName("Image");
    lImage->SetNumberOfComponents(1);
    lImage->SetNumberOfTuples(xModel->GetNumberOfBoundaryPoints());
    pMedial->GetPointData()->AddArray(lImage);
    }

  // Get the internals of the medial model
  MedialIterationContext *xGrid = xModel->GetIterationContext();
  MedialAtom *xAtoms = xModel->GetAtomArray();
  size_t nAtoms = xModel->GetNumberOfAtoms();

  // Compute certain derived terms
  SolutionData soldat(xGrid, xAtoms);
  soldat.ComputeIntegrationWeights();
  BoundaryCurvaturePenalty termBC(xModel);
  termBC.ComputeEnergy(&soldat);

  CodeTimer timer;

  // Add all the points
  for(MedialBoundaryPointIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    BoundaryAtom &B = GetBoundaryPoint(it, xAtoms);
    MedialAtom &A = xAtoms[it.GetAtomIndex()];
    lPoints->InsertNextPoint(B.X[0], B.X[1], B.X[2]);
    lNormals->SetTuple3(i, B.N[0], B.N[1], B.N[2]);
    lMeanCurv->SetTuple1(i, B.curv_mean);
    lGaussCurv->SetTuple1(i, B.curv_gauss);
    lCurvPen->SetTuple1(i, 4 * B.curv_mean * B.curv_mean - 2 * B.curv_gauss);

    SMLVec3d NX = (it.GetBoundarySide() > 0 ? 1.0:1.0) *  vnl_cross_3d(A.Xu, A.Xv);
    SMLVec3d NY = vnl_cross_3d(B.X_i[0], B.X_i[1]);

    SMLVec3d Kvec = termBC.GetCurvatureVector(i);
    double H2 = Kvec.squared_magnitude() / 
      (soldat.xBoundaryWeights[i] * soldat.xBoundaryWeights[i]);
    lSqrMeanCrv->SetTuple1(i, H2);

    // Compute the diffeomorphic penalty
    double diffpen = 0;
    timer.Start();
    for(size_t ia = 0; ia < nAtoms; ia++)
      {
      double d = (xAtoms[ia].X - B.X).squared_magnitude();
      if(d < xAtoms[ia].F)
        {
        diffpen += (xAtoms[ia].F - d) / xAtoms[ia].F;
        }
      }
    timer.Stop();
    lDiffPenalty->SetTuple1(i, sqrt(diffpen));

    // Sample the image along the middle
    if(flagImage)
      {
      itk::Index<3> idx;
      idx[0] = xAtoms[it.GetAtomIndex()].uIndex;
      idx[1] = xAtoms[it.GetAtomIndex()].vIndex;
      idx[2] = it.GetBoundarySide() ? 0 : xImage->GetImageSize(2) - 1;
      lImage->SetTuple1(i, xImage->GetInternalImage()->GetPixel(idx));
      }

    }

  cout << "Diff Pen Time: " << timer.Read() << endl;

  // Add all the quads
  for(MedialBoundaryTriangleIterator itt(xGrid); !itt.IsAtEnd(); ++itt)
    {
    // Insert the cells
    vtkIdType xTri[3];
    xTri[0] = itt.GetBoundaryIndex(0);
    xTri[1] = itt.GetBoundaryIndex(1);
    xTri[2] = itt.GetBoundaryIndex(2);
    pMedial->InsertNextCell(VTK_TRIANGLE, 3, xTri);

    // Set the Jacobian value
    MedialAtom &A0 = xAtoms[itt.GetAtomIndex(0)];
    MedialAtom &A1 = xAtoms[itt.GetAtomIndex(1)];
    MedialAtom &A2 = xAtoms[itt.GetAtomIndex(2)];

    // Compute the Xu and Xv vectors and the normal
    SMLVec3d XU = A1.X - A0.X,  XV = A2.X - A0.X;
    SMLVec3d NX = vnl_cross_3d(XU, XV);
    double gX2 = dot_product(NX, NX);

    size_t z = itt.GetBoundarySide();
    SMLVec3d YU = A1.xBnd[z].X - A0.xBnd[z].X;
    SMLVec3d YV = A2.xBnd[z].X - A0.xBnd[z].X;
    SMLVec3d NY = vnl_cross_3d(YU, YV);
    double J = dot_product(NY, NX) / gX2;
    lCellJac->SetTuple1(itt.GetIndex(), J);
    }

  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();

  lNormals->Delete();
  lPoints->Delete();
  pMedial->Delete();
  lImage->Delete();
  lMeanCurv->Delete();
  lGaussCurv->Delete();
  lCurvPen->Delete();
  lCellJac->Delete();
  lSqrMeanCrv->Delete();
}

/*
void ExportIntensityFieldToVTK(
  MedialAtomGrid *xGrid, 
  MedialAtom *xAtoms, 
  ITKImageWrapper<float> *imgField,
  
  const char *file)
{
  // Add the points to the poly data
  vtkPoints *lPoints = vtkPoints::New();
  lPoints->Allocate(xGrid->GetNumberOfBoundaryPoints());
  
  // Allocate the array of normals
  vtkFloatArray *lNormals = vtkFloatArray::New();
  lNormals->SetNumberOfComponents(3);
  lNormals->SetNumberOfTuples(xGrid->GetNumberOfBoundaryPoints());

  // Allocate the polydata
  vtkPolyData *pMedial = vtkPolyData::New();
  pMedial->Allocate(xGrid->GetNumberOfBoundaryQuads());
  pMedial->SetPoints(lPoints);
  pMedial->GetPointData()->SetNormals(lNormals);

  // Add all the points
  MedialBoundaryPointIterator *it = xGrid->NewBoundaryPointIterator();
  for(; !it->IsAtEnd(); ++(*it))
    {
    size_t i = it->GetIndex();
    BoundaryAtom &B = GetBoundaryPoint(it, xAtoms);
    lPoints->InsertNextPoint(B.X[0], B.X[1], B.X[2]);
    lNormals->SetTuple3(i, B.N[0], B.N[1], B.N[2]);
    }
  delete it;

  // Add all the quads
  MedialBoundaryQuadIterator *itq = xGrid->NewBoundaryQuadIterator();
  for(; !itq->IsAtEnd(); ++(*itq))
    {
    vtkIdType xQuad[4];
    xQuad[0] = itq->GetBoundaryIndex(0, 0);
    xQuad[1] = itq->GetBoundaryIndex(0, 1);
    xQuad[2] = itq->GetBoundaryIndex(1, 1);
    xQuad[3] = itq->GetBoundaryIndex(1, 0);
    pMedial->InsertNextCell(VTK_QUAD, 4, xQuad);
    }
  delete itq;

  vtkPolyDataWriter *fltWriter = vtkPolyDataWriter::New();
  fltWriter->SetInput(pMedial);
  fltWriter->SetFileName(file);
  fltWriter->SetFileTypeToBinary();
  fltWriter->Update();
  fltWriter->Delete();

  lNormals->Delete();
  lPoints->Delete();
  pMedial->Delete();
}
*/

/**
// An attempt to play with tetgen
#include "tetgen.h"

void TestTetGen(GenericMedialModel *xModel)
{
  // Get the internals of the medial model
  MedialIterationContext *xGrid = xModel->GetIterationContext();
  MedialAtom *xAtoms = xModel->GetAtomArray();

  // Create the tetgen structure
  tetgenio tgin, tgout;

  // Populate the node array
  size_t nv = xGrid->GetNumberOfBoundaryPoints();
  tgin.numberofpoints = nv + 8;
  tgin.pointlist = new REAL[tgin.numberofpoints * 3];

  // Also while we do this, compute the bounding box
  REAL bbmin[3], bbmax[3];

  // Add all the points as nodes
  for(MedialBoundaryPointIterator it(xGrid); !it.IsAtEnd(); ++it)
    {
    size_t i = it.GetIndex();
    BoundaryAtom &B = GetBoundaryPoint(it, xAtoms);
    for(size_t k = 0; k < 3; k++)
      {
      tgin.pointlist[3 * i + k] = B.X[k];
      bbmin[k] = (i == 0) ? B.X[k] : min(bbmin[k], B.X[k]);
      bbmax[k] = (i == 0) ? B.X[k] : max(bbmax[k], B.X[k]);
      }
    }

  // Expand the bounding box by 50%
  double w = max(bbmax[0]-bbmin[0], max(bbmax[1]-bbmin[1],bbmax[2]-bbmin[2]));
  for(size_t k = 0; k < 3; k++)
    {
    double c = 0.5 * (bbmin[k] + bbmax[k]);
    bbmin[k] = c - w;
    bbmax[k] = c + w;
    cout << "BB: " << bbmin[k] << " " << bbmax[k] << endl;
    }

  // Add the cube points
  size_t iv = nv * 3;
  tgin.pointlist[iv++] = bbmin[0]; tgin.pointlist[iv++] = bbmin[1]; tgin.pointlist[iv++] = bbmin[2];
  tgin.pointlist[iv++] = bbmax[0]; tgin.pointlist[iv++] = bbmin[1]; tgin.pointlist[iv++] = bbmin[2];
  tgin.pointlist[iv++] = bbmax[0]; tgin.pointlist[iv++] = bbmax[1]; tgin.pointlist[iv++] = bbmin[2];
  tgin.pointlist[iv++] = bbmin[0]; tgin.pointlist[iv++] = bbmax[1]; tgin.pointlist[iv++] = bbmin[2];
  tgin.pointlist[iv++] = bbmin[0]; tgin.pointlist[iv++] = bbmin[1]; tgin.pointlist[iv++] = bbmax[2];
  tgin.pointlist[iv++] = bbmax[0]; tgin.pointlist[iv++] = bbmin[1]; tgin.pointlist[iv++] = bbmax[2];
  tgin.pointlist[iv++] = bbmax[0]; tgin.pointlist[iv++] = bbmax[1]; tgin.pointlist[iv++] = bbmax[2];
  tgin.pointlist[iv++] = bbmin[0]; tgin.pointlist[iv++] = bbmax[1]; tgin.pointlist[iv++] = bbmax[2];

  // Populate the facet array
  size_t nf = xGrid->GetNumberOfBoundaryTriangles();
  tgin.numberoffacets = nf + 6;
  tgin.facetlist = new tetgenio::facet[tgin.numberoffacets];
  tgin.facetmarkerlist = new int[tgin.numberoffacets];

  // Generate the facets
  for(MedialBoundaryTriangleIterator itt(xGrid); !itt.IsAtEnd(); ++itt)
    {
    // Work with the current facet
    tgin.facetmarkerlist[itt.GetIndex()] = -1;
    tetgenio::facet &f = tgin.facetlist[itt.GetIndex()];

    // Initialize the poly
    f.numberofpolygons = 1;
    f.polygonlist = new tetgenio::polygon[1];
    f.numberofholes = 0;
    f.holelist = NULL;

    tetgenio::polygon &p = f.polygonlist[0];
    p.numberofvertices = 3;
    p.vertexlist = new int[3];
    p.vertexlist[0] = itt.GetBoundaryIndex(0);
    p.vertexlist[1] = itt.GetBoundaryIndex(1);
    p.vertexlist[2] = itt.GetBoundaryIndex(2);
    }

  // Generate a single hole (use one of the medial atoms)
  tgin.numberofholes = 1;
  tgin.holelist = new REAL[3];
  tgin.holelist[0] = xAtoms[0].X[0];
  tgin.holelist[1] = xAtoms[0].X[1];
  tgin.holelist[2] = xAtoms[0].X[2];

  // Generate the facets for the cube
  size_t cubeface[6][4] = {
      {0,1,2,3}, {2,1,5,6}, {1,0,4,5}, {0,3,7,4}, {3,2,6,7}, {7,6,5,4}};
  for(size_t j = 0; j < 6; j++)
    {
    tgin.facetmarkerlist[nf + j] = 0;
    tetgenio::facet &f = tgin.facetlist[nf + j];

    // Initialize the poly
    f.numberofpolygons = 1;
    f.polygonlist = new tetgenio::polygon[1];
    f.numberofholes = 0;
    f.holelist = NULL;

    tetgenio::polygon &p = f.polygonlist[0];
    p.numberofvertices = 4;
    p.vertexlist = new int[4];
    for(size_t k = 0; k < 4; k++)
      {
      p.vertexlist[k] = cubeface[j][k] + xGrid->GetNumberOfBoundaryPoints();
      }
    }

  // Write the objects we created
  tgin.save_nodes("tetin");
  tgin.save_poly("tetin");

  // Run the tetrahedralization
  tetrahedralize("pq1.412", &tgin, &tgout);

  tgout.save_nodes("tetout");
  tgout.save_elements("tetout");
  tgout.save_faces("tetout");
}
*/
