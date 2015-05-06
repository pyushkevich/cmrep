#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkCellLocator.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkLODActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkBYUReader.h>
#include <vtkSTLReader.h>
#include <vtkTriangleFilter.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>

#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkCommand.h>
#include <itkImage.h>

#include "DrawTriangles.h"
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_cross.h>

using namespace std;
using namespace itk;

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

int usage()
{
  cout << "meshdiff - compare two VTK meshes" << endl;
  cout << "usage: " << endl;
  cout << "   meshdiff [options] mesh1.vtk mesh2.vtk" << endl;
  cout << "options: " << endl;
  cout << "   -s X.XX             Size of the voxel for overlap/distance measurements" << endl;
  cout << "   -d file.vtk         Save pointwise mesh-to-mesh distance to a file" << endl;
  cout << "   -d-no-attr          Don't copy attributes to output of -d" << endl;
  return -1;
}

typedef itk::Image<unsigned char, 3> ByteImageType;
typedef itk::Image<float, 3> FloatImageType;

vtkPolyData *ReadVTKData(string fn)
{
  vtkPolyData *p1 = NULL;

  // Choose the reader based on extension

  if(fn.rfind(".byu") == fn.length() - 4)
    {
    vtkBYUReader *reader = vtkBYUReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else if(fn.rfind(".vtk") == fn.length() - 4)
    {
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(fn.c_str());
    reader->Update();
    p1 = reader->GetOutput();
    }
  else
    {
    cout << "Could not find a reader for " << fn << endl;
    cout << fn.rfind(".byu") << endl;
    cout << fn.rfind(".vtk") << endl;
    return NULL;
    }

  // Convert the model to triangles
  vtkTriangleFilter *tri = vtkTriangleFilter::New();
  tri->SetInputData(p1);
  cout << "Converting to triangles ..." << endl;
  tri->Update();
  vtkPolyData *pd = tri->GetOutput();

  return pd;
}

void ProgressCommand(Object *, const EventObject &, void *)
{ 
  cout << "." << flush; 
}

void ComputeDistanceTransform(ByteImageType *binary, FloatImageType::Pointer &dist)
{
  CStyleCommand::Pointer xCommand = CStyleCommand::New();
  xCommand->SetCallback(ProgressCommand);
  
  // Compute the forward distance transform
  typedef DanielssonDistanceMapImageFilter<
    ByteImageType, FloatImageType> DistanceMapType;
  DistanceMapType::Pointer fltDistance = DistanceMapType::New();
  fltDistance->SetInput(binary);
  fltDistance->SetSquaredDistance(false);
  fltDistance->SetInputIsBinary(true);
  fltDistance->SetUseImageSpacing(true);
  fltDistance->AddObserver(ProgressEvent(),xCommand);
  fltDistance->Update();
  FloatImageType::Pointer d1 = fltDistance->GetOutput();

  // Create an inverted image
  ByteImageType::Pointer imgInvert = ByteImageType::New();
  imgInvert->SetRegions(binary->GetBufferedRegion());
  imgInvert->SetSpacing(binary->GetSpacing());
  imgInvert->SetOrigin(binary->GetOrigin());
  imgInvert->Allocate();
  ImageRegionIterator<ByteImageType> it1(imgInvert, binary->GetBufferedRegion());
  ImageRegionConstIterator<ByteImageType> it2(binary, binary->GetBufferedRegion());
  for( ; !it1.IsAtEnd(); ++it1, ++it2)
    it1.Set( it2.Value() > 0 ? 0 : 255 );

  // Compute the inverse distance transform
  DistanceMapType::Pointer fltDistanceInv = DistanceMapType::New();
  fltDistanceInv->SetInput(imgInvert);
  fltDistanceInv->SetSquaredDistance(false);
  fltDistanceInv->SetInputIsBinary(true);
  fltDistanceInv->SetUseImageSpacing(true);
  fltDistanceInv->AddObserver(ProgressEvent(),xCommand);
  fltDistanceInv->Update();
  FloatImageType::Pointer d2 = fltDistanceInv->GetOutput();
  
  // Compute the Sum of the two images
  typedef SubtractImageFilter<
    FloatImageType,FloatImageType,FloatImageType> AdderType;
  AdderType::Pointer fltAdder = AdderType::New();
  fltAdder->SetInput1(d1);
  fltAdder->SetInput2(d2);
  fltAdder->Update();

  dist = fltAdder->GetOutput();

  /*
  typedef ImageFileWriter<FloatImageType> Junk;
  Junk::Pointer junk = Junk::New();
  junk->SetInput(dist);
  junk->SetFileName("junk.img.gz");
  junk->Update();
  */
}

void IntegrateDistance(vtkPolyData *pd, FloatImageType *dist, 
  double &xOneNorm, double &xTwoNorm, double &xInfNorm)
{
  // Set up an interpolator
  typedef itk::LinearInterpolateImageFunction<FloatImageType, double> FuncType;
  FuncType::Pointer fnInterp = FuncType::New();
  fnInterp->SetInputImage(dist);

  // Clear the norms
  xOneNorm = xTwoNorm = xInfNorm = 0.0;
  
  // Go over all the points in the image
  for(unsigned int i = 0; i < pd->GetNumberOfPoints(); i++)
    {
    // Get the point
    FuncType::PointType P;
    P[0] = pd->GetPoint(i)[0];
    P[1] = pd->GetPoint(i)[1];
    P[2] = pd->GetPoint(i)[2];
    double d = fabs(fnInterp->Evaluate(P));
    
    if(d > xInfNorm) xInfNorm = d;
    xOneNorm += d;
    xTwoNorm += d * d;
    }

  xOneNorm /= pd->GetNumberOfPoints();
  xTwoNorm = sqrt( xTwoNorm / pd->GetNumberOfPoints() );
}

void ScanConvertPolyData(vtkPolyData *pd, double *b0, double *b1, int *res, double xVox, ByteImageType *img)
{
  // Create a vertex table from the polydata
  unsigned int nt = pd->GetNumberOfPolys(), it = 0;
  double **vtx = new double*[nt*3];  

  vtkCellArray *poly = pd->GetPolys();
  vtkIdType npts;
  vtkIdType *pts;
  for(poly->InitTraversal();poly->GetNextCell(npts,pts);)
    {
    for(unsigned int i=0;i<3;i++)
      {
      vtkFloatingPointType *x = pd->GetPoints()->GetPoint(pts[i]);
      vtx[it] = (double *) malloc(3*sizeof(double));
      for(unsigned int j=0;j<3;j++)
        {
        vtx[it][j] = (x[j] - b0[j]) / xVox;
        }

      ++it;
      }      
    }

  // Create a ITK image to store the results
  ByteImageType::RegionType region;
  region.SetSize(0,res[0]); region.SetSize(1,res[1]); region.SetSize(2,res[2]);
  img->SetRegions(region);
  img->Allocate();
  img->FillBuffer(0);

  // Convert the polydata to an image
  drawBinaryTrianglesFilled(img->GetBufferPointer(), res, vtx, nt);

  // Set the origin and spacing of the image
  ByteImageType::SpacingType xSpacing;
  ByteImageType::PointType xOrigin;
  for(unsigned int d = 0; d < 3; d++)
    {
    xOrigin[d] = b0[d];
    xSpacing[d] = xVox;
    }
  img->SetOrigin(xOrigin);
  img->SetSpacing(xSpacing);
 
  /*
  static int igg = 1;
  typedef ImageFileWriter<ByteImageType> Junk;
  Junk::Pointer junk = Junk::New();
  junk->SetInput(img);
  if (igg==1)
    {
    junk->SetFileName("chunk1.img.gz");
    igg = 0;
    }
  else
    junk->SetFileName("chunk2.img.gz");
  junk->Update();
  */
}

inline vtkFloatingPointType TriangleArea(
  const vnl_vector_fixed<vtkFloatingPointType,3> &A, 
  const vnl_vector_fixed<vtkFloatingPointType,3> &B, 
  const vnl_vector_fixed<vtkFloatingPointType,3> &C)
{
  return 0.5 * vnl_cross_3d(B - A, C - A).magnitude();
}


void ComputeAreaElement(vtkPolyData *poly, vnl_vector<vtkFloatingPointType> &elt)
{
  // For each triangle in the polydata compute its area
  vtkIdType nCells = poly->GetNumberOfCells();
  
  // Initialize the area elt array
  elt.set_size(poly->GetNumberOfPoints());
  elt.fill(0.0);
  
  for(vtkIdType iCell = 0; iCell < nCells; iCell++)
    {
    // Get the points in this cell
    vtkIdType nPoints, *xPoints;
    poly->GetCellPoints(iCell, nPoints, xPoints);
    
    // Only triangles are admitted
    if(nPoints != 3)
      { 
      cerr << "Irregular face (n = " << nPoints << ") detected!" << endl;
      break;
      }

    // Get the three points
    vnl_vector_fixed<vtkFloatingPointType, 3> X0(poly->GetPoint(xPoints[0]));
    vnl_vector_fixed<vtkFloatingPointType, 3> X1(poly->GetPoint(xPoints[1]));
    vnl_vector_fixed<vtkFloatingPointType, 3> X2(poly->GetPoint(xPoints[2]));

    // Compute the area
    double xArea = TriangleArea(X0, X1, X2);
    if(xArea < 0)
      {
      cerr << "Negative area returned at cell " << iCell << endl;
      break;
      }

    // Add the area to all points
    elt[xPoints[0]] += xArea; elt[xPoints[1]] += xArea; elt[xPoints[2]] += xArea;
    }
}

void ComputeExactMeshToMeshDistance(vtkPolyData *source, vtkPolyData *target, vnl_vector<vtkFloatingPointType> &dist)
{
  // Create a cell locator from the target mesh
  vtkCellLocator *xLocator = vtkCellLocator::New();
  xLocator->SetDataSet(target);
  xLocator->BuildLocator();

  // Expand the distance array to match the number of points in the source
  dist.set_size(source->GetNumberOfPoints());

  // Compute each distance
  for(size_t i = 0; i < dist.size(); i++)
    {
    vtkFloatingPointType xClosest[3];
    vtkIdType cellId;
    int subId;
    xLocator->FindClosestPoint(source->GetPoint(i), xClosest, cellId, subId, dist[i]);
    }

  // Take square root
  dist = dist.apply(sqrt);
}

int main(int argc, char **argv)
{
  // Parameter values
  double bb[2][3] = {{0.,0.,0.},{128.,128.,128.}};
  double xVox = 1.0;
  int res[3];
  bool flagNoAttr = false;

  // Input specifications
  string fn1, fn2, fnMesh = "";

  // Check the parameters
  if(argc < 3) return usage();

  // Get the file names
  fn1 = argv[argc-2];
  fn2 = argv[argc-1];

  // Parse the optional parameters
  try 
    {
    for(int i=1;i<argc-2;i++)
      {
      string arg = argv[i];
      if(arg == "-s")
        {
        xVox = atof(argv[++i]);
        }
      else if(arg == "-d")
        {
        fnMesh = argv[++i];
        }
      else if(arg == "-d-no-attr")
        {
        flagNoAttr = true;
        }
      }
    }
  catch(...) 
    {
    cout << "Error parsing command line options!" << endl;
    return usage();
    }

  // The real program begins here

  // Read the appropriate meshes
  vtkPolyData *p1 = ReadVTKData(fn1);
  vtkPolyData *p2 = ReadVTKData(fn2);

  // Get the extents of the data
  vtkFloatingPointType *b1 = p1->GetBounds();  
  vtkFloatingPointType *b2 = p2->GetBounds();  

  // Compute the bounding box
  for(size_t i = 0; i < 3; i++) 
    {
    bb[0][i] = (b1[2*i] < b2[2*i] ? b1[2*i] : b2[2*i]) - 4;
    bb[1][i] = (b1[2*i+1] < b2[2*i+1] ? b1[2*i+1] : b2[2*i+1]) + 4;
    cout << "Bounds[" << i << "]: " << bb[0][i] << " to " << bb[1][i] << endl;
    res[i] = (int)(ceil((bb[1][i] - bb[0][i]) / xVox));
    }

  // Scan convert the mesh to images
  ByteImageType::Pointer i1 = ByteImageType::New();
  ByteImageType::Pointer i2 = ByteImageType::New();
  ScanConvertPolyData(p1, bb[0], bb[1], res, xVox, i1);
  ScanConvertPolyData(p2, bb[0], bb[1], res, xVox, i2);

  // Compute the volume overlap between images
  unsigned long nUnion = 0, nIntersect = 0, nReference = 0, nModel = 0;
  typedef ImageRegionConstIterator<ByteImageType> IterType;
  IterType it1(i1, i1->GetBufferedRegion());
  IterType it2(i2, i2->GetBufferedRegion());
  for( ; !it1.IsAtEnd() ; ++it1, ++it2)
    {
    if(it1.Value() > 0 && it2.Value() > 0)
      nIntersect++;
    if(it1.Value() > 0 || it2.Value() > 0)
      nUnion++;
    if(it2.Value() > 0)
      nReference++;
    if(it1.Value() > 0)
      nModel++;
    }

  /*
  // Compute the distance transform for each image
  FloatImageType::Pointer d1 = FloatImageType::New();
  FloatImageType::Pointer d2 = FloatImageType::New();
  ComputeDistanceTransform(i1, d1);
  ComputeDistanceTransform(i2, d2);

  // Integrate the distance transform over each mesh
  double xOneNorm[2], xTwoNorm[2], xInfNorm[2];
  IntegrateDistance(p1, d2, xOneNorm[0], xTwoNorm[0], xInfNorm[0]);
  IntegrateDistance(p2, d1, xOneNorm[1], xTwoNorm[1], xInfNorm[1]);
  */

  // Compute exact pointwise distances from mesh to mesh
  vnl_vector<vtkFloatingPointType> d1, d2, ae1, ae2;
  ComputeExactMeshToMeshDistance(p1, p2, d1);
  ComputeExactMeshToMeshDistance(p2, p1, d2);
  ComputeAreaElement(p1, ae1);
  ComputeAreaElement(p2, ae2);

  // Save the pointwise mesh distance if asked to
  if(fnMesh != "")
    {
    if(flagNoAttr)
      {
      while(p1->GetPointData()->GetNumberOfArrays() > 0)
        p1->GetPointData()->RemoveArray(p1->GetPointData()->GetArrayName(0));
      }
    vtkFloatArray *dist = vtkFloatArray::New();
    dist->SetName("Distance");
    dist->SetNumberOfComponents(1);
    dist->SetNumberOfTuples(d1.size());
    p1->GetPointData()->AddArray(dist);
    for(size_t q=0; q < d1.size(); q++)
      dist->SetTuple1(q, d1[q]);
    vtkPolyDataWriter *pdw = vtkPolyDataWriter::New();
    pdw->SetFileName(fnMesh.c_str());
    pdw->SetInput(p1);
    pdw->Update();
    }

  // Scale the area by the total area to get a weight
  ae1 /= ae1.one_norm(); ae2 /= ae2.one_norm();

  // Compute the maximum and average distances (this should be a surface integral...)
  double xOneNorm[2], xTwoNorm[2], xInfNorm[2];
  xOneNorm[0] = dot_product(d1, ae1); 
  xTwoNorm[0] = sqrt( dot_product(d1.apply(vnl_math_sqr), ae1) );
  xInfNorm[0] = d1.inf_norm();
  
  xOneNorm[1] = dot_product(d2, ae2); 
  xTwoNorm[1] = sqrt( dot_product(d2.apply(vnl_math_sqr), ae2) );
  xInfNorm[1] = d2.inf_norm();
  
  // Report our findings
  cout << endl;
  cout << "Unbiased overlap: " << nIntersect * 1.0 / nUnion << endl;
  cout << "Biased overlap: " << nIntersect * 1.0 / nReference << endl;
  cout << "Dice similarity index: " << 2.0 * nIntersect / (nReference + nModel) << endl;
  cout << "Reverse biased overlap between meshes is " << nIntersect * 1.0 / nModel << endl;
  cout << "Mesh 1 Volume is " << nModel * xVox * xVox * xVox << endl;
  cout << "Mesh 2 Volume is " << nReference * xVox * xVox * xVox << endl;
  cout << "Distance from mesh 1 to mesh 2 is " << xOneNorm[0] << "; " << xTwoNorm[0] << "; " << xInfNorm[0] << endl;
  cout << "Distance from mesh 2 to mesh 1 is " << xOneNorm[1] << "; " << xTwoNorm[1] << "; " << xInfNorm[1] << endl;
  cout << "RESULT: " 
    << nIntersect * 1.0 / nUnion << "\t"
    << nIntersect * 1.0 / nReference << "\t"
    << 2.0 * nIntersect / (nReference + nModel) << "\t"
    << xOneNorm[0] << "\t" << xOneNorm[1] << "\t"
    << xTwoNorm[0] << "\t" << xTwoNorm[1] << "\t"
    << xInfNorm[0] << "\t" << xInfNorm[1] << endl;
}
