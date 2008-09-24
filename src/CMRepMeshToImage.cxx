#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkTriangleFilter.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_det.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientedRASImage.h>

using namespace itk;
using namespace std;

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif



typedef vnl_vector_fixed<vtkFloatingPointType, 3> Vec;
typedef itk::OrientedRASImage<float, 3> ImageType;
typedef itk::OrientedRASImage<unsigned char,3> ByteImage;

void MapPoint(const Vec &x, ImageType *img, double *vout)
{
  ImageType::PointType pt;
  pt[0] = x[0]; pt[1] = x[1]; pt[2] = x[2];
  ContinuousIndex<double,3> idx;
  img->TransformPhysicalPointToContinuousIndex(pt, idx);
  vout[0] = idx[0]; vout[1] = idx[1]; vout[2] = idx[2];
}

void VertSet(double **vtx, double **X, size_t t, size_t i0, size_t i1, size_t i2, size_t side)
{
  if(side == 1)
    {
    vtx[3*t  ] = X[i0]; vtx[3*t+1] = X[i1]; vtx[3*t+2] = X[i2];
    }
  else
    {
    vtx[3*t  ] = X[i2]; vtx[3*t+1] = X[i1]; vtx[3*t+2] = X[i0];
    }
}

void ScanTetrahedron(double *v1, double *v2, double *v3, double *v4, ByteImage *img)
{
  // Matrices for the point inside test
  typedef vnl_matrix_fixed<double, 4, 4> Mat;
  Mat D0, D1, D2, D3, D4;
  D0.fill(1);
  for(size_t col = 0; col < 3; col++)
    {
    D0(0,col) = v1[col]; 
    D0(1,col) = v2[col]; 
    D0(2,col) = v3[col]; 
    D0(3,col) = v4[col]; 
    }
  D1 = D0; D2 = D0; D3 = D0; D4 = D0;
  double det4 = vnl_det(D4);

  // Get the image region to iterate
  double xmin = D0.get_column(0).min_value(), xmax = D0.get_column(0).max_value();
  double ymin = D0.get_column(1).min_value(), ymax = D0.get_column(1).max_value();
  double zmin = D0.get_column(2).min_value(), zmax = D0.get_column(2).max_value();
  ByteImage::IndexType idx;
  idx[0] = (long) ceil(xmin); idx[1] = (long) ceil(ymin); idx[2] = (long) ceil(zmin);
  ByteImage::SizeType sz;
  sz[0] = (unsigned long) ( ceil(xmax) - ceil(xmin) );
  sz[1] = (unsigned long) ( ceil(ymax) - ceil(ymin) );
  sz[2] = (unsigned long) ( ceil(zmax) - ceil(zmin) );
  ByteImage::RegionType rgn(idx, sz);
  rgn.Crop(img->GetBufferedRegion());

  // Iterate over the region
  ImageRegionIteratorWithIndex<ByteImage> it(img, rgn);
  for( ; !it.IsAtEnd(); ++it)
    {
    ByteImage::IndexType idx = it.GetIndex();
    for(size_t col = 0; col < 3; col++)
      {
      D0(0,col) = idx[col];
      D1(1,col) = idx[col];
      D2(2,col) = idx[col];
      D3(3,col) = idx[col];
      }
    double det0 = vnl_det(D0);
    double det1 = vnl_det(D1);
    double det2 = vnl_det(D2);
    double det3 = vnl_det(D3);
    if(det0 > 0 && det1 > 0 && det2 > 0 && det3 > 0 && det4 > 0)
      it.Set(1);
    else if(det0 < 0 && det1 < 0 && det2 < 0 && det3 < 0 && det4 < 0)
      it.Set(1);
    }
}

int main(int argc, char *argv[])
{
  cout << "Usage: cmrep_fillmesh mesh.med.vtk reference.img output.img" << endl;
  cout << "  creates a binary image from a cm-rep mesh " << endl;
  if(argc < 4) return -1;

  // Load the mesh and the input image
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  vtkPolyData *mesh = reader->GetOutput();
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[2]);
  fltReader->Update();
  ImageType::Pointer img = fltReader->GetOutput();

  // Create a matching pair of images (increment and full)
  ByteImage::Pointer iOut = ByteImage::New();
  iOut->SetRegions(img->GetBufferedRegion());
  iOut->SetOrigin(img->GetOrigin());
  iOut->SetSpacing(img->GetSpacing());
  iOut->Allocate();
  iOut->FillBuffer(0);

  ByteImage::Pointer iTmp = ByteImage::New();
  iTmp->SetRegions(iOut->GetBufferedRegion());
  iTmp->Allocate();
  iTmp->FillBuffer(0);

  // Get the resolution
  int res[3];
  res[0] = iOut->GetBufferedRegion().GetSize()[0];
  res[1] = iOut->GetBufferedRegion().GetSize()[1];
  res[2] = iOut->GetBufferedRegion().GetSize()[2];

  // Create an array of vertices for each 'wedge'
  double **vtx = new double *[42];
  double **X = new double*[10];
  for(size_t q = 0; q < 10; q++)
    X[q] = new double[3];

  // Get the spoke arrays, radius
  vtkDataArray *spoke[2], *rad;
  spoke[0] = mesh->GetPointData()->GetArray("Spoke1");
  spoke[1] = mesh->GetPointData()->GetArray("Spoke2");
  rad = mesh->GetPointData()->GetArray("Radius Function");

  if(spoke[0] == NULL || spoke[1] == NULL || rad == NULL)
    {
    cerr << "Spoke/radius arrays missing in mesh" << endl;
    return -1;
    }

  vtkPolyData *dump = vtkPolyData::New();
  vtkPoints *dumppts = vtkPoints::New();
  vtkIntArray *dumparr = vtkIntArray::New();
  dumparr->SetName("wedge");
  dump->Allocate();
  dump->SetPoints(dumppts);
  dump->GetPointData()->AddArray(dumparr);

  // Iterate over the cells in the image
  size_t nt = mesh->GetNumberOfPolys();
  vtkCellArray *tri = mesh->GetPolys();
  vtkIdType npts; vtkIdType *pts;
  typedef vnl_vector_fixed<vtkFloatingPointType, 3> Vec;
  size_t itri = 0;
  size_t iwedge = 0;
  for(tri->InitTraversal(); tri->GetNextCell(npts,pts); itri++)
    {
    for(size_t side = 0; side < 2; side++)
      {
      // Get the vertices
      Vec x[3], b[3], f[3];
      for(size_t i=0;i<3;i++)
        {
        x[i] = Vec(mesh->GetPoint(pts[i]));
        Vec spk(spoke[side]->GetTuple(pts[i]));
        b[i] = x[i] + spk;
        }

      // Add the vertices to each of the wedge's quad faces
      Vec q[3];
      q[0] = 0.25 * (x[1] + b[1] + x[2] + b[2]);
      q[1] = 0.25 * (x[2] + b[2] + x[0] + b[0]);
      q[2] = 0.25 * (x[0] + b[0] + x[1] + b[1]);
      Vec C = (q[0] + q[1] + q[2]) / 3.0;


      // Convert to image coordinates
      MapPoint(x[0], img, X[0]); MapPoint(x[1], img, X[1]); MapPoint(x[2], img, X[2]);
      MapPoint(b[0], img, X[3]); MapPoint(b[1], img, X[4]); MapPoint(b[2], img, X[5]);
      MapPoint(q[0], img, X[6]); MapPoint(q[1], img, X[7]); MapPoint(q[2], img, X[8]);
      MapPoint(C, img, X[9]);

      ScanTetrahedron(X[0], X[1], X[2], X[9], iOut);
      ScanTetrahedron(X[3], X[5], X[4], X[9], iOut);

      ScanTetrahedron(X[0], X[8], X[1], X[9], iOut);
      ScanTetrahedron(X[1], X[8], X[4], X[9], iOut);
      ScanTetrahedron(X[4], X[8], X[3], X[9], iOut);
      ScanTetrahedron(X[3], X[8], X[0], X[9], iOut);

      ScanTetrahedron(X[1], X[6], X[2], X[9], iOut);
      ScanTetrahedron(X[2], X[6], X[5], X[9], iOut);
      ScanTetrahedron(X[5], X[6], X[4], X[9], iOut);
      ScanTetrahedron(X[4], X[6], X[1], X[9], iOut);

      ScanTetrahedron(X[2], X[7], X[0], X[9], iOut);
      ScanTetrahedron(X[0], X[7], X[3], X[9], iOut);
      ScanTetrahedron(X[3], X[7], X[5], X[9], iOut);
      ScanTetrahedron(X[5], X[7], X[2], X[9], iOut);

      // Create an array of faces
      VertSet(vtx, X, 0, 0, 1, 2, side);
      VertSet(vtx, X, 1, 3, 5, 4, side);

      VertSet(vtx, X, 2, 0, 8, 1, side);
      VertSet(vtx, X, 3, 1, 8, 4, side);
      VertSet(vtx, X, 4, 4, 8, 3, side);
      VertSet(vtx, X, 5, 3, 8, 0, side);

      VertSet(vtx, X, 6, 1, 6, 2, side);
      VertSet(vtx, X, 7, 2, 6, 5, side);
      VertSet(vtx, X, 8, 5, 6, 4, side);
      VertSet(vtx, X, 9, 4, 6, 1, side);

      VertSet(vtx, X,10, 2, 7, 0, side);
      VertSet(vtx, X,11, 0, 7, 3, side);
      VertSet(vtx, X,12, 3, 7, 5, side);
      VertSet(vtx, X,13, 5, 7, 2, side);

      for(size_t r = 0; r < 14; r++)
        {
        vtkIdType ids[3];
        ids[0] = dump->GetPoints()->InsertNextPoint(vtx[r*3+0]);
        ids[1] = dump->GetPoints()->InsertNextPoint(vtx[r*3+1]);
        ids[2] = dump->GetPoints()->InsertNextPoint(vtx[r*3+2]);
        dump->InsertNextCell(VTK_TRIANGLE, 3, ids);
        dumparr->InsertNextTuple1(iwedge);
        dumparr->InsertNextTuple1(iwedge);
        dumparr->InsertNextTuple1(iwedge);
        }
      iwedge++;

      // for(size_t r= 0; r<42;r++)
      //  cout << vtx[r][0] << " " << vtx[r][1] << " "  << vtx[r][2] << endl;

      // drawBinaryTrianglesSheetFilled(iOut->GetBufferPointer(), res, vtx, 14);
      // drawBinaryTrianglesFilled(iOut->GetBufferPointer(), res, vtx, 14);
      // unsigned char *ptOut = iOut->GetBufferPointer();
      // unsigned char *ptTmp = iTmp->GetBufferPointer();
      // for(size_t j = 0; j < iTmp->GetBufferedRegion().GetNumberOfPixels(); j++)
      //   {
      //  if(ptTmp[j] > 0)
      //    {
      //    ptOut[j] = 1; ptTmp[j] = 0;
      //    }
      //  }
      }
      cout << "." << flush;
      if((itri + 1) % 80 == 0)
        cout << " " << itri+1 << "/" << nt << endl;
    }

  cout << endl;

  typedef itk::ImageFileWriter<ByteImage> WriterType;
  WriterType::Pointer fltWrite= WriterType::New();
  fltWrite->SetFileName(argv[3]);
  fltWrite->SetInput(iOut);
  fltWrite->Update();
}
