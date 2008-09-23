#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkTriangleFilter.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_det.h>
#include <vnl/vnl_inverse.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientedRASImage.h>
#include <itkInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>

using namespace itk;
using namespace std;

#ifndef vtkFloatingPointType
#define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif



typedef vnl_vector_fixed<vtkFloatingPointType, 3> Vec;
typedef itk::Image<float, 3> ImageType;
typedef itk::Image<unsigned char,3> ByteImage;

typedef itk::InterpolateImageFunction<ImageType, double> InterpType;

void MapPoint(const Vec &x, ImageType *img, Vec &vout)
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

void ScanTetrahedron(
  Vec *X, Vec *Y, size_t i1, size_t i2, size_t i3, size_t i4, 
  ImageType::Pointer itrg, InterpType::Pointer isrc)
{
  // Matrices for the point inside test
  typedef vnl_matrix_fixed<double, 4, 4> Mat;
  Mat D0, D1, D2, D3, D4;
  D0.fill(1);
  for(size_t col = 0; col < 3; col++)
    {
    D0(0,col) = X[i1][col]; 
    D0(1,col) = X[i2][col]; 
    D0(2,col) = X[i3][col]; 
    D0(3,col) = X[i4][col]; 
    }
  D1 = D0; D2 = D0; D3 = D0; D4 = D0;
  double det4 = vnl_det(D4);

  // Get the image region to iterate
  double xmin = D0.get_column(0).min_value(), xmax = D0.get_column(0).max_value();
  double ymin = D0.get_column(1).min_value(), ymax = D0.get_column(1).max_value();
  double zmin = D0.get_column(2).min_value(), zmax = D0.get_column(2).max_value();
  ImageType::IndexType idx;
  idx[0] = (long) ceil(xmin); idx[1] = (long) ceil(ymin); idx[2] = (long) ceil(zmin);
  ImageType::SizeType sz;
  sz[0] = (unsigned long) ( ceil(xmax) - ceil(xmin) );
  sz[1] = (unsigned long) ( ceil(ymax) - ceil(ymin) );
  sz[2] = (unsigned long) ( ceil(zmax) - ceil(zmin) );
  ImageType::RegionType rgn(idx, sz);

  // Proceed only if region fits inside
  if(rgn.Crop(itrg->GetBufferedRegion()))
    {
    // Set up interpolators for the deformation fields
    vnl_matrix_fixed<double, 4, 3> VM, Q;
    for(size_t col = 0; col < 3; col++)
      {
      VM[0][col] = Y[i1][col];
      VM[1][col] = Y[i2][col];
      VM[2][col] = Y[i3][col];
      VM[3][col] = Y[i4][col];
      }
    Q = vnl_inverse(D4) * VM;

    // Iterate over the region    
    for(ImageRegionIteratorWithIndex<ImageType> it(itrg, rgn); 
      !it.IsAtEnd(); ++it)
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
      if(
        (det0 > 0 && det1 > 0 && det2 > 0 && det3 > 0 && det4 > 0) ||
        (det0 < 0 && det1 < 0 && det2 < 0 && det3 < 0 && det4 < 0))
        {
        // Interpolate the position of the voxel inside the tetrahedron
        Vec v = D0.get_row(0) * Q;

        // Sample the source image at this location
        InterpType::ContinuousIndexType cindex;
        cindex[0] = v[0]; cindex[1] = v[1]; cindex[2] = v[2];

        if(isrc->IsInsideBuffer(cindex))
          it.Set(isrc->EvaluateAtContinuousIndex(cindex));
        }
      }
    }
}

int usage()
{
  cout << "cmrep_warpimage: " << endl;
  cout << "  Use cm-rep correspondence to warp source image to target image" << endl;
  cout << "usage:" << endl;
  cout << "  cmrep_warpimage [opts] cmrep_trg cmrep_src image_trg image_src image_out" << endl;
  cout << "parameters:" << endl;
  cout << "  cmrep_trg  : VTK meshes (cmr_fit mesh dir, defX.med.vtk) of 'target' cm-rep" << endl;
  cout << "  cmrep_src  : VTK meshes (cmr_fit mesh dir, defX.med.vtk) of 'source' cm-rep" << endl;
  cout << "  image_trg             : Image that defines voxel space for the output" << endl;
  cout << "  image_src             : Image that will be sampled (resliced)" << endl;
  cout << "  image_out             : Output image" << endl;
  cout << "options:" << endl;
  cout << "  -x Xmax Xstep   : xi coord. sampling (def: 1 0.1)" << endl;
  cout << "                    specifies whether the warp field extends past " << endl;
  cout << "                    the model's boundary (i.e., 2 0.1). The step " << endl;
  cout << "                    size affects the size of tetrahedra used to " << endl;
  cout << "                    the warp field" << endl;
  cout << "  -i model        : Interpolation model (0 = NN, 1 = trilinear, 3 = cubic)" << endl;
  cout << "                    Default is 1 (linear)" << endl;
  cout << "  -f              : flip atom normal, i.e., map spoke1 in cmrep1 to " << endl;
  cout << "                    spoke2 in cmrep2" << endl;
  cout << "  -b value        : Background value that is used to fill in voxels that " << endl;
  cout << "                    are outside of the target cm-rep" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 6) return usage();

  // Get the main parameters
  string fn_image_out = argv[argc-1];
  string fn_image_src = argv[argc-2];
  string fn_image_trg = argv[argc-3];
  string fn_cmrep_src = argv[argc-4];
  string fn_cmrep_trg = argv[argc-5];

  // Get the options
  double xiMax = 1, xiStep = 0.1;
  float bkg = 0.0f;
  int interp = 1;
  bool flip = false;

  // Parse the options
  for(size_t iopt = 1; iopt < argc-5; iopt++)
    {
    if(!strcmp(argv[iopt], "-f"))
      {
      flip = true;
      }
    else if(!strcmp(argv[iopt], "-x"))
      {
      xiMax = atof(argv[++iopt]);
      xiStep = atof(argv[++iopt]);
      }
    else if(!strcmp(argv[iopt], "-i"))
      {
      interp = atoi(argv[++iopt]);
      }
    else if(!strcmp(argv[iopt], "-b"))
      {
      bkg = atof(argv[++iopt]);
      }
    else
      {
      cerr << "Unknown option " << argv[iopt] << endl;
      return usage();
      }
    }

  // Load the two cm-rep meshes
  vtkPolyDataReader *r1 = vtkPolyDataReader::New();
  r1->SetFileName(fn_cmrep_trg.c_str());
  r1->Update();
  vtkPolyData *m1 = r1->GetOutput();

  vtkPolyDataReader *r2 = vtkPolyDataReader::New();
  r2->SetFileName(fn_cmrep_src.c_str());
  r2->Update();
  vtkPolyData *m2 = r2->GetOutput();

  // Read the reference image
  typedef itk::Image<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(fn_image_trg.c_str());
  fltReader->Update();
  ImageType::Pointer img = fltReader->GetOutput();

  // Read the source image
  ReaderType::Pointer fltReaderSrc = ReaderType::New();
  fltReaderSrc->SetFileName(fn_image_src.c_str());
  fltReaderSrc->Update();
  ImageType::Pointer isrc = fltReaderSrc->GetOutput();

  // Create output image
  ImageType::Pointer ires = ImageType::New();
  ires->SetRegions(img->GetBufferedRegion());
  ires->SetOrigin(img->GetOrigin());
  ires->SetSpacing(img->GetSpacing());
  ires->Allocate();
  ires->FillBuffer(bkg);

  // Create an interpolator for the output image
  InterpType::Pointer ifun;
  if(interp == 0)
    {
    ifun = itk::NearestNeighborInterpolateImageFunction<ImageType, double>::New();
    }
  else if(interp == 1)
    {
    ifun = itk::LinearInterpolateImageFunction<ImageType, double>::New();
    }
  else
    {
    itk::BSplineInterpolateImageFunction<ImageType>::Pointer
      ispline =itk::BSplineInterpolateImageFunction<ImageType>::New();
    ispline->SetSplineOrder(interp);
    ifun = ispline.GetPointer();
    }
  ifun->SetInputImage(isrc);

  // Get the resolution
  int res[3];
  res[0] = ires->GetBufferedRegion().GetSize()[0];
  res[1] = ires->GetBufferedRegion().GetSize()[1];
  res[2] = ires->GetBufferedRegion().GetSize()[2];

  // Get the spoke arrays, radius
  vtkDataArray *spoke1[2], *spoke2[2], *rad1, *rad2;
  spoke1[0] = m1->GetPointData()->GetArray("Spoke1");
  spoke1[1] = m1->GetPointData()->GetArray("Spoke2");
  rad1 = m1->GetPointData()->GetArray("Radius Function");
  
  if(!flip)
    {
    spoke2[0] = m2->GetPointData()->GetArray("Spoke1");
    spoke2[1] = m2->GetPointData()->GetArray("Spoke2");
    }
  else
    {
    spoke2[0] = m2->GetPointData()->GetArray("Spoke2");
    spoke2[1] = m2->GetPointData()->GetArray("Spoke1");
    }

  rad2 = m2->GetPointData()->GetArray("Radius Function");

  if(spoke1[0] == NULL || spoke1[1] == NULL || rad1 == NULL)
    {
    cerr << "Spoke/radius arrays missing in cmrep1" << endl;
    return -1;
    }
  if(spoke2[0] == NULL || spoke2[1] == NULL || rad2 == NULL)
    {
    cerr << "Spoke/radius arrays missing in cmrep2" << endl;
    return -1;
    }
  if(m1->GetNumberOfPoints() != m2->GetNumberOfPoints())
    {
    cerr << "CM-rep meshes have different number of points " << endl;
    return -1;
    }
  if(m1->GetNumberOfPolys() != m2->GetNumberOfPolys())
    {
    cerr << "CM-rep meshes have different number of triangles " << endl;
    return -1;
    }

  // Iterate over the cells in the image
  size_t nt = m1->GetNumberOfPolys();
  vtkCellArray *tri = m1->GetPolys();
  vtkIdType npts; vtkIdType *pts;
  typedef vnl_vector_fixed<vtkFloatingPointType, 3> Vec;
  size_t itri = 0;
  size_t iwedge = 0;
  for(tri->InitTraversal(); tri->GetNextCell(npts,pts); itri++)
    {
    for(size_t side = 0; side < 2; side++)
      {
      // Get the medial vertices and the spokes
      Vec x1[3], s1[3];
      Vec x2[3], s2[3];
      for(size_t i=0;i<3;i++)
        {
        x1[i] = Vec(m1->GetPoint(pts[i]));
        x2[i] = Vec(m2->GetPoint(pts[i]));
        s1[i] = Vec(spoke1[side]->GetTuple(pts[i]));
        s2[i] = Vec(spoke2[side]->GetTuple(pts[i]));
        }

      // Loop over the wedges
      for(double xi = 0; xi < xiMax; xi += xiStep)
        {
        
        // Get the extents of the wedge in each cm-rep
        Vec a1[3], b1[3], a2[3], b2[3];
        for(size_t i=0;i<3;i++)
          {
          a1[i] = x1[i] + xi * s1[i];
          a2[i] = x2[i] + xi * s2[i];
          b1[i] = x1[i] + (xi + xiStep) * s1[i];
          b2[i] = x2[i] + (xi + xiStep) * s2[i];
          }

        // Break quad faces of wedges by inserting vertices
        Vec q1[3], q2[3];
        q1[0] = 0.25 * (a1[1] + b1[1] + a1[2] + b1[2]);
        q1[1] = 0.25 * (a1[2] + b1[2] + a1[0] + b1[0]);
        q1[2] = 0.25 * (a1[0] + b1[0] + a1[1] + b1[1]);
        q2[0] = 0.25 * (a2[1] + b2[1] + a2[2] + b2[2]);
        q2[1] = 0.25 * (a2[2] + b2[2] + a2[0] + b2[0]);
        q2[2] = 0.25 * (a2[0] + b2[0] + a2[1] + b2[1]);

        // Compute centroid of each wedge
        Vec C1 = (q1[0] + q1[1] + q1[2]) / 3.0;
        Vec C2 = (q2[0] + q2[1] + q2[2]) / 3.0;

        // Convert reference cm-rep to image coordinates
        Vec X[10];
        MapPoint(a1[0], img, X[0]); MapPoint(a1[1], img, X[1]); MapPoint(a1[2], img, X[2]);
        MapPoint(b1[0], img, X[3]); MapPoint(b1[1], img, X[4]); MapPoint(b1[2], img, X[5]);
        MapPoint(q1[0], img, X[6]); MapPoint(q1[1], img, X[7]); MapPoint(q1[2], img, X[8]);
        MapPoint(C1, img, X[9]);

        // Convert source cm-rep to image coordinates
        Vec Y[10];
        MapPoint(a2[0], isrc, Y[0]); MapPoint(a2[1], isrc, Y[1]); MapPoint(a2[2], isrc, Y[2]);
        MapPoint(b2[0], isrc, Y[3]); MapPoint(b2[1], isrc, Y[4]); MapPoint(b2[2], isrc, Y[5]);
        MapPoint(q2[0], isrc, Y[6]); MapPoint(q2[1], isrc, Y[7]); MapPoint(q2[2], isrc, Y[8]);
        MapPoint(C2, isrc, Y[9]);

        // Interpolate the warp for each tetrahedron 
        ScanTetrahedron(X, Y, 0, 1, 2, 9, ires, ifun);
        ScanTetrahedron(X, Y, 3, 5, 4, 9, ires, ifun);

        ScanTetrahedron(X, Y, 0, 8, 1, 9, ires, ifun);
        ScanTetrahedron(X, Y, 1, 8, 4, 9, ires, ifun);
        ScanTetrahedron(X, Y, 4, 8, 3, 9, ires, ifun);
        ScanTetrahedron(X, Y, 3, 8, 0, 9, ires, ifun);

        ScanTetrahedron(X, Y, 1, 6, 2, 9, ires, ifun);
        ScanTetrahedron(X, Y, 2, 6, 5, 9, ires, ifun);
        ScanTetrahedron(X, Y, 5, 6, 4, 9, ires, ifun);
        ScanTetrahedron(X, Y, 4, 6, 1, 9, ires, ifun);

        ScanTetrahedron(X, Y, 2, 7, 0, 9, ires, ifun);
        ScanTetrahedron(X, Y, 0, 7, 3, 9, ires, ifun);
        ScanTetrahedron(X, Y, 3, 7, 5, 9, ires, ifun);
        ScanTetrahedron(X, Y, 5, 7, 2, 9, ires, ifun);
        }
      }
      cout << "." << flush;
      if((itri + 1) % 80 == 0)
        cout << " " << itri+1 << "/" << nt << endl;
    }

  cout << endl;

  // Write the output image out
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fltWrite = WriterType::New();
  fltWrite->SetFileName(fn_image_out.c_str());
  fltWrite->SetInput(ires);
  fltWrite->Update();
}
