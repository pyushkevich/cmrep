#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkTriangleFilter.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_det.h>
#include <vnl/vnl_inverse.h>
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

void MapPoint(const Vec &x, ImageType *img, Vec &vout)
{
  ImageType::PointType pt;
  pt[0] = x[0]; pt[1] = x[1]; pt[2] = x[2];
  ContinuousIndex<double,3> idx;
  img->TransformRASPhysicalPointToContinuousIndex(pt, idx);
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

void ScanTetrahedron(Vec *X, Vec *V, size_t i1, size_t i2, size_t i3, size_t i4, ImageType::Pointer warp[3])
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
  if(rgn.Crop(warp[0]->GetBufferedRegion()))
    {
    // Set up interpolators for the deformation fields
    vnl_matrix_fixed<double, 4, 3> VM, Q;
    for(size_t col = 0; col < 3; col++)
      {
      VM[0][col] = V[i1][col];
      VM[1][col] = V[i2][col];
      VM[2][col] = V[i3][col];
      VM[3][col] = V[i4][col];
      }
    Q = vnl_inverse(D4) * VM;

    // Iterate over the region
    ImageRegionIteratorWithIndex<ImageType> itx(warp[0], rgn);
    ImageRegionIteratorWithIndex<ImageType> ity(warp[1], rgn);
    ImageRegionIteratorWithIndex<ImageType> itz(warp[2], rgn);
    for( ; !itx.IsAtEnd(); ++itx, ++ity, ++itz)
      {
      ByteImage::IndexType idx = itx.GetIndex();
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
        itx.Set(v[0]);
        ity.Set(v[1]);
        itz.Set(v[2]);
        }
      }
    }
}

int usage()
{
  cout << "cmrep_genwarp: " << endl;
  cout << "  create a warp field between two images based on cm-rep normalization" << endl;
  cout << "usage:" << endl;
  cout << "  cmrep_to_warp [opts] cmrep1 cmrep2 ref_image warp_name" << endl;
  cout << "parameters:" << endl;
  cout << "  cmrep1, cmrep2  : VTK meshes (cmr_fit mesh dir, defX.med.vtk) " << endl;
  cout << "  ref_image       : Image based on which warp field is evaluated" << endl;
  cout << "options:" << endl;
  cout << "  -e TEXT         : extension for warp files (def: nii.gz)" << endl;
  cout << "  -x Xmax Xstep   : xi coord. sampling (def: 1 0.1)" << endl;
  cout << "                    specifies whether the warp field extends past " << endl;
  cout << "                    the model's boundary (i.e., 2 0.1). The step " << endl;
  cout << "                    size affects the size of tetrahedra used to " << endl;
  cout << "                    the warp field" << endl;
  cout << "  -f              : flip atom normal, i.e., map spoke1 in cmrep1 to " << endl;
  cout << "                    spoke2 in cmrep2" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 5) return usage();

  // Get the main parameters
  string fnwarp = argv[argc-1];
  string fnimage = argv[argc-2];
  string fncmrep2 = argv[argc-3];
  string fncmrep1 = argv[argc-4];

  // Get the options
  string fnext = "nii.gz";
  double xiMax = 1, xiStep = 0.1;
  bool flip = false;

  // Parse the options
  for(size_t iopt = 1; iopt < argc-4; iopt++)
    {
    if(!strcmp(argv[iopt], "-e"))
      {
      fnext = argv[++iopt];
      }
    else if(!strcmp(argv[iopt], "-f"))
      {
      flip = true;
      }
    else if(!strcmp(argv[iopt], "-x"))
      {
      xiMax = atof(argv[++iopt]);
      xiStep = atof(argv[++iopt]);
      }
    else
      {
      cerr << "Unknown option " << argv[iopt] << endl;
      return usage();
      }
    }

  // Load the two cm-rep meshes
  vtkPolyDataReader *r1 = vtkPolyDataReader::New();
  r1->SetFileName(fncmrep1.c_str());
  r1->Update();
  vtkPolyData *m1 = r1->GetOutput();

  vtkPolyDataReader *r2 = vtkPolyDataReader::New();
  r2->SetFileName(fncmrep2.c_str());
  r2->Update();
  vtkPolyData *m2 = r2->GetOutput();

  // Read the reference image
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(fnimage.c_str());
  fltReader->Update();
  ImageType::Pointer img = fltReader->GetOutput();

  // Create warp image (for now only the forward warp)
  ImageType::Pointer warp[3], winverse[3];
  for(size_t d = 0; d < 3; d++)
    {
    warp[d] = ImageType::New();
    warp[d]->SetRegions(img->GetBufferedRegion());
    warp[d]->SetOrigin(img->GetOrigin());
    warp[d]->SetSpacing(img->GetSpacing());
    warp[d]->Allocate();
    warp[d]->FillBuffer(0.0);
    }

  // Get the resolution
  int res[3];
  res[0] = warp[0]->GetBufferedRegion().GetSize()[0];
  res[1] = warp[0]->GetBufferedRegion().GetSize()[1];
  res[2] = warp[0]->GetBufferedRegion().GetSize()[2];

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

        // Convert cm-rep1 to image coordinates
        Vec X[10], V[10];
        MapPoint(a1[0], img, X[0]); MapPoint(a1[1], img, X[1]); MapPoint(a1[2], img, X[2]);
        MapPoint(b1[0], img, X[3]); MapPoint(b1[1], img, X[4]); MapPoint(b1[2], img, X[5]);
        MapPoint(q1[0], img, X[6]); MapPoint(q1[1], img, X[7]); MapPoint(q1[2], img, X[8]);
        MapPoint(C1, img, X[9]);

        // For each of these ten vertices, compute the displacement
        V[0] = a2[0] - a1[0]; V[1] = a2[1] - a1[1]; V[2] = a2[2] - a1[2];
        V[3] = b2[0] - b1[0]; V[4] = b2[1] - b1[1]; V[5] = b2[2] - b1[2];
        V[6] = q2[0] - q1[0]; V[7] = q2[1] - q1[1]; V[8] = q2[2] - q1[2];
        V[9] = C2-C1;

        // Interpolate the warp for each tetrahedron 
        ScanTetrahedron(X, V, 0, 1, 2, 9, warp);
        ScanTetrahedron(X, V, 3, 5, 4, 9, warp);

        ScanTetrahedron(X, V, 0, 8, 1, 9, warp);
        ScanTetrahedron(X, V, 1, 8, 4, 9, warp);
        ScanTetrahedron(X, V, 4, 8, 3, 9, warp);
        ScanTetrahedron(X, V, 3, 8, 0, 9, warp);

        ScanTetrahedron(X, V, 1, 6, 2, 9, warp);
        ScanTetrahedron(X, V, 2, 6, 5, 9, warp);
        ScanTetrahedron(X, V, 5, 6, 4, 9, warp);
        ScanTetrahedron(X, V, 4, 6, 1, 9, warp);

        ScanTetrahedron(X, V, 2, 7, 0, 9, warp);
        ScanTetrahedron(X, V, 0, 7, 3, 9, warp);
        ScanTetrahedron(X, V, 3, 7, 5, 9, warp);
        ScanTetrahedron(X, V, 5, 7, 2, 9, warp);
        }
      }
      cout << "." << flush;
      if((itri + 1) % 80 == 0)
        cout << " " << itri+1 << "/" << nt << endl;
    }

  cout << endl;

  // Write the warp fields out
  string suffix = "xvec." + fnext;
  for(size_t d = 0; d < 3; d++)
    {
    suffix[0] = d + 'x';
    string fname = fnwarp + suffix;

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer fltWrite= WriterType::New();
    fltWrite->SetFileName(fname.c_str());
    fltWrite->SetInput(warp[d]);
    fltWrite->Update();
    }
}
