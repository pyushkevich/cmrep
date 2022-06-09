#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkTriangleFilter.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_det.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientedRASImage.h>

using namespace itk;
using namespace std;


typedef vnl_vector_fixed<double, 3> Vec;
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

void VertSet(double **vtx, Vec *X, size_t t, size_t i0, size_t i1, size_t i2, size_t side)
{
  if(side == 1)
    {
    vtx[3*t  ] = X[i0].data_block(); 
    vtx[3*t+1] = X[i1].data_block(); 
    vtx[3*t+2] = X[i2].data_block();
    }
  else
    {
    vtx[3*t  ] = X[i2].data_block(); 
    vtx[3*t+1] = X[i1].data_block(); 
    vtx[3*t+2] = X[i0].data_block();
    }
}


void ScanTetrahedron(
  Vec *X, double *val, 
  size_t i1, size_t i2, size_t i3, size_t i4, 
  ImageType *image)
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

  itk::Index<3> idx;
  idx[0] = (long) ceil(xmin); idx[1] = (long) ceil(ymin); idx[2] = (long) ceil(zmin);

  itk::Size<3> sz;
  sz[0] = (unsigned long) ( ceil(xmax) - ceil(xmin) );
  sz[1] = (unsigned long) ( ceil(ymax) - ceil(ymin) );
  sz[2] = (unsigned long) ( ceil(zmax) - ceil(zmin) );

  itk::ImageRegion<3> rgn(idx, sz);

  // Proceed only if region fits inside
  if(rgn.Crop(image->GetBufferedRegion()))
    {
    // Set up interpolators for the deformation fields
    vnl_vector_fixed<double, 4> VM, Q;
    VM[0] = val[i1]; VM[1] = val[i2]; VM[2] = val[i3]; VM[3] = val[i4];
    Q = vnl_inverse(D4) * VM;

    // Iterate over the region
    ImageRegionIteratorWithIndex<ImageType> it(image, rgn);
    for( ; !it.IsAtEnd(); ++it)
      {
      itk::Index<3> idx = it.GetIndex();
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
        double v = dot_product(D0.get_row(0), Q);
        it.Set(v);
        }
      }
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
  cout << "Usage: cmrep_fillmesh [options]  mesh.med.vtk reference.img output.img" << endl;
  cout << "  creates a binary image from a cm-rep mesh " << endl;
  cout << "Options: " << endl;
  cout << "  -a ArrayName      : interpolates an array from the VTK mesh instead" << endl;
  cout << "  -d                : interpolates the depth coordinate (xi)" << endl;
  cout << "  -x xiMax xiStep   : specifies whether the warp field extends past " << endl;
  cout << "                      the model's boundary (i.e., 2 0.1)" << endl;

  if(argc < 4) return -1;

  // Is an array specified?
  std::string arrayName;
  bool extractDepth = false;
  double xiMax = 1.0, xiStep = 0.1;
  for(int i = 1; i < argc-3; i++)
    {
    if(strcmp(argv[i], "-a") == 0)
      arrayName = argv[++i];

    else if(strcmp(argv[i], "-d") == 0)
      extractDepth = true;

    else if(strcmp(argv[i], "-x") == 0)
      {
      xiMax = atof(argv[++i]);
      xiStep = atof(argv[++i]);
      }
    }

  // Load the mesh and the input image
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(argv[argc-3]);
  reader->Update();

  vtkPolyData *mesh = reader->GetOutput();
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[argc-2]);
  fltReader->Update();
  ImageType::Pointer img = fltReader->GetOutput();

  // Create a matching pair of images (increment and full)
  ImageType::Pointer iOut = ImageType::New();
  iOut->SetRegions(img->GetBufferedRegion());
  iOut->SetOrigin(img->GetOrigin());
  iOut->SetDirection(img->GetDirection());
  iOut->SetSpacing(img->GetSpacing());
  iOut->Allocate();
  iOut->FillBuffer(0.0);

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
  Vec X[10];

  // Get the spoke arrays, radius
  bool haveSpokes = false, haveCurvedSpokes = false;

  vtkDataArray *spoke[2], *rad, *cspoke[2];

  spoke[0] = mesh->GetPointData()->GetArray("Spoke1");
  spoke[1] = mesh->GetPointData()->GetArray("Spoke2");
  rad = mesh->GetPointData()->GetArray("Radius Function");

  if(spoke[0] != NULL && spoke[1] != NULL && rad != NULL)
    {
    haveSpokes = true;
    }
  else
    {
    // The new version of this code accepts curved spokes, which are a list of 
    // points (streamline) for each vertex in the mesh
    cspoke[0] = mesh->GetPointData()->GetArray("CurvedSpoke1");
    cspoke[1] = mesh->GetPointData()->GetArray("CurvedSpoke2");

    if(cspoke[0] != NULL && cspoke[1] != NULL)
      haveCurvedSpokes = true;
    }

  if(!haveSpokes && !haveCurvedSpokes)
    {
    cerr << "Spoke/radius arrays missing in mesh" << endl;
    return -1;
    }

  // If asked for an array, retrieve it
  vtkDataArray *userArray = NULL;
  if(arrayName.length())
    userArray = mesh->GetPointData()->GetArray(arrayName.c_str());

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
  vtkIdType npts; const vtkIdType *pts;
  typedef vnl_vector_fixed<double, 3> Vec;
  size_t itri = 0;
  size_t iwedge = 0;

  // Generate an array of sample points for each vertex
  typedef vnl_vector_fixed<double, 3> Vec3;
  typedef std::vector<Vec> StreamLine;
  typedef std::vector<StreamLine> StreamField;

  // Interpolate the spokes into a list of streamlines - to unify with curved spokes
  StreamField streams[2];
  for(int side = 0; side < 2; side++)
    {
    for(int v = 0; v < mesh->GetNumberOfPoints(); v++)
      {
      // Streamline for this vertex for this side
      StreamLine sl;

      // Get the coordinate of the vertex itself
      Vec x0; x0[0] = mesh->GetPoint(v)[0]; x0[1] = mesh->GetPoint(v)[1]; x0[2] = mesh->GetPoint(v)[2];
      sl.push_back(x0);

      // If spokes, interpolate spokes
      if(haveSpokes)
        {
        Vec spk(spoke[side]->GetTuple(v));

        // Sample along the spoke
        for(double xi = 0; xi < xiMax; xi += xiStep)
          {
          Vec px = x0 + (xi + xiStep) * spk;
          sl.push_back(px);
          }
        }
      else if(haveCurvedSpokes)
        {
        double *spoke_coord = cspoke[side]->GetTuple(v);
        for(int i = 0; i < cspoke[side]->GetNumberOfComponents(); i+=3)
          {
          Vec px;
          px[0] = spoke_coord[i]; px[1] = spoke_coord[i+1]; px[2] = spoke_coord[i+2];
          sl.push_back(px);
          }
        }

        streams[side].push_back(sl);
      }
    }

  for(tri->InitTraversal(); tri->GetNextCell(npts,pts); itri++)
    {
    for(size_t side = 0; side < 2; side++)
      {
      // Store pointers to three streamlines
      StreamLine &sl0 = streams[side][pts[0]];
      StreamLine &sl1 = streams[side][pts[1]];
      StreamLine &sl2 = streams[side][pts[2]];
      
      // Loop over all the streamlines
      for(int i = 0; i < sl0.size() - 1; i++)
        {
        Vec x[3], b[3];
        x[0] = sl0[i];   x[1] = sl1[i];   x[2] = sl2[i];
        b[0] = sl0[i+1]; b[1] = sl1[i+1]; b[2] = sl2[i+1];

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

        // Create an array of values to interpolate
        double val[10];
        if(userArray)
          {
          val[0] = userArray->GetTuple1(pts[0]);
          val[1] = userArray->GetTuple1(pts[1]);
          val[2] = userArray->GetTuple1(pts[2]);
          val[3] = userArray->GetTuple1(pts[0]);
          val[4] = userArray->GetTuple1(pts[1]);
          val[5] = userArray->GetTuple1(pts[2]);
          }
        else if(extractDepth)
          {
          double xi = xiStep * i;
          val[0] = val[1] = val[2] = side ? -xi : xi;
          val[3] = val[4] = val[5] = side ? -(xi + xiStep) : (xi + xiStep);
          }
        else
          {
          val[0] = val[1] = val[2] = 1.0;
          val[3] = val[4] = val[5] = 1.0;
          }

        // Interpolate the other values
        val[6] = 0.25 * (val[1] + val[4] + val[2] + val[5]);
        val[7] = 0.25 * (val[2] + val[5] + val[0] + val[3]);
        val[8] = 0.25 * (val[0] + val[3] + val[1] + val[4]);
        val[9] = (val[6] + val[7] + val[8]) / 3.0;

        ScanTetrahedron(X, val, 0, 1, 2, 9, iOut);
        ScanTetrahedron(X, val, 3, 5, 4, 9, iOut);

        ScanTetrahedron(X, val, 0, 8, 1, 9, iOut);
        ScanTetrahedron(X, val, 1, 8, 4, 9, iOut);
        ScanTetrahedron(X, val, 4, 8, 3, 9, iOut);
        ScanTetrahedron(X, val, 3, 8, 0, 9, iOut);

        ScanTetrahedron(X, val, 1, 6, 2, 9, iOut);
        ScanTetrahedron(X, val, 2, 6, 5, 9, iOut);
        ScanTetrahedron(X, val, 5, 6, 4, 9, iOut);
        ScanTetrahedron(X, val, 4, 6, 1, 9, iOut);

        ScanTetrahedron(X, val, 2, 7, 0, 9, iOut);
        ScanTetrahedron(X, val, 0, 7, 3, 9, iOut);
        ScanTetrahedron(X, val, 3, 7, 5, 9, iOut);
        ScanTetrahedron(X, val, 5, 7, 2, 9, iOut);

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
        }

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

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer fltWrite= WriterType::New();
  fltWrite->SetFileName(argv[argc-1]);
  fltWrite->SetInput(iOut);
  fltWrite->Update();
}
