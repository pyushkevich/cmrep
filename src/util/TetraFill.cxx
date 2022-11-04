#include <iostream>
#include <string>
#include <sstream>
#include <set>
#include <limits>

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_det.h>
#include "ReadWriteVTK.h"

#include <itkVectorImage.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

using namespace std;

typedef vnl_matrix_fixed<double, 4, 4> MatrixType;
typedef itk::VectorImage<double, 3> ImageType;


int usage()
{
  cout << "tetfill - Fills the interior of a tetrahedral mesh" << endl;
  cout << "usage: " << endl;
  cout << "   tetfill [options] tetinput.vtk refimage.nii outimage.nii" << endl;
  cout << "options: " << endl;
  cout << "  -c Name           Use cell array 'Name' to fill the tetrahedra" << endl;
  cout << "  -s source.vtk     Read the cell array data from source.vtk" << endl;
  cout << "  -b Value          Background value (default: 0)" << endl;
  cout << "  -N                Background value is nan" << endl;
  return -1;
}

void ScanTetrahedron(
  double *v1, double *v2, double *v3, double *v4, 
  ImageType *img, unsigned int ncomp, double *fillvalue)
{
  // Create a fill pixel
  ImageType::PixelType fillpixel(fillvalue, ncomp);

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
  ImageType::IndexType idx;
  idx[0] = (long) ceil(xmin); idx[1] = (long) ceil(ymin); idx[2] = (long) ceil(zmin);
  ImageType::SizeType sz;
  sz[0] = (unsigned long) ( ceil(xmax) - ceil(xmin) );
  sz[1] = (unsigned long) ( ceil(ymax) - ceil(ymin) );
  sz[2] = (unsigned long) ( ceil(zmax) - ceil(zmin) );
  ImageType::RegionType rgn(idx, sz);
  rgn.Crop(img->GetBufferedRegion());

  // Iterate over the region
  itk::ImageRegionIteratorWithIndex<ImageType> it(img, rgn);
  for( ; !it.IsAtEnd(); ++it)
    {
    ImageType::IndexType idx = it.GetIndex();
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
      it.Set(fillpixel);
    else if(det0 < 0 && det1 < 0 && det2 < 0 && det3 < 0 && det4 < 0)
      it.Set(fillpixel);
    }
}

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 4) return usage();

  // Parse the optional parameters
  int ch;
  string inCellArray = "";
  string sourceCellArray = "";
  double inBackground = 0.0;
  while ((ch = getopt(argc, argv, "c:b:s:N")) != -1)
    {
    switch(ch)
      {
    case 'b': inBackground = atof(optarg); break;
    case 'c': inCellArray = optarg; break;
    case 's': sourceCellArray = optarg; break;
    case 'N': inBackground = std::numeric_limits<double>::quiet_NaN(); break;
    default: return usage();
      }
    }

  if(optind + 3 != argc) return usage();
  const char * fnTetra = argv[optind++];
  const char * fnRef = argv[optind++];
  const char * fnOutput = argv[optind++];

  // Read the reference image
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(fnRef);
  reader->Update();
  ImageType::Pointer ref = reader->GetOutput();

  // Allocate the output image
  ImageType::Pointer out = ImageType::New();
  out->SetRegions(ref->GetBufferedRegion());
  out->SetOrigin(ref->GetOrigin());
  out->SetSpacing(ref->GetSpacing());
  out->SetDirection(ref->GetDirection());
  
  // Read the tetrahedral mesh
  vtkUnstructuredGridReader *ugreader = vtkUnstructuredGridReader::New();
  ugreader->SetFileName(fnTetra);
  ugreader->Update();
  vtkUnstructuredGrid *tet = ugreader->GetOutput();

  // Read the data array
  vtkDataArray *da = NULL;

  if(inCellArray.length()) 
    {
    if(sourceCellArray.length())
      {
      vtkUnstructuredGridReader *sourceugreader = vtkUnstructuredGridReader::New();
      sourceugreader->SetFileName(sourceCellArray.c_str());
      sourceugreader->Update();
      vtkUnstructuredGrid *sourcetet = sourceugreader->GetOutput();
      da = sourcetet->GetCellData()->GetArray(inCellArray.c_str());
      if (da->GetNumberOfTuples()!=tet->GetNumberOfCells())
        {
        cout << "No. of cells in source mesh should equal that in input mesh" << endl;
        return usage();
        }
      }
    else
      {
      cout << "\"" << inCellArray << "\"" << endl;
      da = tet->GetCellData()->GetArray(inCellArray.c_str());
      if (da->GetNumberOfTuples()!=tet->GetNumberOfCells())
        {
        cout << "No. of cells in input mesh is not equal to isize of cell data array" << endl;
        return usage();
        }

      }
    }

  // Output dimensionality
  unsigned int ncomp = da ? da->GetNumberOfComponents() : 1;
  out->SetNumberOfComponentsPerPixel(ncomp);
  out->Allocate();

  // Background pixel
  ImageType::PixelType bkg_pixel(ncomp);
  bkg_pixel.Fill(inBackground);
  out->FillBuffer(bkg_pixel);
  std::cout << "Background value: " << inBackground << std::endl;

  // Scan convert each of the tetrahedra
  for(size_t ic = 0; ic < (size_t) tet->GetNumberOfCells(); ic++)
    {
    // Get the cell ID
    vtkCell *cell = tet->GetCell(ic);

    // Get the point coordinates
    double p[4][3];
    tet->GetPoint(cell->GetPointId(0), p[0]);
    tet->GetPoint(cell->GetPointId(1), p[1]);
    tet->GetPoint(cell->GetPointId(2), p[2]);
    tet->GetPoint(cell->GetPointId(3), p[3]);

    // Remap each point to image index
    double q[4][3];
    for(size_t i=0; i<4; i++)
      {
      itk::Point<double,3> pt; 
      itk::ContinuousIndex<double, 3> ci;
      pt[0] = p[i][0]; pt[1] = p[i][1]; pt[2] = p[i][2];

      // Map RAS physical point to continuous index
      pt[0] = -pt[0]; pt[1] = -pt[1];
      ref->TransformPhysicalPointToContinuousIndex(pt, ci);
      q[i][0] = ci[0]; q[i][1] = ci[1]; q[i][2] = ci[2];
      }

    // Create the fill value
    double default_fill_value = 1.0;
    double *fillvec = da ? da->GetTuple(ic) : &default_fill_value;

    ScanTetrahedron(q[0], q[1], q[2], q[3], out, ncomp, fillvec);
    }

  // Write the output image
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(fnOutput);
  writer->SetInput(out);
  writer->Update();
}




    




