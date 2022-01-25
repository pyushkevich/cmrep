#include <iostream>
#include <string>
#include <sstream>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include "ReadWriteVTK.h"
#include <vtkTriangleFilter.h>

#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>

#include <itk_to_nifti_xform.h>
#include "DrawTriangles.h"

using namespace std;
using namespace itk;

typedef Image<float, 3> RefImageType;

int usage()
{
  cout << "mesh2img - fill in a 3D mesh, producing a binary image" << endl;
  cout << "description:" << endl;
  cout << "   Given an 3D mesh, this command fill scan-convert the faces in " << endl;
  cout << "   the mesh to a 3D image, optionally filling the interior (-f option)"  << endl;
  cout << "   The user must specify the corners of the image in the coordinate"  << endl;
  cout << "   space of the mesh using the -o and -s options, as well as the resolution"  << endl;
  cout << "   of the image using the -r option. " << endl;
  cout << "   To find out the spatial extent of the mesh, call the program without" << endl;
  cout << "   the -o and -s parameters. The extent will be printed out." << endl;
  cout << "usage: " << endl;
  cout << "   mesh2img [options] output.nii.gz" << endl;
  cout << "options: " << endl;
  cout << "   -i <mesh>                   Input mesh file (VTK, STL, etc)" << endl;
  cout << "   -f                          Fill the interior of the mesh also" << endl;
  cout << "   -o NN NN NN                 Image origin in x,y,z in mesh space" << endl;
  cout << "   -s NN NN NN                 Image size in x,y,z in mesh space" << endl;
  cout << "   -r NN NN NN                 Image resolution in x,y,z in pixels" << endl;
  cout << "   -a NN NN NN XX              Automatic bounding box computation. First three " << endl;
  cout << "                               parameters are the voxel size, the last is the margin in mm" << endl;
  cout << "   -ref file                   reference file for the origin, size, resolution, etc." << endl;
  cout << "   -c number                   set segmentation mask value (default: 1) " << endl;
  return -1;
}

int main(int argc, char **argv)
{
  // Parameter values
  double bb[2][3] = {{0.,0.,0.},{128.,128.,128.}};
  double xAutoSpacing[3] = {0,0,0};
  double xAutoMargin = 0;
  int res[3] = {128,128,128};
  bool doFloodFill = false, flagAutoSize = false;

  bool useRef = false;
  int segColor = 1;  // default value for segmentation mask

  // Input specifications
  enum MeshType {NONE, BYU, STL, VTK} iMeshType = NONE;
  string fnOutput, fnInput;

  // for Reference image by wook
  string fnRef;
  RefImageType::Pointer ref;

  // Check the parameters
  if(argc < 2) return usage();

  // Get the file names
  fnOutput = argv[argc-1];

  // Parse the optional parameters
  try 
    {
    for(int i=1;i<argc-1;i++)
      {
      string arg = argv[i];
      if(arg == "-i")
        {
        fnInput = argv[++i];
        }
      else if(arg == "-vtk")
        {
        fnInput = argv[++i];
        iMeshType = VTK;
        }
      else if(arg == "-stl")
        {
        fnInput = argv[++i];
        iMeshType = STL;
        }
      else if(arg == "-byu")
        {
        fnInput = argv[++i];
        iMeshType = BYU;
        }
      else if(arg == "-ref")
        {
        fnRef = argv[++i];
        useRef = true;
        flagAutoSize = true;
        }
      else if(arg == "-o")
        {
        bb[0][0] = atof(argv[++i]);
        bb[0][1] = atof(argv[++i]);
        bb[0][2] = atof(argv[++i]);
        }
      else if(arg == "-s")
        {
        bb[1][0] = atof(argv[++i]);
        bb[1][1] = atof(argv[++i]);
        bb[1][2] = atof(argv[++i]);
        }
      else if(arg == "-r")
        {
        res[0] = atoi(argv[++i]);
        res[1] = atoi(argv[++i]);
        res[2] = atoi(argv[++i]);
        }      
      else if(arg == "-a")
        {
        xAutoSpacing[0] = atof(argv[++i]);
        xAutoSpacing[1] = atof(argv[++i]);
        xAutoSpacing[2] = atof(argv[++i]);
        xAutoMargin = atof(argv[++i]);
        flagAutoSize = true;
        }      
      else if(arg == "-f")
        {
        doFloodFill = true;
        }
      else if(arg == "-c")
        {
        segColor = atoi(argv[++i]);
        }
      else
        {
        cout << "Unknown argument " << arg << endl;
        return usage();
        }
      }
    }
  catch(...) 
    {
    cout << "Error parsing command line options!" << endl;
    return usage();
    }

  // If no file given, return
  if(fnInput.size() == 0)
    {
    cerr << "No input mesh specified!" << endl;
    return usage();
    }    	

  //cout << endl << "useRef ="  << useRef << endl;
  if(useRef)
    {
    cout << "Reference file: "<< fnRef << endl;
    // Read the reference image
    typedef itk::ImageFileReader<RefImageType> RefReaderType;
    RefReaderType::Pointer Refreader = RefReaderType::New();
    Refreader->SetFileName(fnRef.c_str());
    cout << "Reading the Reference file ..." << endl;
    Refreader->Update();
    ref = Refreader->GetOutput();
    }

  // Print the parameters
  cout << endl << "Parameters:" << endl;
  if(!flagAutoSize)
    {
    cout << "   Image resolution  : " << res[0] << " by " << res[1] << " by " << res[2] << endl;
    cout << "   Image box origin  : {" << bb[0][0] << ", " << bb[0][1] << ", " << bb[0][2] << "}" << endl;
    cout << "   Image box size    : {" << bb[1][0] << ", " << bb[1][1] << ", " << bb[1][2] << "}" << endl << endl;
    }
  else if( (flagAutoSize) && (useRef)) 
    {
    //      cout << " else if (flagAutoSize && useRef) " << endl;
    bb[0][0] = ref->GetOrigin()[0];
    bb[0][1] = ref->GetOrigin()[1];
    bb[0][2] = ref->GetOrigin()[2];
    cout << "Origin: ("<<bb[0][0]<<", "<<bb[0][1]<<", "<<bb[0][2]<<")"<<endl;

    bb[1][0] = ref->GetRequestedRegion().GetSize()[0];
    bb[1][1] = ref->GetRequestedRegion().GetSize()[1];
    bb[1][2] = ref->GetRequestedRegion().GetSize()[2];
    cout << "Size: ("<<bb[1][0]<<", "<<bb[1][1]<<", "<<bb[1][2]<<")"<<endl;

    xAutoSpacing[0] = ref->GetSpacing()[0];
    xAutoSpacing[1] = ref->GetSpacing()[1];
    xAutoSpacing[2] = ref->GetSpacing()[2];

    cout << "Spacing X="<< xAutoSpacing[0]
      << "  Spacing Y="<< xAutoSpacing[1]
      << "  Spacing Z="<< xAutoSpacing[2] << endl << endl;

    // Compute the resolution
    // Actually res[] is not used as resolution, but used as image size
    // by wook 4/6/2011
    res[0] = (int) bb[1][0];
    res[1] = (int) bb[1][1];
    res[2] = (int) bb[1][2];

    cout << "res X="<< res[0]
      << "  res Y="<< res[1]
      << "  res Z="<< res[2] << endl << endl;

    cout << "segColor = "<<segColor << endl;
    }
  else
    {
    cout << "   Image spacing     : " << xAutoSpacing[0] << " by " << xAutoSpacing[1] 
      << " by " << xAutoSpacing[2] << endl;
    cout << "   Margin            : " << xAutoMargin << endl;
    }

  // The real program begins here

  // Read the appropriate mesh type
  vtkPolyData *polyInput = ReadVTKData(fnInput.c_str());

  // Convert the model to triangles
  vtkTriangleFilter *tri = vtkTriangleFilter::New();
  tri->SetInputData(polyInput);
  cout << "Converting to triangles ..." << endl;
  tri->Update();
  vtkPolyData *pd = tri->GetOutput();

  // Get the extents of the data
  cout << "STL mesh bounds: " << endl;
  cout << "   X : " << pd->GetBounds()[0] << " to " << pd->GetBounds()[1] << endl;
  cout << "   Y : " << pd->GetBounds()[2] << " to " << pd->GetBounds()[3] << endl;
  cout << "   Z : " << pd->GetBounds()[4] << " to " << pd->GetBounds()[5] << endl;

  if((flagAutoSize) && (!useRef))
    {
    // Compute the bounding box and resolution automatically
    bb[1][0] = 2.0 * xAutoMargin + pd->GetBounds()[1] - pd->GetBounds()[0];
    bb[1][1] = 2.0 * xAutoMargin + pd->GetBounds()[3] - pd->GetBounds()[2];
    bb[1][2] = 2.0 * xAutoMargin + pd->GetBounds()[5] - pd->GetBounds()[4];
    bb[0][0] = pd->GetBounds()[0] - xAutoMargin;
    bb[0][1] = pd->GetBounds()[2] - xAutoMargin;
    bb[0][2] = pd->GetBounds()[4] - xAutoMargin;

    // Compute the resolution
    res[0] = (int) ceil(bb[1][0] / xAutoSpacing[0]);
    res[1] = (int) ceil(bb[1][1] / xAutoSpacing[1]);
    res[2] = (int) ceil(bb[1][2] / xAutoSpacing[2]);

    // Rescompute the size
    bb[1][0] = res[0] * xAutoSpacing[0];
    bb[1][1] = res[1] * xAutoSpacing[1];
    bb[1][2] = res[2] * xAutoSpacing[2];

    // Report the bounds
    cout << "   Image resolution  : " << res[0] << " by " << res[1] << " by " << res[2] << endl;
    cout << "   Image box origin  : {" << bb[0][0] << ", " << bb[0][1] << ", " << bb[0][2] << "}" << endl;
    cout << "   Image box size    : {" << bb[1][0] << ", " << bb[1][1] << ", " << bb[1][2] << "}" << endl << endl;
    }

  // Create a ITK image to store the results
  typedef itk::Image<unsigned char,3> ImageType;
  ImageType::Pointer img = ImageType::New();
  ImageType::RegionType region;


  if(useRef)
    {
    img->SetRegions(ref->GetBufferedRegion());
    img->SetOrigin(ref->GetOrigin());
    img->SetSpacing(ref->GetSpacing());
    img->SetDirection(ref->GetDirection());
    img->Allocate();
    img->FillBuffer(0);
    }
  else
    {
    // Allocate the image
    region.SetSize(0,res[0]); region.SetSize(1,res[1]); region.SetSize(2,res[2]);
    img->SetRegions(region);
    img->Allocate();
    img->FillBuffer(0);


    // Make it RAS!
    ImageType::DirectionType dir;
    dir(0,0) = -1; dir(1,1) = -1; dir(2,2) = 1;
    img->SetDirection(dir);

    // Set the origin and spacing of the image
    ImageType::SpacingType xSpacing;
    ImageType::PointType xOrigin;

    // Make origin RAS
    xOrigin[0] = -bb[0][0];
    xOrigin[1] = -bb[0][1];
    xOrigin[2] = bb[0][2];
    img->SetOrigin(xOrigin);

    // Set spacing
    for(unsigned int d = 0; d < 3; d++)
      {
      if(flagAutoSize)
        xSpacing[d] = xAutoSpacing[d];
      else
        xSpacing[d] = bb[1][d] / res[d];
      }
    img->SetSpacing(xSpacing);
    }

  // Create a vertex table from the polydata
  unsigned int nt = pd->GetNumberOfPolys(), it = 0;
  double **vtx = new double*[nt*3];  

  vtkCellArray *poly = pd->GetPolys();
  vtkIdType npts;
  const vtkIdType *pts;
  for(poly->InitTraversal();poly->GetNextCell(npts,pts);)
    {
    for(unsigned int i=0;i<3;i++)
      {
      double *x = pd->GetPoints()->GetPoint(pts[i]);
      vtx[it] = (double *) malloc(3*sizeof(double));

      // Use ITK to map from spatial to voxel coordinates
      Point<double,3> pt;
      for(unsigned int j=0;j<3;j++)
        if(j < 2)
          pt[j] = -x[j];
        else
          pt[j] = x[j];

      ContinuousIndex<double,3> idx;
      img->TransformPhysicalPointToContinuousIndex(pt, idx);
      for(unsigned int j=0;j<3;j++)
        vtx[it][j] = idx[j];

      ++it;
      }      
    }


  // Convert the polydata to an image
  if(doFloodFill)
    {
    cout << "Scan converting triangles and filling the interior ..." << endl;
    drawBinaryTrianglesFilled(img->GetBufferPointer(), res, vtx, nt, segColor);
    }
  else
    {
    cout << "Scan converting triangles ..." << endl;
    drawBinaryTrianglesSheetFilled(img->GetBufferPointer(), res, vtx, nt, segColor);
    }  

  cout << "zzz" << endl;

  cout << "ORIGIN: " << img->GetOrigin() << endl;

  // Save the image to disk
  typedef ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(img);
  writer->SetFileName(fnOutput.c_str());

  cout << "Writing the output image ..." << endl;
  writer->Update();
}
