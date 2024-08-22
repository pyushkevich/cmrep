#include <iostream>
#include <string>

#include <vnl/vnl_inverse.h>

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkSubdivideTetra.h>
#include <itkImageFunction.h>
#include <itkImageFileReader.h>
#include <itkInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include "itkGaussianInterpolateImageFunction.h"
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include "ReadWriteVTK.h"

using namespace std;
using namespace itk;

int tet_sample_image_usage()
{
  cout << "tetsample - samples an image at points on a VTK mesh" << endl;
  cout << "usage: " << endl;
  cout << "   warpmesh [options] mesh.vtk out.vtk image scalarname" << endl;
  cout << "options: " << endl;
  cout << "   -i spec    :  Image coordinate specification. It can be " << endl;
  cout << "                 voxel space (ijk), scaled+offset voxel space (ijkos) " << endl;
  cout << "                 lps or ras. " << endl;
  cout << "   -m spec    : Mesh coordinate specification" << endl;
  cout << "   -r interp  : Interpolation Type (nn, linear, gauss)" << endl;
  cout << "   -s sigma   : Sigma for Gaussian interpolator" << endl;
  cout << "   -a alpha   : Alpha for Gaussian interpolator" << endl;
  return -1;
}

enum Coord { IJK, IJKOS, LPS, RAS, ANTS };
enum Interp { NN, LINEAR, GAUSS };

bool scan_coord(const char *text, Coord &val)
{
  if(!strcmp(text, "ijk")) val = IJK;
  else if(!strcmp(text, "ijkos")) val = IJKOS;
  else if(!strcmp(text, "lps")) val = LPS;
  else if(!strcmp(text, "ras")) val = RAS;
  else if(!strcmp(text, "ants")) val = ANTS;
  else return false;
  return true;
}

bool get_interp(const char *text, Interp &val)
{
  if(!strcmp(text, "nn")) val = NN;
  else if(!strcmp(text, "linear")) val = LINEAR;
  else if(!strcmp(text, "gauss")) val = GAUSS;
  else return false;
  return true;
}

bool get_interp_params(const char *text, double &val)
{
  if ((val=atof( text ))==0.0)
  return false;
  else return true;
}


/** 
 * This static function constructs a NIFTI matrix from the ITK direction
 * cosines matrix and Spacing and Origin vectors
 */
vnl_matrix_fixed<double,4,4> ConstructNiftiSform(
  vnl_matrix<double> m_dir, 
  vnl_vector<double> v_origin,
  vnl_vector<double> v_spacing)
{
  // Set the NIFTI/RAS transform
  vnl_matrix<double> m_ras_matrix;
  vnl_diag_matrix<double> m_scale, m_lps_to_ras;
  vnl_vector<double> v_ras_offset;

  // Compute the matrix
  m_scale.set(v_spacing);
  m_lps_to_ras.set(vnl_vector<double>(3, 1.0));
  m_lps_to_ras[0] = -1;
  m_lps_to_ras[1] = -1;
  m_ras_matrix = m_lps_to_ras * m_dir * m_scale;

  // Compute the vector
  v_ras_offset = m_lps_to_ras * v_origin;

  // Create the larger matrix
  vnl_vector<double> vcol(4, 1.0);
  vcol.update(v_ras_offset);

  vnl_matrix_fixed<double,4,4> m_sform;
  m_sform.set_identity();
  m_sform.update(m_ras_matrix);
  m_sform.set_column(3, vcol);
  return m_sform;
}

vnl_matrix_fixed<double,4,4> ConstructVTKtoNiftiTransform(
  vnl_matrix<double> m_dir, 
  vnl_vector<double> v_origin,
  vnl_vector<double> v_spacing)
{
  vnl_matrix_fixed<double,4,4> vox2nii = ConstructNiftiSform(m_dir, v_origin, v_spacing);
  vnl_matrix_fixed<double,4,4> vtk2vox; 
  vtk2vox.set_identity();
  for(size_t i = 0; i < 3; i++)
    {
    vtk2vox(i,i) = 1.0 / v_spacing[i];
    vtk2vox(i,3) = - v_origin[i] / v_spacing[i];
    }
  return vox2nii * vtk2vox;
}

struct TetSampleParam
{
  string fnMeshIn;
  string fnMeshOut;
  Coord mesh_coord, warp_coord;
  string fnImage;
  string scalarName;
  Interp interptype;
  double interpsigma;
  double interpalpha;
  TetSampleParam()
    {
    fnMeshIn = "";
    fnMeshOut = "";
    fnImage = "";
    scalarName = "";
    mesh_coord = RAS;
    warp_coord = RAS;
    interptype = NN;
    interpsigma = 3.0;
    interpalpha = 4.0;
    }
};

template <class X> X mymode(X *data, int size)
{
  int t, w;
  X md, oldmd;
  int count, oldcount;

  oldmd = 0;
  oldcount = 0;
  for(t=0; t<size; t++) {
    md = data[t];
    count = 1;
    for(w = t+1; w < size; w++) 
      if(md==data[w]) count++;
    if(count > oldcount) {
      oldmd = md;
      oldcount = count;
    }
  }
  return oldmd;
}

template <class ImageType>
class MeshImageSampler
{
public:
  typedef typename itk::InterpolateImageFunction<ImageType,double> FuncType;

  MeshImageSampler(FuncType *func, TetSampleParam &parm)
    {
    this->func = func;
    this->parm = parm;

    // Get the image
    const ImageType *sampim = func->GetInputImage();

    // Set up the transforms
    ijk2ras = ConstructNiftiSform(
      sampim->GetDirection().GetVnlMatrix().as_ref(),
      sampim->GetOrigin().GetVnlVector(),
      sampim->GetSpacing().GetVnlVector());

    vtk2ras = ConstructVTKtoNiftiTransform(
      sampim->GetDirection().GetVnlMatrix().as_ref(),
      sampim->GetOrigin().GetVnlVector(),
      sampim->GetSpacing().GetVnlVector());

    vnl_matrix_fixed<double, 4, 4> lps2ras;
    lps2ras.set_identity();
    lps2ras(0,0) = -1; lps2ras(1,1) = -1;

    ras2ijk = vnl_inverse(ijk2ras);
    ras2vtk = vnl_inverse(vtk2ras);
    ras2lps = vnl_inverse(lps2ras);

    // Store the active transform
    if(parm.mesh_coord == RAS)
      mesh2ras.set_identity();
    else if(parm.mesh_coord == LPS)
      mesh2ras = lps2ras;
    else if(parm.mesh_coord == IJKOS)
      mesh2ras = vtk2ras;
    else 
      mesh2ras = ijk2ras;

    }

  float SampleImage(double *xin) 
    {
    vnl_vector_fixed<double, 4> x_mesh, x_ras, x_ijk, v_ras;

    // Get the point (in whatever format that it's stored)
    x_mesh[0] = xin[0]; x_mesh[1] = xin[1]; x_mesh[2] = xin[2]; x_mesh[3] = 1.0;

    // Map the point into RAS coordinates
    x_ras = mesh2ras * x_mesh;

    // Map the point to IJK coordinates (continuous index)
    x_ijk = ras2ijk * x_ras;
    typename FuncType::ContinuousIndexType idx;
    idx[0] = x_ijk[0]; idx[1] = x_ijk[1]; idx[2] = x_ijk[2];

    // Interpolate the image at the point
    return func->EvaluateAtContinuousIndex(idx);
    }

private:
  typedef vnl_matrix_fixed<double, 4, 4> Mat44;
  Mat44 ijk2ras, vtk2ras, lps2ras, ras2ijk, ras2vtk, ras2lps, mesh2ras;
  TetSampleParam parm;
  typename FuncType::Pointer func;
};

/**
 * The actual method is templated over the VTK data type (unstructured/polydata)
 */
template <class TMeshType>
int TetSample(TetSampleParam &parm)
{
  // Read warp field
  typedef itk::Image<double, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::InterpolateImageFunction<ImageType> FuncType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NNFuncType;
  typedef itk::LinearInterpolateImageFunction<ImageType> LinearFuncType;
  typedef itk::GaussianInterpolateImageFunction<ImageType> GaussianFuncType;

  FuncType::Pointer func;
  if (parm.interptype == NN)
    {
    NNFuncType:: Pointer typedfunc = NNFuncType::New();
    func = typedfunc;
    }
  else if (parm.interptype == LINEAR) 
    {
    LinearFuncType:: Pointer typedfunc = LinearFuncType::New();
    func = typedfunc;
    }
  else if (parm.interptype == GAUSS) 
    {
    GaussianFuncType:: Pointer typedfunc = GaussianFuncType::New();
    double sigma[3];
    for (size_t i=0; i<3; i++)
      sigma[i] = parm.interpsigma;
    typedfunc->SetParameters(sigma, parm.interpalpha);
    cout << "Setting sigma=" << parm.interpsigma << " , alpha="<< parm.interpalpha << endl; 
    func = typedfunc;
    }
  else
    std::cerr << "Unknown interpolator" << std::endl;
  // Read each of the images in turn
  ImageType::Pointer sampim;
  ReaderType::Pointer reader;
  
  reader = ReaderType::New();
  reader->SetFileName(parm.fnImage.c_str());
  reader->Update();
  sampim = reader->GetOutput();
  func->SetInputImage(sampim);

  // Read the mesh
  TMeshType *mesh = ReadMesh<TMeshType>(parm.fnMeshIn.c_str());

  // Set up the mesh/image sampler
  MeshImageSampler<ImageType> sampler(func, parm);

  // Create the sample point arrays
  vtkFloatArray *sampledScalarPoint = vtkFloatArray::New();
  sampledScalarPoint->SetName(parm.scalarName.c_str());
  sampledScalarPoint->Allocate(mesh->GetNumberOfPoints());

  // Sample the image at vertices
  for(int k = 0; k < mesh->GetNumberOfPoints(); k++)
    {
    // Get the point (in whatever format that it's stored)
    float v_image = sampler.SampleImage(mesh->GetPoint(k));
    sampledScalarPoint->InsertNextValue(v_image);
    }

  // Create an array for cell integration
  vtkDoubleArray *sampledScalarCell = vtkDoubleArray::New();
  sampledScalarCell->SetName(parm.scalarName.c_str());
  sampledScalarCell->Allocate(mesh->GetNumberOfCells());

  // Temporary array for volume values
  vtkDoubleArray *tempvol = vtkDoubleArray::New();
  tempvol->Allocate(mesh->GetNumberOfCells());
  
  // Fill the array
  for(int i = 0; i < mesh->GetNumberOfCells(); i++)
    {
    tempvol->InsertNextValue(0);
    sampledScalarCell->InsertNextValue(0);
    }

  // Subdivide the tetrahedra in the mesh N times
  const size_t nsub = 3;
  TMeshType *subs[nsub];
  subs[0] = mesh;
  int celldiv = 1;
  for(size_t isub = 1; isub < nsub; isub++)
    {
    vtkSubdivideTetra *st = vtkSubdivideTetra::New();
    st->SetInputData(subs[isub-1]);
    st->Update();
    subs[isub] = st->GetOutput();
    celldiv *= 12;
    }

  // Get the finest-level mesh
  TMeshType *msub = subs[nsub-1];

  // Loop over subdivided cells
  for(size_t i = 0; i < (size_t)msub->GetNumberOfCells(); i++)
    {
    // Get the corners
    double x[4][3];
    vtkCell *cell = msub->GetCell(i);
    for(size_t j=0;j<4;j++)
      msub->GetPoint(cell->GetPointId(j), x[j]);

    // Compute cell volume
    double v = fabs(vtkTetra::ComputeVolume(x[0],x[1],x[2],x[3]));

    // Compute cell center
    vnl_vector_fixed<double, 3> tet_center;
    vtkTetra::TetraCenter(x[0],x[1],x[2],x[3],
      tet_center.data_block());

    // Sample image at the cell center
    // cout << "Tetcenter " << tet_center;
    float val = sampler.SampleImage(tet_center.data_block());
    // cout <<"; Value = " << val << endl;

    // Accumulate these values
    int target_id = i / celldiv;
    sampledScalarCell->SetTuple1(target_id, sampledScalarCell->GetTuple1(target_id) + val * v);
    tempvol->SetTuple1(target_id, tempvol->GetTuple1(target_id) + v);
    // cout << "Target " << target_id << " gets " << v << " and " << val * v << endl;
    }

  // Divide by the total volumes
  for(size_t i = 0; i < (size_t) mesh->GetNumberOfCells(); i++)
    {
    sampledScalarCell->SetTuple1(i, sampledScalarCell->GetTuple1(i) / tempvol->GetTuple1(i));
    }

  // Add both arrays
  mesh->GetCellData()->AddArray(sampledScalarCell);
  mesh->GetPointData()->AddArray(sampledScalarPoint);
  mesh->GetCellData()->SetActiveScalars(parm.scalarName.c_str());
  mesh->GetPointData()->SetActiveScalars(parm.scalarName.c_str());

  // Write the mesh
  WriteMesh<TMeshType>(mesh, parm.fnMeshOut.c_str());
  return 0;
}


  

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 4) return tet_sample_image_usage();

  // Parse the optional parameters
  int ch;
  TetSampleParam parm;
  while ((ch = getopt(argc, argv, "m:i:r:s:a:")) != -1)
    {
    switch(ch)
      {
    case 'm': if(!scan_coord(optarg, parm.mesh_coord) || parm.mesh_coord==ANTS)
                return tet_sample_image_usage();
              break;
    case 'i': if(!scan_coord(optarg, parm.warp_coord))
                return tet_sample_image_usage();
              break;
    case 'r': if(!get_interp(optarg, parm.interptype))
                return tet_sample_image_usage();
              break;
    case 's': if(!get_interp_params(optarg, parm.interpsigma))
                return tet_sample_image_usage();
              break;
    case 'a': if(!get_interp_params(optarg, parm.interpalpha))
                return tet_sample_image_usage();
              break;
    default: return tet_sample_image_usage();
      }
    }

  // Parse the filenames
  if(optind + 4 != argc) return tet_sample_image_usage();
  parm.fnMeshIn = argv[optind++];
  parm.fnMeshOut = argv[optind++];
  parm.fnImage = argv[optind++];
  parm.scalarName = argv[optind];

  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(parm.fnMeshIn.c_str());
  reader->OpenVTKFile();
  reader->ReadHeader();

  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    return TetSample<vtkUnstructuredGrid>(parm);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    cerr << "PolyData not supported, sorry!" << endl;
    // return TetSample<vtkPolyData>(parm);
    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}
