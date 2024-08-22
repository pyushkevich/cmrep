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
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <itkImageFileReader.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include "ReadWriteVTK.h"

using namespace std;
using namespace itk;

int warp_mesh_backwards_usage()
{
  cout << "warpmesh - Applies a warp field (Brian's format) to a VTK mesh" << endl;
  cout << "usage: " << endl;
  cout << "   warpmesh [options] mesh.vtk out.vtk warp_images" << endl;
  cout << "   warpmesh [options] mesh.vtk out.vtk affine_matrix" << endl;
  cout << "options: " << endl;
  cout << "   -w spec    : Warp coordinate specification. The warp field gives " << endl;
  cout << "                 a displacement in some coordinate space. It can be " << endl;
  cout << "                 voxel space (ijk), scaled+offset voxel space (ijkos) " << endl;
  cout << "                 lps or ras. Also, 'ants' for ANTS warps " << endl;
  cout << "   -m spec    : Mesh coordinate specification" << endl;
  return -1;
}

enum Coord { IJK, IJKOS, LPS, RAS, ANTS };

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

struct WarpMeshParam
{
  string fnMeshIn;
  string fnMeshOut;
  Coord mesh_coord, warp_coord;
  string fnWarp[3];
  string fnAffine;
  bool flagWarp;
  WarpMeshParam()
    {
    fnMeshIn = "";
    fnMeshOut = "";
    fnAffine = "";
    for(size_t d = 0; d < 3; d++)
      fnWarp[d] = "";
    flagWarp = true;
    mesh_coord = RAS;
    warp_coord = RAS;
    }
};


template<class TMeshType>
int ComputeVolumeOrArea(TMeshType *mesh, vtkFloatArray *cellarray, vtkFloatArray *pointarray, double &totalvol)
{
  int celltype;
  totalvol = 0.0;
  cellarray->Allocate(mesh->GetNumberOfCells());
  pointarray->Allocate(mesh->GetNumberOfPoints());
  for(int j = 0; j < mesh->GetNumberOfPoints(); j++)
     pointarray->InsertNextValue(0.0);
  for(int j = 0; j < mesh->GetNumberOfCells(); j++)
    {
    vtkCell *cell =  mesh->GetCell(j);
    if(cell->GetCellType() == VTK_TETRA)
      {
      double p[4][3];
      mesh->GetPoint(cell->GetPointId(0), p[0]);
      mesh->GetPoint(cell->GetPointId(1), p[1]);
      mesh->GetPoint(cell->GetPointId(2), p[2]);
      mesh->GetPoint(cell->GetPointId(3), p[3]);
      double vol = vtkTetra::ComputeVolume(p[0], p[1], p[2], p[3]);
      cellarray->InsertNextValue(vol);
      totalvol += vol;
      for(size_t ivert = 0; ivert < 4; ivert++)
         pointarray->SetTuple1(cell->GetPointId(ivert), pointarray->GetTuple1(cell->GetPointId(ivert)) + vol/4.0);
      celltype = VTK_TETRA;

      }
    else if (cell->GetCellType() == VTK_TRIANGLE)
      {
      double p[3][3];
      mesh->GetPoint(cell->GetPointId(0), p[0]);
      mesh->GetPoint(cell->GetPointId(1), p[1]);
      mesh->GetPoint(cell->GetPointId(2), p[2]);
      double vol = vtkTriangle::TriangleArea(p[0], p[1], p[2]);
      cellarray->InsertNextValue(vol);
      totalvol += vol;
      for(size_t ivert = 0; ivert < 3; ivert++)
         pointarray->SetTuple1(cell->GetPointId(ivert), pointarray->GetTuple1(cell->GetPointId(ivert)) + vol/3.0);

      celltype = VTK_TRIANGLE;
      }
    else return VTK_EMPTY_CELL;
    }

  return celltype;
}



void ReadMatrix(const char *fname, itk::Matrix<double,4,4> &mat)
  {
  ifstream fin(fname);
  for(size_t i = 0; i < 4; i++)
    for(size_t j = 0; j < 4; j++)
      if(fin.good())
        {
        fin >> mat[i][j];
        }
      else
        {
        throw itk::ExceptionObject("Unable to read matrix");
        }
  fin.close();
  }



/**
 * The actual method is templated over the VTK data type (unstructured/polydata)
 */
template <class TMeshType>
int WarpMesh(WarpMeshParam &parm)
{
  // Read warp field
  typedef itk::Image<double, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::LinearInterpolateImageFunction<ImageType> FuncType;

  // Read each of the images in turn
  ImageType::Pointer warp[3];
  ReaderType::Pointer reader[3];
  FuncType::Pointer func[3];
  itk::Matrix<double,4,4> affine;

  vnl_matrix_fixed<double, 4, 4> ijk2ras, vtk2ras, lps2ras;
  ijk2ras.set_identity();
  vtk2ras.set_identity();
  lps2ras.set_identity();
  lps2ras(0,0) = -1; lps2ras(1,1) = -1;

  if(parm.flagWarp)
    {
    // Read in the warps
    for(size_t d = 0; d < 3; d++)
      {
      reader[d] = ReaderType::New();
      reader[d]->SetFileName(parm.fnWarp[d].c_str());
      reader[d]->Update();
      warp[d] = reader[d]->GetOutput();
      func[d] = FuncType::New();
      func[d]->SetInputImage(warp[d]);
      }

    // Set up the transforms
    ijk2ras = ConstructNiftiSform(
      warp[0]->GetDirection().GetVnlMatrix().as_ref(),
      warp[0]->GetOrigin().GetVnlVector(),
      warp[0]->GetSpacing().GetVnlVector());

    vtk2ras = ConstructVTKtoNiftiTransform(
      warp[0]->GetDirection().GetVnlMatrix().as_ref(),
      warp[0]->GetOrigin().GetVnlVector(),
      warp[0]->GetSpacing().GetVnlVector());
    }
  else
    {
    ReadMatrix(parm.fnAffine.c_str(), affine);
    }

  // Read the mesh
  TMeshType *mesh = ReadMesh<TMeshType>(parm.fnMeshIn.c_str());


  vnl_matrix_fixed<double, 4, 4> ras2ijk = vnl_inverse(ijk2ras);
  vnl_matrix_fixed<double, 4, 4> ras2vtk = vnl_inverse(vtk2ras);
  vnl_matrix_fixed<double, 4, 4> ras2lps = vnl_inverse(lps2ras);

  cout << "RAS transform " << endl;
  cout << ijk2ras << endl;

  // Create the volume array (cell-wise) for future jacobian computation
  vtkFloatArray *jacobian = vtkFloatArray::New();
  jacobian->SetName("Cell Jacobian");
  vtkFloatArray *pointjacobian = vtkFloatArray::New();
  pointjacobian->SetName("Point Jacobian");
  double total_old_vol = 0;
  int celltype = ComputeVolumeOrArea(mesh, jacobian, pointjacobian, total_old_vol);

  // Update the coordinates
  for(int k = 0; k < mesh->GetNumberOfPoints(); k++)
    {
    // Get the point (in whatever format that it's stored)
    vnl_vector_fixed<double, 4> x_mesh, x_ras, x_ijk, v_warp, v_ras;
    vnl_vector_fixed<double, 4> y_ras, y_mesh;
    x_mesh[0] = mesh->GetPoint(k)[0]; x_mesh[1] = mesh->GetPoint(k)[1]; x_mesh[2] = mesh->GetPoint(k)[2];
    x_mesh[3] = 1.0;

    // Map the point into RAS coordinates
    if(parm.mesh_coord == RAS)
      x_ras = x_mesh;
    else if(parm.mesh_coord == LPS)
      x_ras = lps2ras * x_mesh;
    else if(parm.mesh_coord == IJKOS)
      x_ras = vtk2ras * x_mesh;
    else 
      x_ras = ijk2ras * x_mesh;

    if(parm.flagWarp)
      {
      // Map the point to IJK coordinates (continuous index)
      x_ijk = ras2ijk * x_ras;
      FuncType::ContinuousIndexType idx;
      idx[0] = x_ijk[0]; idx[1] = x_ijk[1]; idx[2] = x_ijk[2];

      // Interpolate the warp at the point
      // cout << "Evaluate at index " << idx[0] << " " << idx[1] << " " << idx[2] << endl;
      for(size_t d = 0; d < 3; d++)
        v_warp[d] = func[d]->EvaluateAtContinuousIndex(idx);
      v_warp[3] = 0.0;

      // Compute the displacement in RAS coordinates
      if(parm.warp_coord == RAS)
        v_ras = v_warp;
      else if(parm.warp_coord == LPS)
        v_ras = lps2ras * v_warp;
      else if(parm.warp_coord == IJKOS)
        v_ras = vtk2ras * v_warp;
      else if(parm.warp_coord == ANTS)
        {
        // vector is multiplied by the direction matrix ??? Really???
        v_ras = lps2ras * v_warp;
        }
      else 
        v_ras = ijk2ras * v_warp;

      // Add displacement
      y_ras = x_ras + v_ras;
      }
    else
      {
      vnl_matrix_fixed<double, 4,4> M = affine.GetVnlMatrix();
      y_ras = M * x_ras;
      }

    // Map new coordinate to desired system
    if(parm.mesh_coord == RAS)
      y_mesh = y_ras;
    else if(parm.mesh_coord == LPS)
      y_mesh = ras2lps * y_ras;
    else if(parm.mesh_coord == IJKOS)
      y_mesh = ras2vtk * y_ras;
    else 
      y_mesh = ras2ijk * y_ras;

    // Map the transformation to RAS
    // vnl_vector_fixed<double, 4> w_ras = 
    //  wimg[0]->GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix() * w;

    // vnl_vector_fixed<double, 4> w_ras = w;
      
    /*
     * OLD CODE THT WORKED WITH ANTS
    // Assume transformation is in spacing units
    itk::FixedArray<double, 3> aw, awras;
    aw[0] = w[0];
    aw[1] = w[1];
    aw[2] = w[2];
    wimg[0]->TransformLocalVectorToPhysicalVector(aw, awras);
    p->GetPoints()->SetPoint(k, pt[0] + awras[0], pt[1] + awras[1], pt[2] + awras[2]);
    */

    mesh->GetPoints()->SetPoint(k, y_mesh[0], y_mesh[1], y_mesh[2]);
    }
    
  // Create the volume array (cell-wise) for future jacobian computation
  if(celltype != VTK_EMPTY_CELL)
    {
    // Compute the updated volumes
    vtkFloatArray *newvol = vtkFloatArray::New();
    vtkFloatArray *newpointvol = vtkFloatArray::New();
    if (celltype == VTK_TETRA)
    {
      newvol->SetName("Cell Volume");
      newpointvol->SetName("Point Volume");
    }
    else // VTK_TRIANGLE
    {
      newvol->SetName("Cell Area");
      newpointvol->SetName("Point Area");
    }
    double total_new_vol = 0;
    ComputeVolumeOrArea(mesh, newvol, newpointvol, total_new_vol);

    // Compute the Jacobian
    for(int i = 0; i < mesh->GetNumberOfCells(); i++)
      jacobian->SetTuple1(i, newvol->GetTuple1(i) / jacobian->GetTuple1(i));
    for(int i = 0; i < mesh->GetNumberOfPoints(); i++)
      pointjacobian->SetTuple1(i, newpointvol->GetTuple1(i) / pointjacobian->GetTuple1(i));

    // Add both arrays
    mesh->GetCellData()->AddArray(jacobian);
    mesh->GetCellData()->AddArray(newvol);
    mesh->GetPointData()->AddArray(pointjacobian);
    mesh->GetPointData()->AddArray(newpointvol);

    // Report change in volume/area
    if (celltype == VTK_TETRA)
      printf("VOLUME_STATS: %12.8f %12.8f %12.8f\n", 
        total_old_vol, total_new_vol, total_new_vol/total_old_vol); 
    else
      printf("AREA_STATS: %12.8f %12.8f %12.8f\n", 
        total_old_vol, total_new_vol, total_new_vol/total_old_vol); 
    }

  // Write the mesh
  WriteMesh<TMeshType>(mesh, parm.fnMeshOut.c_str());
  return 0;
}


  

int main(int argc, char **argv)
{
  // Check the parameters
  if(argc < 4) return warp_mesh_backwards_usage();

  // Parse the optional parameters
  int ch;
  WarpMeshParam parm;
  while ((ch = getopt(argc, argv, "m:w:")) != -1)
    {
    switch(ch)
      {
    case 'm': if(!scan_coord(optarg, parm.mesh_coord) || parm.mesh_coord==ANTS)
                return warp_mesh_backwards_usage();
              break;
    case 'w': if(!scan_coord(optarg, parm.warp_coord))
                return warp_mesh_backwards_usage();
              break;
    default: return warp_mesh_backwards_usage();
      }
    }

  // Parse the filenames
  if(optind + 5 != argc && optind + 3 != argc) return warp_mesh_backwards_usage();
  parm.fnMeshIn = argv[optind++];
  parm.fnMeshOut = argv[optind++];
  if(optind + 3 == argc)
    {
    parm.flagWarp = true;
    for(size_t d = 0; d < 3; d++)
      parm.fnWarp[d] = argv[optind++];
    }
  else
    {
    parm.flagWarp = false;
    parm.fnAffine = argv[optind++];
    }

  // Check the data type of the input file
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(parm.fnMeshIn.c_str());
  reader->OpenVTKFile();
  reader->ReadHeader();

  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    return WarpMesh<vtkUnstructuredGrid>(parm);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return WarpMesh<vtkPolyData>(parm);
    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}
