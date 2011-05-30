#include "itkImageFileReader.h"
#include "itkOrientedRASImage.h"
#include "itkVTKImageExport.h"
#include "vtkImageImport.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkMarchingCubes.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include <vtkMatrix4x4.h>
#include "vnl/vnl_matrix_fixed.h"
#include "ReadWriteVTK.h"
#include <iostream>

#include "itk_to_nifti_xform.h"

using namespace std;

int usage()
{
  cout << "Usage: vtklevelset input.img output.vtk threshold" << endl;
  return -1;
}

template<class TImage>
void ConnectITKToVTK(itk::VTKImageExport<TImage> *fltExport,vtkImageImport *fltImport)
{
  fltImport->SetUpdateInformationCallback( fltExport->GetUpdateInformationCallback());
  fltImport->SetPipelineModifiedCallback( fltExport->GetPipelineModifiedCallback());
  fltImport->SetWholeExtentCallback( fltExport->GetWholeExtentCallback());
  fltImport->SetSpacingCallback( fltExport->GetSpacingCallback());
  fltImport->SetOriginCallback( fltExport->GetOriginCallback());
  fltImport->SetScalarTypeCallback( fltExport->GetScalarTypeCallback());
  fltImport->SetNumberOfComponentsCallback( fltExport->GetNumberOfComponentsCallback());
  fltImport->SetPropagateUpdateExtentCallback( fltExport->GetPropagateUpdateExtentCallback());
  fltImport->SetUpdateDataCallback( fltExport->GetUpdateDataCallback());
  fltImport->SetDataExtentCallback( fltExport->GetDataExtentCallback());
  fltImport->SetBufferPointerCallback( fltExport->GetBufferPointerCallback());
  fltImport->SetCallbackUserData( fltExport->GetCallbackUserData());
}

int main(int argc, char *argv[])
{
  if(argc != 4)
    return usage();

  // Read the input image
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[1]);
  fltReader->Update();
  ImageType::Pointer imgInput = fltReader->GetOutput();

  // Get the range of the input image
  float imax = imgInput->GetBufferPointer()[0];
  float imin = imax;
  for(size_t i = 0; i < imgInput->GetBufferedRegion().GetNumberOfPixels(); i++)
    {
    float x = imgInput->GetBufferPointer()[i];
    imax = std::max(imax, x);
    imin = std::min(imin, x);
    }

  float cut = atof(argv[3]);
  cout << "Image Range: [" << imin << ", " << imax << "]" << endl;
  cout << "Taking level set at " << cut << endl;

  // Create an importer and an exporter in VTK
  typedef itk::VTKImageExport<ImageType> ExporterType;
  ExporterType::Pointer fltExport = ExporterType::New();
  fltExport->SetInput(imgInput);
  vtkImageImport *fltImport = vtkImageImport::New();
  ConnectITKToVTK(fltExport.GetPointer(), fltImport);

  // Run marching cubes on the input image
  vtkMarchingCubes *fltMarching = vtkMarchingCubes::New();
  fltMarching->SetInput(fltImport->GetOutput());
  fltMarching->ComputeScalarsOff();
  fltMarching->ComputeGradientsOff();
  fltMarching->ComputeNormalsOn();
  fltMarching->SetNumberOfContours(1);
  fltMarching->SetValue(0,cut);
  fltMarching->Update();

  // Create the transform filter
  vtkTransformPolyDataFilter *fltTransform = vtkTransformPolyDataFilter::New();
  fltTransform->SetInput(fltMarching->GetOutput());
 
  // Compute the transform from VTK coordinates to NIFTI/RAS coordinates
  vnl_matrix_fixed<double, 4, 4> vtk2nii = 
    ConstructVTKtoNiftiTransform(
      imgInput->GetDirection().GetVnlMatrix(),
      imgInput->GetOrigin().GetVnlVector(),
      imgInput->GetSpacing().GetVnlVector());
  
  // Update the VTK transform to match
  vtkTransform *transform = vtkTransform::New();
  transform->SetMatrix(vtk2nii.data_block());
  fltTransform->SetTransform(transform);
  fltTransform->Update();

  // Get final output
  vtkPolyData *mesh = fltTransform->GetOutput();

  // Flip normals if determinant of SFORM is negative
  if(transform->GetMatrix()->Determinant() < 0)
    {
    vtkPointData *pd = mesh->GetPointData();
    vtkDataArray *nrm = pd->GetNormals();
    for(size_t i = 0; i < (size_t)nrm->GetNumberOfTuples(); i++)
      for(size_t j = 0; j < (size_t)nrm->GetNumberOfComponents(); j++)
        nrm->SetComponent(i,j,-nrm->GetComponent(i,j));
    nrm->Modified();
    }

  // Write the output
  WriteVTKData(mesh, argv[2]);
}
