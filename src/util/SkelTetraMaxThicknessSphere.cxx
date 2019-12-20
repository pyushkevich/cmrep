#include <vtkUnstructuredGrid.h>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

/**
 * This helper utility scans an image containing a radius map and 
 * a second vector image containing a center map, and finds the 
 * pixel with the maximum radius value. It then creates a vtkSphere
 * with that radius and center and saves it to VTK
 */
int main(int argc, char *argv[])
{
  if(argc < 4)
    {
    std::cerr << "Usage: " << argv[0] << " fn_img_rad fn_img_ctr fn_sphere" << std::endl;
    return -1;
    }

  const char *fn_img_rad = argv[1];
  const char *fn_img_ctr = argv[2];
  const char *fn_sphere = argv[3];

  typedef itk::VectorImage<double, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer r_rad = ReaderType::New();
  r_rad->SetFileName(fn_img_rad);
  r_rad->Update();

  ReaderType::Pointer r_ctr = ReaderType::New();
  r_ctr->SetFileName(fn_img_ctr);
  r_ctr->Update();

  typedef itk::ImageRegionConstIterator<ImageType> IterType;
  IterType it_rad(r_rad->GetOutput(), r_rad->GetOutput()->GetBufferedRegion());
  IterType it_ctr(r_ctr->GetOutput(), r_ctr->GetOutput()->GetBufferedRegion());

  double rad_max = 0.0;
  ImageType::PixelType ctr_max(r_ctr->GetOutput()->GetNumberOfComponentsPerPixel());
  for(; !it_rad.IsAtEnd(); ++it_rad, ++it_ctr)
    {
    double rad = it_rad.Get()[0];
    if(rad > rad_max)
      {
      ctr_max = it_ctr.Get();
      rad_max = rad;
      }
    }

  // Create a sphere to output
  vtkSphereSource* sphere = vtkSphereSource::New();
  sphere->SetCenter(ctr_max[0], ctr_max[1], ctr_max[2]);
  sphere->SetRadius(rad_max);
  sphere->SetThetaResolution(48);
  sphere->SetPhiResolution(96);
  sphere->Update();

  // Add a data array
  vtkPolyData *pSphere = sphere->GetOutput();
  vtkDoubleArray *daRad = vtkDoubleArray::New();
  daRad->SetNumberOfTuples(pSphere->GetNumberOfPoints());
  daRad->SetNumberOfComponents(1);
  daRad->FillComponent(0, rad_max);
  daRad->SetName("VoronoiRadius");
  pSphere->GetPointData()->AddArray(daRad);

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(pSphere);
  writer->SetFileName(fn_sphere);
  writer->Update();

  std::cout << "Max MIB Radius: " << rad_max << std::endl;
  std::cout << "Max MIB Center: " << ctr_max << std::endl;

}
