#include "ScriptInterface.h"
#include "ITKImageWrapper.h"
#include "itkImage.h"
#include "MedialAtomIterators.h"
#include "GenericMedialModel.h"
#include "MedialAtomGrid.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vtkCellLocator.h"
#include "vtkCell.h"
#include "vtkUnstructuredGrid.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"

vtkUnstructuredGrid * ExportVolumeMeshToVTK(
  GenericMedialModel *xModel, size_t nSamples);

using namespace medialpde;

int usage()
{
  printf(
    "Program: \n"
    "  cmrep_sample_img\n"
    "Description:\n"
    "  Sample image using correspondences in the cm-rep coordinate system\n"
    "Usage:\n"
    "  cmrep_sample_img trg_cmrep trg_ref_image src_cmrep src_image out_image mode\n"
    "Parameters:\n"
    "  trg_cmrep      CM-rep in the space of the output image\n"
    "  trg_ref_image  Image of same dimensions as the output image\n"
    "  src_cmrep      CM-rep in the space of the image to be sampled\n"
    "  src_image      Image to be sampled\n"
    "  out_image      Output image\n"
    "  mode           Mode: 0 for binary/label images, 1 for grayscale\n");

  return -1;
}

int main(int argc, char *argv[])
{
  // Sample an image using cm-reps
  if(argc != 7) return usage();

  // Get the mode
  int mode = atoi(argv[6]);

  // Read the source and target cm-reps
  MedialPDE cmr_target(argv[1]);
  MedialPDE cmr_source(argv[3]);

  // Read the source and target images
  typedef ITKImageWrapper<float> IWrapper;
  IWrapper *img_source = ITKImageWrapperFactory<float>::NewImageWrapper();
  IWrapper *img_target_ref = ITKImageWrapperFactory<float>::NewImageWrapper();
  img_source->LoadFromFile(argv[4]);
  img_target_ref->LoadFromFile(argv[2]);

  // vtkPolyData for holding all the cells in the cm-rep
  vtkUnstructuredGrid *ugSrc = 
    ExportVolumeMeshToVTK( cmr_source.GetMedialModel(), 4);
  vtkUnstructuredGrid *ugTrg = 
    ExportVolumeMeshToVTK( cmr_target.GetMedialModel(), 4);

  vtkCellLocator *loc = vtkCellLocator::New();
  loc->SetDataSet(ugTrg);
  loc->BuildLocator();

  double weights[1000];
  
  // Compute the bounds of the model
  double margin = 1.0;

  ugTrg->ComputeBounds();
  double xbnd[6]; // xmin, xmax, ymin, ymax, zmin, zmax;
  
  ugTrg->GetBounds(xbnd);
  xbnd[0] -= margin; xbnd[1] += margin;
  xbnd[2] -= margin; xbnd[3] += margin;
  xbnd[4] -= margin; xbnd[5] += margin;

  // Look up every voxel in the target image
  typedef itk::Image<float, 3> FloatImage;
  typedef itk::NearestNeighborInterpolateImageFunction<FloatImage,double> NNFunction;
  typedef itk::LinearInterpolateImageFunction<FloatImage,double> LIFunction;
  NNFunction::Pointer fnNN = NNFunction::New();
  fnNN->SetInputImage(img_source->GetInternalImage());
  LIFunction::Pointer fnLI = LIFunction::New();
  fnLI->SetInputImage(img_source->GetInternalImage());

  for(itk::ImageRegionIteratorWithIndex<FloatImage> iit(
      img_target_ref->GetInternalImage(),
      img_target_ref->GetInternalImage()->GetBufferedRegion());
    !iit.IsAtEnd(); ++iit)
    {
    // Get the physical point
    itk::Index<3> idx = iit.GetIndex();
    itk::Point<double, 3> ptRef;
    img_target_ref->GetInternalImage()->TransformIndexToPhysicalPoint(idx, ptRef);
    double x[3] = {ptRef[0], ptRef[1], ptRef[2]};

    iit.Set(0.0);

    // Check bounds
    if(
      xbnd[0] <= x[0] && xbnd[1] >= x[0] &&
      xbnd[2] <= x[1] && xbnd[3] >= x[1] &&
      xbnd[4] <= x[2] && xbnd[5] >= x[2])
      {
      // Find the nearest cell
      double xClosest[3];
      vtkIdType cellid;
      int subid;
      double dist2;
      loc->FindClosestPoint(x, xClosest, cellid, subid, dist2);
      if(dist2 < margin*margin)
        {
        // Find the weights of the point
        vtkCell *cell = ugTrg->GetCell(cellid);
        double pcoords[3];
        int inside = cell->EvaluatePosition(x, xClosest, subid, pcoords, dist2, weights);

        // Interpolate the look up position
        itk::Point<double, 3> ptSrc;
        ptSrc[0] = ptSrc[1] = ptSrc[2] = 0.0;
        for(size_t i = 0; i < cell->GetNumberOfPoints(); i++)
          {
          double w = weights[i];
          size_t k = cell->GetPointId(i);
          ptSrc[0] += w * ugSrc->GetPoint(k)[0];
          ptSrc[1] += w * ugSrc->GetPoint(k)[1];
          ptSrc[2] += w * ugSrc->GetPoint(k)[2];
          }

        // Compute the index
        itk::ContinuousIndex<double, 3> cidxSrc;
        img_source->GetInternalImage()->TransformPhysicalPointToContinuousIndex(ptSrc, cidxSrc);

        // Sample
        double val = 0.0;
        if(mode == 0)
          {
          if(fnNN->IsInsideBuffer(cidxSrc))
            val = fnNN->EvaluateAtContinuousIndex(cidxSrc);
          else
            val = 0.0;
          }
        else
          {
          if(fnLI->IsInsideBuffer(cidxSrc))
            val = fnLI->EvaluateAtContinuousIndex(cidxSrc);
          else
            val = 0.0;
          }

        // Store the pixel value
        iit.Set(val);
        }
      }
    }

  img_target_ref->SaveToFile(argv[5]);
}
