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
#include "itkLinearInterpolateImageFunction.h"
#include "vtkImplicitFunction.h"
#include "vtkCleanPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkCell.h"
#include "vtkContourFilter.h"
#include <vtkDiscreteMarchingCubes.h>
#include <vtkUnsignedShortArray.h>
#include <vtkAppendPolyData.h>

#include "MeshTraversal.h"

#include "itk_to_nifti_xform.h"
#include <vnl/vnl_inverse.h>

#ifdef CMREP_USE_VCG
// stuff to define the mesh
#include "VCGTriMesh.h"

// local optimization
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#endif

using namespace std;

int usage()
{
  cout << "Usage: vtklevelset [options] input.img output.vtk threshold" << endl;
  cout << "Options: " << endl;
  cout << "  -c clipimage.img   : clip output (keep part where clipimage > 0" << endl;
  cout << "  -v                 : export mesh in voxel coordinates (not physical)" << endl;
  cout << "  -k                 : apply clean filter to the mesh" << endl;
  cout << "  -d                 : perform Delaunay edge flipping" << endl;
  cout << "  -2                 : use 2D algorithm" << endl;
  cout << "  -pl                : preserve labels" << endl; // issue #4: added label preserving option
  cout << "  -r <factor|num>    : reduce the size of the final mesh via quadric edge collapse" << endl;
  cout << "                       algorithm from MeshLab by a factor (e.g., 0.1) or to given " << endl;
  cout << "                       number of vertices " << endl;
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

template<class TImage>
class ClipFunction : public vtkImplicitFunction 
{
public:

  vtkTypeMacro(ClipFunction, vtkImplicitFunction);

  // This is used by the clip code
  double EvaluateFunction(double x[3]);

  // This is not really used in the clip code
  void EvaluateGradient(double x[3], double g[3]) {}

  // Set the image
  void SetImage(TImage *image);

private:

  typedef itk::LinearInterpolateImageFunction<TImage> func;



};

int main(int argc, char *argv[])
{
  if(argc < 4)
    return usage();

  // Clip image
  const char *imClip = NULL;
  bool voxelSpace = false;
  bool preserveLabels = false;
  bool flag_clean = false, flag_delaunay = false, flag_2d = false;
  double reduction_factor = 1.0;
  for(int i = 1; i < argc - 3; i++)
    {
    if(!strcmp(argv[i], "-c"))
      {
      imClip = argv[++i];
      }
    else if(!strcmp(argv[i], "-v"))
      {
      voxelSpace = true;
      }
    else if(!strcmp(argv[i], "-d"))
      {
      flag_delaunay = true;
      }
    else if(!strcmp(argv[i], "-k"))
      {
      flag_clean = true;
      }
    else if(!strcmp(argv[i], "-2"))
      {
      flag_2d = true;
      }
    else if(!strcmp(argv[i], "-pl"))
      {
      preserveLabels = true;
      }
    else if(!strcmp(argv[i], "-r"))
      {
      reduction_factor = atof(argv[++i]);
      }
    else
      { 
      cerr << "Unknown option " << argv[i] << endl;
      return -1;
      }
    }

  // Read the input image
  typedef itk::OrientedRASImage<float, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer fltReader = ReaderType::New();
  fltReader->SetFileName(argv[argc-3]);
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

  float cut = atof(argv[argc-1]);
  cout << "Image Range: [" << imin << ", " << imax << "]" << endl;
  cout << "Taking level set at " << cut << endl;

  // Create an importer and an exporter in VTK
  typedef itk::VTKImageExport<ImageType> ExporterType;
  ExporterType::Pointer fltExport = ExporterType::New();
  fltExport->SetInput(imgInput);
  vtkImageImport *fltImport = vtkImageImport::New();
  ConnectITKToVTK(fltExport.GetPointer(), fltImport);

  // Run marching cubes on the input image
  vtkPolyData *pipe_tail;
  if(preserveLabels)
    {
      std::cout << "Preserving Labels Mode" << std::endl;

      // Append filter for assembling labels
      vtkAppendPolyData *fltAppend = vtkAppendPolyData::New();

      // Extracting one label at a time and assigning label value
      for (float i = cut; i <= imax; i += 1.0)
        {
        float lbl = floor(i);

        std::cout << "  -- Processing Label: " << lbl << std::endl;

        // Extract one label
        vtkDiscreteMarchingCubes *fltDMC = vtkDiscreteMarchingCubes::New();
        fltDMC->SetInputConnection(fltImport->GetOutputPort());
        fltDMC->ComputeGradientsOff();
        fltDMC->ComputeScalarsOff();
        fltDMC->SetNumberOfContours(1);
        fltDMC->ComputeNormalsOn();
        fltDMC->SetValue(0, lbl);
        fltDMC->Update();

        vtkPolyData *labelMesh = fltDMC->GetOutput();

        // Set scalar values for the label
        vtkUnsignedShortArray *scalar = vtkUnsignedShortArray::New();
        scalar->SetNumberOfComponents(1);
        for (vtkIdType i = 0; i < labelMesh->GetNumberOfPoints(); ++i)
          {
          scalar->InsertNextTuple1(lbl);
          }
        scalar->SetName("Label");
        labelMesh->GetPointData()->SetScalars(scalar);
        fltAppend->AddInputData(labelMesh);
        }

    fltAppend->Update();
    pipe_tail = fltAppend->GetOutput();
    }
  else if(!flag_2d)
    {
    vtkMarchingCubes *fltMarching = vtkMarchingCubes::New();
    fltMarching->SetInputConnection(fltImport->GetOutputPort());
    fltMarching->ComputeScalarsOff();
    fltMarching->ComputeGradientsOff();
    fltMarching->ComputeNormalsOn();
    fltMarching->SetNumberOfContours(1);
    fltMarching->SetValue(0,cut);
    fltMarching->Update();
    pipe_tail = fltMarching->GetOutput();
    }
  else
    {
    vtkContourFilter *fltContour = vtkContourFilter::New();
    fltContour->SetInputConnection(fltImport->GetOutputPort());
    fltContour->ComputeScalarsOff();
    fltContour->ComputeGradientsOff();
    fltContour->ComputeNormalsOn();
    fltContour->SetNumberOfContours(1);
    fltContour->SetValue(0,cut);
    fltContour->Update();
    pipe_tail = fltContour->GetOutput();
    }

  // If the clean option is requested, use it
  if(flag_clean || flag_delaunay) 
    {
    vtkTriangleFilter *trifi = vtkTriangleFilter::New();
    trifi->SetInputData(pipe_tail);
    trifi->PassLinesOff();
    trifi->PassVertsOff();
    trifi->Update();

    vtkCleanPolyData *clean = vtkCleanPolyData::New();
    clean->SetInputConnection(trifi->GetOutputPort());
    clean->PointMergingOn();
    clean->SetTolerance(0.0);
    clean->Update();

    vtkTriangleFilter *trifi2 = vtkTriangleFilter::New();
    trifi2->SetInputConnection(clean->GetOutputPort());
    trifi2->PassLinesOff();
    trifi2->PassVertsOff();
    trifi2->Update();
    pipe_tail = trifi2->GetOutput();
    }

  if(flag_delaunay)
    {

    typedef vnl_vector_fixed<double, 3> Vec;
    Vec *X = new Vec[pipe_tail->GetNumberOfPoints()];
    for(size_t i = 0; i < pipe_tail->GetNumberOfPoints(); i++)
      for(size_t k = 0; k < 3; k++)
        X[i][k] = pipe_tail->GetPoint(i)[k];

    // Generate the mesh object
    TriangleMesh bmesh;
    TriangleMeshGenerator gen(&bmesh, pipe_tail->GetNumberOfPoints());
    for(int i = 0; i < pipe_tail->GetNumberOfCells(); i++)
      {
      vtkCell *c = pipe_tail->GetCell(i);
      gen.AddTriangle(c->GetPointId(0), c->GetPointId(1), c->GetPointId(2));
      }
    gen.GenerateMesh();

    // Perform Delaunay
    bmesh.MakeDelaunay(X);

    // Set new cell data on the mesh
    pipe_tail->DeleteCells();
    pipe_tail->Allocate(bmesh.triangles.size());
    for(int i = 0; i < bmesh.triangles.size(); i++)
      {
      Triangle &tri = bmesh.triangles[i];
      vtkIdType pts[3];
      pts[0] = tri.vertices[0];
      pts[1] = tri.vertices[1];
      pts[2] = tri.vertices[2];
      pipe_tail->InsertNextCell(VTK_TRIANGLE, 3, pts);
      }

    pipe_tail->BuildCells();
    pipe_tail->BuildLinks();
    }
  
  // Create the transform filter
  vtkTransformPolyDataFilter *fltTransform = vtkTransformPolyDataFilter::New();
  fltTransform->SetInputData(pipe_tail);
 
  // Compute the transform from VTK coordinates to NIFTI/RAS coordinates
  typedef vnl_matrix_fixed<double, 4, 4> Mat44;
  Mat44 vtk2out;
  Mat44 vtk2nii = ConstructVTKtoNiftiTransform(
    imgInput->GetDirection().GetVnlMatrix().as_ref(),
    imgInput->GetOrigin().GetVnlVector(),
    imgInput->GetSpacing().GetVnlVector());

  // If we actually asked for voxel coordinates, we need to fix that
  if(voxelSpace)
    {
    Mat44 vox2nii = ConstructNiftiSform(
      imgInput->GetDirection().GetVnlMatrix().as_ref(),
      imgInput->GetOrigin().GetVnlVector(),
      imgInput->GetSpacing().GetVnlVector());
    
    Mat44 nii2vox = vnl_inverse(vox2nii);
    vtk2out = nii2vox * vtk2nii;
    }
  else
    {
    vtk2out = vtk2nii;
    }
  
  // Update the VTK transform to match
  vtkTransform *transform = vtkTransform::New();
  transform->SetMatrix(vtk2out.data_block());
  fltTransform->SetTransform(transform);
  fltTransform->Update();

  // Get final output
  vtkSmartPointer<vtkPolyData> mesh = fltTransform->GetOutput();

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

  if(reduction_factor != 1.0)
    {
#ifdef CMREP_USE_VCG
    // Convert to a VCG complex
    VCGTriMesh tri_mesh;
    tri_mesh.ImportFromVTK(mesh);
    tri_mesh.CleanMesh();

    vcg::tri::TriEdgeCollapseQuadricParameter qparams;
    qparams.QualityThr = 0.3;
    qparams.PreserveBoundary=true;
    qparams.BoundaryQuadricWeight=1.0;
    qparams.PreserveTopology=true;
    qparams.QualityWeight=false;
    qparams.NormalCheck=true;
    qparams.OptimalPlacement=true;
    qparams.QualityQuadric=true;
    qparams.QualityQuadricWeight=0.001;

    // decimator initialization
    vcg::LocalOptimization<VCGTriMesh::Mesh> DeciSession(tri_mesh.GetMesh(), &qparams);

    // specialization
    typedef vcg::tri::BasicVertexPair<VCGTriMesh::Vertex> VertexPair;
    class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric<
        VCGTriMesh::Mesh, VertexPair,
        MyTriEdgeCollapse,
        vcg::tri::QInfoStandard<VCGTriMesh::Vertex>  >
    {
    public:
      typedef vcg::tri::TriEdgeCollapseQuadric<
          VCGTriMesh::Mesh, VertexPair,
          MyTriEdgeCollapse,
          vcg::tri::QInfoStandard<VCGTriMesh::Vertex>  > TECQ;
      typedef  VCGTriMesh::Mesh::VertexType::EdgeType EdgeType;
      inline MyTriEdgeCollapse(const VertexPair &p, int i, vcg::BaseParameterClass *pp) :TECQ(p,i,pp){}
    };

    // Target number of vertices
    int n_target = reduction_factor < 1.0
        ? (int) (tri_mesh.GetMesh().FN() * reduction_factor)
        : (int) reduction_factor;

    int t1=clock();
    DeciSession.Init<MyTriEdgeCollapse>();
    int t2=clock();
    printf("Simplifying to reduce from %7i to %7i faces\n", tri_mesh.GetMesh().FN(), n_target);

    DeciSession.SetTargetSimplices(n_target);
    DeciSession.SetTimeBudget(2.0f);
    // DeciSession.SetTargetOperations(100000);
    // if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);

    double TargetError = std::numeric_limits<double>::max();
    while( DeciSession.DoOptimization() && tri_mesh.GetMesh().FN() > n_target && DeciSession.currMetric < TargetError)
      printf("Current Mesh size %7i heap sz %9i err %9g \n", tri_mesh.GetMesh().FN(), int(DeciSession.h.size()),DeciSession.currMetric);

    int t3=clock();
    printf("mesh (%d,%d) Error %g \n", tri_mesh.GetMesh().VN(), tri_mesh.GetMesh().FN(), DeciSession.currMetric);
    printf("\nCompleted in (%5.3f+%5.3f) sec\n", float(t2-t1)/CLOCKS_PER_SEC, float(t3-t2)/CLOCKS_PER_SEC);

    DeciSession.Finalize<MyTriEdgeCollapse>();

    // Postprocess by cleaning and fixing normals
    tri_mesh.CleanMesh();
    tri_mesh.RecomputeNormals();

    tri_mesh.ExportToVTK(mesh);
#else
    cout << "Mesh reduction is not implemented - need to compile with VCG" << endl;
    return -1;
#endif
    }

  // Write the output
  WriteVTKData(mesh, argv[argc-2]);
}
