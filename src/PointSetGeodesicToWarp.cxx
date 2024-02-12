#include "PointSetHamiltonianSystem.h"
#include "lddmm_data.h"
#include <iostream>
#include "util/ReadWriteVTK.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "FastLinearInterpolator.h"
#include "itkMultiThreaderBase.h"

#include <cstdarg>

using namespace std;


int usage()
{
  cout << "lmtowarp: Generate warp field from landmark geodesic shooting" << endl;
  cout << "Usage:" << endl;
  cout << "  lmtowarp [options]" << endl;
  cout << "Required options:" << endl;
  cout << "  -r image.nii       : Reference image (defines warp space)" << endl;
  cout << "  -o warp.nii        : Output deformation field" << endl;
  cout << "  -m mesh.vtk        : Mesh with the InitialMomentum array defined" << endl;
  cout << "  -s sigma           : kernel standard deviation" << endl;
  cout << "Additional options:" << endl;
  cout << "  -d dim             : problem dimension (3)" << endl;
  cout << "  -n N               : number of time steps (100)" << endl;
  cout << "  -g mask.nii        : limit warp computation to a masked region" << endl;
  cout << "  -B                 : use brute force warp computation (not approximation)" << endl;
  cout << "  -t n_threads       : limit number of concurrent threads to n_threads" << endl;
  cout << "Mesh warping mode:" << endl;
  cout << "  -M in.vtk out.vtk  : additional meshes to apply warp to" << endl;
  cout << "                       when using -M, -r/-o are optional" << endl;
  cout << "Animation:" << endl;
  cout << "  -a k               : save an animation frame every k-th timepoint" << endl;
  cout << "                       output files must have %03d pattern in them" << endl;

  return -1;
}

void check(bool condition, const char *format,...)
{
  if(!condition)
    {
    char buffer[256];
    va_list args;
    va_start (args, format);
    vsprintf (buffer,format, args);
    va_end (args);

    cerr << buffer << endl;
    exit(-1);
    }
}


struct WarpGenerationParameters
{
  typedef std::pair<string, string> MeshPair;

  string fnReference, fnMesh, fnOutWarp, fnMask;
  double sigma;
  unsigned int dim;
  unsigned int N;
  unsigned int anim_freq;
  unsigned int n_threads;
  bool brute;

  std::list<MeshPair> fnWarpMeshes;

  WarpGenerationParameters():
    sigma(0.0), dim(3), N(100), anim_freq(0), n_threads(0), brute(false) {}
};

template <class TPixel, unsigned int VDim>
class PointSetGeodesicToWarp
{
public:
  typedef PointSetHamiltonianSystem<double, VDim> HSystem;
  typedef typename HSystem::Vector Vector;
  typedef typename HSystem::Matrix Matrix;


  typedef LDDMMData<TPixel, VDim> LDDMM;
  typedef typename LDDMM::ImageType ImageType;
  typedef typename LDDMM::VectorImageType VectorImageType;
  typedef typename LDDMM::ImagePointer ImagePointer;
  typedef typename LDDMM::VectorImagePointer VectorImagePointer;

  typedef typename ImageType::PointType PointType;
  typedef FastLinearInterpolator<VectorImageType,TPixel,VDim> FastInterpolator;
  typedef itk::ContinuousIndex<TPixel,VDim> ContIndexType;

  static int run(const WarpGenerationParameters &param);

  static VectorImagePointer brute_force_method(
      const WarpGenerationParameters &param,
      HSystem &hsys, int tdir, int tStart, int tEnd);

private:

  static void UpdateAndWriteMesh(
    vtkPolyData *mesh, const Matrix &x, const Matrix &v, const Matrix &x0,
    const std::string filePattern,
    int k);
};


template <class TPixel, unsigned int VDim>
void
PointSetGeodesicToWarp<TPixel, VDim>
::UpdateAndWriteMesh(
  vtkPolyData *mesh, const Matrix &x, const Matrix &v, const Matrix &x0,
  const std::string filePattern,
  int k)
{
  // Velocity array
  vtkDoubleArray *arr_v = vtkDoubleArray::New();
  arr_v->SetNumberOfComponents(VDim);
  arr_v->SetNumberOfTuples(v.rows());
  arr_v->SetName("Velocity");
  mesh->GetPointData()->AddArray(arr_v);

  // Initial position array
  vtkDoubleArray *arr_x0 = vtkDoubleArray::New();
  arr_x0->SetNumberOfComponents(VDim);
  arr_x0->SetNumberOfTuples(v.rows());
  arr_x0->SetName("InitialPosition");
  mesh->GetPointData()->AddArray(arr_x0);

  // Assign new points to the mesh
  for(int i = 0; i < x.rows(); i++)
    {
    double x_out[VDim];
    for(int a = 0; a < VDim; a++) 
      {
      x_out[a] = x(i,a);
      arr_v->SetComponent(i,a,v(i,a));
      arr_x0->SetComponent(i,a,x0(i,a));
      }
    mesh->GetPoints()->SetPoint(i,x_out);
    }

  // Write the mesh
  char buffer[2048];
  sprintf(buffer, filePattern.c_str(), k);
  WriteVTKData(mesh, buffer);
}

template <class TPixel, unsigned int VDim>
typename PointSetGeodesicToWarp<TPixel, VDim>::VectorImagePointer
PointSetGeodesicToWarp<TPixel, VDim>
::brute_force_method(
    const WarpGenerationParameters &param,
    HSystem &hsys, int tdir, int tStart, int tEnd)
{
  // Read the reference image
  ImagePointer imRef = LDDMM::img_read(param.fnReference.c_str());

  // Read the optional mask image. If not supplied, create a mask of all ones
  // to simplify the code below (fewer ifs)
  ImagePointer imMask;
  if(param.fnMask.size())
    {
    imMask = LDDMM::img_read(param.fnMask.c_str());
    itkAssertOrThrowMacro(imMask->GetBufferedRegion() == imRef->GetBufferedRegion(),
                          "Region mismatch between mask and reference");
    }
  else
    {
    imMask = LDDMM::new_img(imRef);
    imMask->FillBuffer(1.0);
    }

  // Count the number of masked pixels
  unsigned int nv_masked = 0;
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IterType;
  for(IterType it(imMask, imMask->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    nv_masked += it.Value() >= 1.0 ? 1 : 0;

  // Create an array to store the RAS coordinates of the points
  vnl_matrix<TPixel> X(nv_masked, VDim);

  // Initialize all the points inside the mask
  unsigned int ix = 0;
  itk::Point<double, VDim> pt_lps;
  for(IterType it(imMask, imMask->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(it.Value() >= 1.0)
      {
      // Compute the RAS coordinate of the current vertex position (applying displacement phi)
      imMask->TransformIndexToPhysicalPoint(it.GetIndex(), pt_lps);
      for(unsigned int a = 0; a < VDim; a++)
        X(ix,a) = a < 2 ? -pt_lps[a] : pt_lps[a];
      ++ix;
      }
    }

  // F for the Gaussian
  double gaussian_f = -1.0 / (2 * param.sigma * param.sigma);

  // Cutoff where the Gaussian is < 1e-6
  double d2_cutoff = 27.63102 * param.sigma * param.sigma;

  // Delta-t
  double dt = tdir * 1.0 / (param.N - 1);

  // Perform iteration over time
  for(int t = tStart; t != tEnd; t+=tdir)
    {
    // Get all landmarks
    const Matrix &qt = hsys.GetQt(t);
    const Matrix &pt = hsys.GetPt(t);

    // Multithread the vertices
    itk::MultiThreaderBase::Pointer mt = itk::MultiThreaderBase::New();
    itk::ImageRegion<1> full_region({{0}}, {{nv_masked}});
    mt->ParallelizeImageRegion<1>(
          full_region,
          [&qt, &pt, &X, dt, d2_cutoff, gaussian_f](const itk::ImageRegion<1> &thread_region)
      {
        int i0 = thread_region.GetIndex(0);
        int i1 = i0 + thread_region.GetSize(0);
        for(int i = i0; i < i1; i++)
          {
          // Velocity at this index
          TPixel xi[VDim], vi[VDim];
          for(unsigned int a = 0; a < VDim; a++)
            {
            xi[a] = X(i,a);
            vi[a] = 0.0;
            }

          // Loop over all q nodes
          for(unsigned int j = 0; j < qt.rows(); j++)
            {
            // Compute distance to that point
            double delta, d2 = 0;
            for(unsigned int a = 0; a < VDim; a++)
              {
              delta = xi[a] - qt(j,a);
              d2 += delta * delta;
              }

            // Only proceed if distance is below cutoff
            if(d2 < d2_cutoff)
              {
              // Take the exponent
              double g = exp(gaussian_f * d2);

              // Scale momentum by exponent
              for(unsigned int a = 0; a < VDim; a++)
                vi[a] += g * pt(j,a);
              }
            }

          // Use Euler's method to get position at next timepoint, encoding it as a LPS
          // displacement
          for(unsigned int a = 0; a < VDim; a++)
            X(i,a) += vi[a] * dt;
          }
      }, nullptr);

    std::cout << "." << std::flush;
    }

  // Form the phi image
  VectorImagePointer imPhi = LDDMM::new_vimg(imRef);
  ix = 0;
  IterType it(imMask, imMask->GetBufferedRegion());
  itk::ImageRegionIteratorWithIndex<VectorImageType> itPhi(imPhi, imMask->GetBufferedRegion());
  for(; !it.IsAtEnd(); ++it,++itPhi)
    {
    if(it.Value() >= 1.0)
      {
      // Compute the RAS coordinate of the current vertex position (applying displacement phi)
      imMask->TransformIndexToPhysicalPoint(it.GetIndex(), pt_lps);
      for(unsigned int a = 0; a < VDim; a++)
        {
        double x_lps = a < 2 ? -X(ix,a) : X(ix,a);
        double phi_lps = x_lps - pt_lps[a];
        itPhi.Value()[a] = phi_lps;
        }
      ++ix;
      }
    }

  return imPhi;
}

template <class TPixel, unsigned int VDim>
int
PointSetGeodesicToWarp<TPixel, VDim>
::run(const WarpGenerationParameters &param)
{
  // Read the VTK mesh containing the points
  vtkPolyData *mesh = ReadVTKData(param.fnMesh);
  if(!mesh)
    {
    cerr << "Failed to read mesh from " << param.fnMesh << endl;
    return -1;
    }

  // Read the momentum field
  vtkDataArray *arr_p0 = mesh->GetPointData()->GetArray("InitialMomentum");
  if(!arr_p0 || arr_p0->GetNumberOfComponents() != VDim)
    {
    cerr << "Failed to read initial momentum from " << param.fnMesh << endl;
    return -1;
    }

  // Count the number of non-null entries
  vector<unsigned int> index;
  for(unsigned int i = 0; i < arr_p0->GetNumberOfTuples(); i++)
    {
    bool has_value = 1;
    for(unsigned int a = 0; a < VDim; a++)
      {
      if(std::isnan(arr_p0->GetComponent(i,a)))
        {
        has_value = 0;
        break;
        }
      }
    if(has_value)
      index.push_back(i);
    }

  // Populate the q0 and p0 arrays
  unsigned int k = index.size();
  Matrix q0(k, VDim), p0(k,VDim), q1(k, VDim), p1(k,VDim);

  for(unsigned int i = 0; i < k; i++)
    {
    for(unsigned int a = 0; a < VDim; a++)
      {
      q0(i,a) = mesh->GetPoint(index[i])[a];
      p0(i,a) = arr_p0->GetComponent(i,a);
      }
    }

  // We have read the mesh successfully
  // Create the hamiltonian system
  HSystem hsys(q0, param.sigma, param.N, 0, param.n_threads);

  // Flow without gradients - we have streamlines
  hsys.FlowHamiltonian(p0, q1, p1);

  // Direction of interpolation
  int tdir = 1;
  int tStart = (tdir > 0) ? 0 : param.N - 1;
  int tEnd = (tdir > 0) ? param.N : 0;

  // F for the Gaussian
  double gaussian_f = -1.0 / (2 * param.sigma * param.sigma);

  // Cutoff where the Gaussian is < 1e-6
  double d2_cutoff = 27.63102 * param.sigma * param.sigma;

  // Delta-t
  double dt = tdir * 1.0 / (param.N - 1);

  // Apply the warp to other meshes
  std::list<WarpGenerationParameters::MeshPair>::const_iterator it;
  for(it = param.fnWarpMeshes.begin(); it != param.fnWarpMeshes.end(); ++it)
    {
    // Read the VTK mesh
    vtkPolyData *mesh_to_warp = ReadVTKData(it->first);
    if(!mesh_to_warp)
      {
      cerr << "Failed to read mesh from " << param.fnMesh << endl;
      return -1;
      }

    // Create the arrays to track the point positions
    Matrix m_x(mesh_to_warp->GetNumberOfPoints(), VDim);

    // Initialize
    for(int i = 0; i < m_x.rows(); i++)
      for(int a = 0; a < VDim; a++)
        m_x(i,a) = mesh_to_warp->GetPoint(i)[a];

    // Store the initial array
    Matrix m_x0 = m_x;

    // Some print-out
    std::cout << "Warping mesh " << it->first << " : " << std::flush;

    // Iterate over time
    for(int t = tStart; t != tEnd; t+=tdir)
      {
      // Get all landmarks
      const Matrix &qt = hsys.GetQt(t);
      const Matrix &pt = hsys.GetPt(t);

      // Create a vector for the velocity
      Matrix v(m_x.rows(), VDim);

      // Compute the velocity for each landmark
      for(int i = 0; i < m_x.rows(); i++)
        {
        // Local arrays for speed
        double xi[VDim];
        for(int a = 0; a < VDim; a++) 
          {
          xi[a] = m_x(i,a);
          v(i,a) = 0;
          }

        // Loop over all q nodes
        for(int j = 0; j < qt.rows(); j++)
          {
          // Compute distance to that point
          double delta, d2 = 0;
          for(int a = 0; a < VDim; a++) 
            {
            delta = xi[a] - qt(j,a);
            d2 += delta * delta;
            }

          // Only proceed if distance is below cutoff
          if(d2 < d2_cutoff)
            {
            // Take the exponent
            double g = exp(gaussian_f * d2);

            // Scale momentum by exponent
            for(int a = 0; a < VDim; a++) 
              v(i,a) += g * pt(j,a);
            }
          }

        // Use Euler's method to get position at next timepoint
        for(int a = 0; a < VDim; a++) 
          m_x(i,a) += v(i,a) * dt;
        }

      // If this is an animation frame, save it
      if((param.anim_freq > 0 && 0 == (t + 1) % param.anim_freq) || t + 1 == tEnd)
        UpdateAndWriteMesh(mesh_to_warp, m_x, v, m_x0, it->second, t + 1);

      std::cout << "." << std::flush;
      }

    // Done for this mesh - just have to save it
    std::cout << std::endl;
    }

  // Compute warp if that was requested
  if(param.fnReference.size() && param.fnOutWarp.size())
    {
    if(param.brute)
      {
      VectorImagePointer imPhi = brute_force_method(param, hsys, tdir, tStart, tEnd);
      LDDMM::vimg_write(imPhi, param.fnOutWarp.c_str());
      }
    else
      {
      // Read the reference image
      ImagePointer imRef;
      LDDMM::img_read(param.fnReference.c_str(), imRef);

      // If mask specified, read the mask
      ImagePointer imMask;
      if(param.fnMask.size())
        imMask = LDDMM::img_read(param.fnMask.c_str());

      // Image for splatting
      VectorImagePointer
          imSplat = LDDMM::new_vimg(imRef),
          imVelocity = LDDMM::new_vimg(imRef),
          imLagragean = LDDMM::new_vimg(imRef),
          imPhi = LDDMM::new_vimg(imRef);

      // Compute the scaling factor for the Gaussian smoothing
      double scaling = 1.0;
      for(unsigned int a = 0; a < VDim; a++)
        {
        scaling *= sqrt(2 * vnl_math::pi) * param.sigma / imRef->GetSpacing()[a];
        }

      // Iterate over time
      // for(unsigned int t = 0; t < param.N; t++)
      for(int t = tStart; t != tEnd; t+=tdir)
        {
        // Get all landmarks
        const Matrix &qt = hsys.GetQt(t);
        const Matrix &pt = hsys.GetPt(t);

        // Splatting approach
        // Splat each landmark onto the splat image
        imSplat->FillBuffer(typename VectorImageType::PixelType(0.0));

        // Create an interpolator for splatting
        FastInterpolator flint(imSplat);

        for(int i = 0; i < k; i++)
          {
          // Map landmark point to continuous index
          ContIndexType cix;
          PointType point;
          typename VectorImageType::InternalPixelType vec;
          for(unsigned int a = 0; a < VDim; a++)
            {
            // Here we have to correct for the fact that the landmark coordinates are
            // assumed to be in RAS space, but point is in LPS space. We therefore have
            // to apply the LPS/RAS transform before splatting
            if(a < 2)
              {
              point[a] = -qt(i,a);
              vec[a] = -pt(i,a);
              }
            else
              {
              point[a] = qt(i,a);
              vec[a] = pt(i,a);
              }

            }
          imRef->TransformPhysicalPointToContinuousIndex(point, cix);

          // Splat landmark's momentum at the point location
          flint.Splat(cix.GetVnlVector().data_block(), &vec);
          }

        // Now all the momenta have been splatted. The next step is to smooth with a Gaussian
        typename LDDMM::Vec vec_sigma; vec_sigma.Fill(param.sigma);
        LDDMM::vimg_smooth(imSplat, imVelocity, vec_sigma);

        // Accumulate the velocity into phi - using imSplat as a temporary
        LDDMM::vimg_scale_in_place(imVelocity, scaling);

        // Euler interpolation
        LDDMM::interp_vimg(imVelocity, imPhi, 1.0, imSplat, false, true);
        LDDMM::vimg_add_scaled_in_place(imPhi, imSplat, tdir * 1.0 / (param.N - 1));

        cout << "." << flush;

        // Save the warp if needed
        if(t + 1 == tEnd || (param.anim_freq && 0 == (t + 1) % param.anim_freq))
          {
          char buffer[2048];
          sprintf(buffer, param.fnOutWarp.c_str(), t + 1);

          // Save the warp
          LDDMM::vimg_write(imPhi, buffer);
          }
        }
      }
    }

  return 0;
}



int main(int argc, char *argv[])
{
  WarpGenerationParameters param;

  if(argc < 2)
    return usage();

  // Process parameters
  for(int i = 1; i < argc; i++)
    {
    string arg = argv[i];
    if(arg == "-r" && i < argc-1)
      {
      param.fnReference = argv[++i];
      }
    else if(arg == "-m" && i < argc-1)
      {
      param.fnMesh = argv[++i];
      }
    else if(arg == "-o" && i < argc-1)
      {
      param.fnOutWarp = argv[++i];
      }
    else if(arg == "-g" && i < argc-1)
      {
      param.fnMask = argv[++i];
      }
    else if(arg == "-s" && i < argc-1)
      {
      param.sigma = atof(argv[++i]);
      }
    else if(arg == "-d" && i < argc-1)
      {
      param.dim = (unsigned int) atoi(argv[++i]);
      }
    else if(arg == "-n" && i < argc-1)
      {
      param.N = (unsigned int) atoi(argv[++i]);
      }
    else if(arg == "-a" && i < argc-1)
      {
      param.anim_freq = (unsigned int) atoi(argv[++i]);
      }
    else if(arg == "-B")
      {
      param.brute = true;
      }
    else if(arg == "-t")
      {
      param.n_threads = (unsigned int) atoi(argv[++i]);
      }
    else if(arg == "-M" && i < argc-2)
      {
      string fn1 = std::string(argv[++i]);
      string fn2 = std::string(argv[++i]);
      WarpGenerationParameters::MeshPair pair = std::make_pair(fn1,fn2);
      param.fnWarpMeshes.push_back(pair);
      }
    else if(arg == "-h")
      {
      return usage();
      }
    else
      {
      cerr << "Unknown option " << arg << endl;
      return -1;
      }
    }

  check(param.sigma > 0, "Missing or negative sigma parameter");
  check(param.N > 0 && param.N < 10000, "Incorrect N parameter");
  check(param.dim >= 2 && param.dim <= 3, "Incorrect N parameter");
  check(param.fnMesh.length(), "Missing target filename");

  if(param.fnWarpMeshes.size() == 0)
    {
    check(param.fnReference.length(), "Missing template filename");
    check(param.fnOutWarp.length(), "Missing output filename");
    }

  // Set the number of threads if not specified
  if(param.n_threads == 0)
    param.n_threads = std::thread::hardware_concurrency();
  else
    itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(param.n_threads);

  // Specialize by dimension
  if(param.dim == 2)
    return PointSetGeodesicToWarp<float,2>::run(param);
  else
    return PointSetGeodesicToWarp<float,3>::run(param);


  return 0;
}
