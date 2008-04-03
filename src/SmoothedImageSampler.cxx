#include "SmoothedImageSampler.h"
#include <iostream>

SmoothedImageSampler
::SmoothedImageSampler(
  double sigma,           // Standard deviation of gaussian
  double alpha,           // The number of sigmas at which to cut off the Gaussian
  float *img,             // Pointer to raw image values
  double *bb_start,       // Start of the bounding box
  double *bb_end)
{
  // Store the parameters
  this->sigma = sigma;
  this->alpha = alpha;
  this->img = img;
  for(unsigned int d = 0; d < 3; d++)
    {
    this->bb_start[d] = bb_start[d];
    this->bb_end[d] = bb_end[d];
    }

  // Number of interpixel ticks in the bounding box
  ntx = (int)(bb_end[0] - bb_start[0] + 0.5);
  nty = (int)(bb_end[1] - bb_start[1] + 0.5);
  ntz = (int)(bb_end[2] - bb_start[2] + 0.5);

  // Arrays to store the erf values in x, y and z
  dx = new double[ntx]; dy = new double[nty]; dz = new double[ntz];

  // Derivative arrays to store the erf values in x, y and z
  ddx = new double[ntx]; ddy = new double[nty]; ddz = new double[ntz];

  // Set up the stides
  stride_z = nty * ntx; 
  stride_y = ntx;

  // Get the sigma scale factors
  sfx = 1.0 / (sqrt(2.0) * sigma);
  sfy = 1.0 / (sqrt(2.0) * sigma);
  sfz = 1.0 / (sqrt(2.0) * sigma);

  // Compute the cutoffs alpha * sigma
  cutx = sigma * alpha;
  cuty = sigma * alpha;
  cutz = sigma * alpha;
}

SmoothedImageSampler
::~SmoothedImageSampler()
{
  delete dx; delete dy; delete dz;
  delete ddx; delete ddy; delete ddz;
}

double  
SmoothedImageSampler
::Sample(const double *X, double *grad_f)
{
  // The bound variables for x, y, z
  int ix0, ix1, iy0, iy1, iz0, iz1;

  // Compute the ERF difference arrays
  bool inside = true;
  inside = inside & 
    compute_erf_array_with_deriv(dx, ddx, ix0, ix1, bb_start[0], ntx, cutx, X[0], sfx);
  inside = inside & 
    compute_erf_array_with_deriv(dy, ddy, iy0, iy1, bb_start[1], nty, cuty, X[1], sfy);
  inside = inside & 
    compute_erf_array_with_deriv(dz, ddz, iz0, iz1, bb_start[2], ntz, cutz, X[2], sfz);

  // If we ain't inside, return 0
  if(!inside)
    {
    grad_f[0] = 0.0;
    grad_f[1] = 0.0;
    grad_f[2] = 0.0;
    return 0.0;
    }


  // printf("Range is [%d %d] [%d %d] [%d %d]\n", ix0, ix1, iy0, iy1, iz0, iz1);

  // Get a pointer to the output value
  double sum_wf = 0.0, sum_w = 0.0;
  double sum_wfx = 0.0, sum_wx = 0.0;
  double sum_wfy = 0.0, sum_wy = 0.0;
  double sum_wfz = 0.0, sum_wz = 0.0;


  // Loop over the voxels in the region identified
  for(int iz = iz0; iz < iz1; iz++)
    {
    double wz = dz[iz];
    double dwdz = ddz[iz];
    int oz = stride_z * iz;

    for(int iy = iy0; iy < iy1; iy++)
      {
      double wyz = wz * dy[iy];
      double dwdy_wz = wz * ddy[iy];
      double wy_dwdz = dwdz * dy[iy];
      int oyz = oz + stride_y * iy;

      for(int ix = ix0; ix < ix1; ix++)
        {
        int oxyz = oyz + ix;
        double w = wyz * dx[ix];
        double wx_dwdy_wz = dwdy_wz * dx[ix];
        double wx_wy_dwdz = wy_dwdz * dx[ix];
        double dwdx_wy_wz = wyz * ddx[ix];

        // Sample the image at this location
        double ival = img[oxyz];

        // Accumulate the function
        // printf("Pixel %10d (%4d %4d %4d) with value %8f weighted %7.3f\n", oxyz, ix, iy, iz, ival, w);
        sum_wf += w * ival;
        sum_w += w;

        // Compute the derivatives 
        sum_wfx += dwdx_wy_wz * ival;
        sum_wx += dwdx_wy_wz;
        
        // Compute the derivatives 
        sum_wfy += wx_dwdy_wz * ival;
        sum_wy += wx_dwdy_wz;

        // Compute the derivatives 
        sum_wfz += wx_wy_dwdz * ival;
        sum_wz += wx_wy_dwdz;
        }
      }
    }

  // Scaling factor for speed
  double inv_sum_w = 1.0 / sum_w;

  // Set the derivative values
  double f = sum_wf * inv_sum_w;
  grad_f[0] = (sum_wfx - f * sum_wx) * inv_sum_w;
  grad_f[1] = (sum_wfy - f * sum_wy) * inv_sum_w;
  grad_f[2] = (sum_wfz - f * sum_wz) * inv_sum_w;

  // f = dy[48];
  // grad_f[1] = ddy[49];
  // grad_f[0] = grad_f[2] = 0;

  // Set the output value
  return f;
}


double 
SmoothedImageSampler
::Sample(const double *X)
{
  // The bound variables for x, y, z
  int ix0, ix1, iy0, iy1, iz0, iz1;

  // Compute the ERF difference arrays
  bool inside = true;
  inside = inside & 
    compute_erf_array_with_deriv(dx, ddx, ix0, ix1, bb_start[0], ntx, cutx, X[0], sfx);
  inside = inside & 
    compute_erf_array_with_deriv(dy, ddy, iy0, iy1, bb_start[1], nty, cuty, X[1], sfy);
  inside = inside & 
    compute_erf_array_with_deriv(dz, ddz, iz0, iz1, bb_start[2], ntz, cutz, X[2], sfz);

  // If we ain't inside, return 0
  if(!inside) 
    {
    return 0;
    }

  // printf("Range is [%d %d] [%d %d] [%d %d]\n", ix0, ix1, iy0, iy1, iz0, iz1);

  // Get a pointer to the output value
  double sum_wf = 0.0, sum_w = 0.0;

  // Loop over the voxels in the region identified
  for(int iz = iz0; iz < iz1; iz++)
    {
    double wz = dz[iz];
    int oz = stride_z * iz;

    for(int iy = iy0; iy < iy1; iy++)
      {
      double wyz = wz * dy[iy];
      int oyz = oz + stride_y * iy;

      for(int ix = ix0; ix < ix1; ix++)
        {
        int oxyz = oyz + ix;
        double w = wyz * dx[ix];

        // Sample the image at this location
        double ival = img[oxyz];

        // Accumulate the function
        // printf("Pixel %10d (%4d %4d %4d) with value %8f weighted %7.3f\n", oxyz, ix, iy, iz, ival, w);
        sum_wf += w * ival;
        sum_w += w;
        }
      }
    }

  // Scaling factor for speed
  double inv_sum_w = 1.0 / sum_w;

  // Set the output value
  return sum_wf / sum_w;
}


bool SmoothedImageSampler::compute_erf_array(
  double *dx_erf,         // The output array of erf(p+i+1) - erf(p+i)
  int &k0, int &k1,       // The range of integration 0 <= k0 < k1 <= n
  double b,               // Lower bound of the bounding box
  int n,                  // Size of the bounding box in steps
  double cut,             // The distance at which to cut off
  double p,               // the value p
  double sfac)            // scaling factor 1 / (Sqrt[2] sigma)
{
  // x is the offset of the point wrt the bounding box
  double x = p - b;

  // These variables define the range over which to integrate
  double u0 = x - cut, u1 = x + cut;

  // Squeeze these variables to the allowed range
  if(u0 < 0.0) u0 = 0.0; if(u1 > n) u1 = n;

  // Determine the range of voxels over which to integrate
  k0 = (int) floor(u0); k1 = (int) ceil(u1);

  // Squeeze the range so that it is valid
  if(k0 >= k1) return false;

  // The weight of each voxel is the integral of exp(t^2 / 2) from
  // the beginning of the voxel to the end of the voxel. But for the
  // first and last voxel, the integral is from u0 to the end of the
  // voxel and from the end of the voxel to u1, correspondingly.

  // Start at the first voxel, at u0
  double e_last = erf((u0 - x) * sfac), e_now;

  // Iterate over the middle voxels
  double t = (k0 - x) * sfac;
  for(int i = k0; i < k1-1; i++)
    {
    t += sfac;
    e_now = erf(t);
    dx_erf[i] = e_now - e_last;
    e_last = e_now;
    }

  // Handle the last voxel
  e_now = erf((u1 - x) * sfac);
  dx_erf[k1-1] = e_now - e_last;

  return true;
}


bool SmoothedImageSampler::compute_erf_array_with_deriv(
  double *dx_erf,         // The output array of erf(p+i+1) - erf(p+i)
  double *dx_erf_deriv,   // The derivative of the output array wrt x
  int &k0, int &k1,       // The range of integration 0 <= k0 < k1 <= n
  double b,               // Lower bound of the bounding box
  int n,                  // Size of the bounding box in steps
  double cut,             // The distance at which to cut off
  double p,               // the value p
  double sfac)            // scaling factor 1 / (Sqrt[2] sigma)
{
  // x is the offset of the point wrt the bounding box
  double x = p - b;

  // Scaling factor for derivative computations
  double dscale = - sfac * 1.128379167;

  // These variables define the range over which to integrate
  double u0 = x - cut, u1 = x + cut;

  // Squeeze these variables to the allowed range
  if(u0 < 0.0) u0 = 0.0; if(u1 > n) u1 = n;

  // Determine the range of voxels over which to integrate
  k0 = (int) floor(u0); k1 = (int) ceil(u1);

  // Squeeze the range so that it is valid
  if(k0 >= k1) return false;

  // The weight of each voxel is the integral of exp(t^2 / 2) from
  // the beginning of the voxel to the end of the voxel. But for the
  // first and last voxel, the integral is from u0 to the end of the
  // voxel and from the end of the voxel to u1, correspondingly.

  // Start at the first voxel, at u0
  double t0 = (u0 - x) * sfac;
  double e_last = erf(t0), e_now;
  double d_last = dscale * exp(-t0 * t0), d_now;
  d_last = 0.0;

  // Iterate over the middle voxels
  double t = (k0 - x) * sfac;
  for(int i = k0; i < k1-1; i++)
    {
    t += sfac;

    e_now = erf(t);
    d_now = dscale * exp(-t * t);

    dx_erf[i] = e_now - e_last;
    dx_erf_deriv[i] = d_now - d_last;

    e_last = e_now;
    d_last = d_now;
    }

  // Handle the last voxel
  double t1 = (u1 - x) * sfac;
  e_now = erf(t1);
  d_now = dscale * exp(-t1 * t1);
  d_now = 0.0;

  dx_erf[k1-1] = e_now - e_last;
  dx_erf_deriv[k1-1] = d_now - d_last;
  
  // dx_erf[k0] = erf(t0);
  // dx_erf_deriv[k0] = dscale * exp(-t0 * t0);

  return true;
}
