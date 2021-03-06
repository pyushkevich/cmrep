#ifndef __CoarseToFineMappingDescriptor_h_
#define __CoarseToFineMappingDescriptor_h_

#include "Registry.h"
#include "OptimizationParameters.h"
#include <vnl/vnl_vector.h>
#include "MedialException.h"

class GenericMedialModel;

class CoarseToFineMappingDescriptor {
public:
  virtual ~CoarseToFineMappingDescriptor() {}
  typedef vnl_vector<size_t> MaskVector;

  /** 
   * Get the number of levels that a model supports. Levels are a fairly
   * abstract concept and will vary substantially depending on the model
   * representation. For subdivision surfaces, each level of the subdivision is
   * a coarse-to-fine level. For Fourier surfaces, levels are generated by
   * subtracting one from the number of coefficients in both u and v.
   */
  virtual size_t GetNumberOfLevels() const = 0;

  /** 
   * Get the mask corresponding to a given level. Levels are indexed as arrays
   * in Mathematica: 1 .. N or -1 to -N. The coarsest level is 1 (or -N) and the
   * finest level is -1 (or N).
   */
  virtual MaskVector GetMaskForLevel(int level) const = 0;

  /**
   * Get the mask specified by a registry folder. This will be very specific to
   * the individual models.
   */
  virtual MaskVector GetMask(const CoarseToFineSettings *settings) const = 0;

};

/**
 * This is the c2f mapping descriptor for a Fourier (cosine) surface. The levels
 * here are pretty simple. The surface is defined by coefficients C[i,j]. If
 * there are M x N coefficients and M > N, the levels, from finest to coarsest are
 *
 * M x N, (M-1) x (N-1), ... (M - N + 2) x 2
 */
class FourierCoarseToFineMappingDescriptor : public CoarseToFineMappingDescriptor
{
public:
  virtual ~FourierCoarseToFineMappingDescriptor() {}

  FourierCoarseToFineMappingDescriptor(size_t m, size_t n) 
    {
    this->m = m;
    this->n = n;
    this->nlevels = std::min(n, m) - 1;
    }

  size_t GetNumberOfLevels() const
    { return nlevels; }
    
  MaskVector GetMaskForLevel(int level) const
    {
    assert(-((int)nlevels) <= level && level >= (int) nlevels && level != 0);
    if(level > 0)
      {
      size_t d = (size_t) (nlevels - level);
      return this->GetMask(m-d, n-d, m-d, n-d);
      }
    else
      {
      size_t d = (size_t) (1 - level);
      return GetMask(m-d, n-d, m-d, n-d);
      }
    }

  MaskVector GetMask(const CoarseToFineSettings *ctfSettings) const
    {
    size_t nxu = 0, nxv = 0, nru = 0, nrv = 0;

    // Cast the settings to the right format
    try {
      const FourierCoarseToFineSettings *mySettings = 
        static_cast<const FourierCoarseToFineSettings *>(ctfSettings);

      // Read the registry, check for overruns
      nxu = mySettings->nxu;
      nxv = mySettings->nxv;
      nru = mySettings->nru;
      nrv = mySettings->nrv;
    }
    catch(...) {
      throw ModelIOException(
        "Coarse-to-fine settings don't match Fourier-based medial model");
    }

    // Check that the numbers match the model
    if(nxu > m || nru > m || nxv > n || nrv > n)
      throw ModelIOException(
        "Coarse-to-fine settings exceed number of parameters "
        "in the Fourier model");

    // One more check (do we need this?)
    if(nxu * nxv + nru * nrv == 0)
      throw ModelIOException(
        "There are no parameters selected in the coarse-to-fine settings");

    // Return the corresponding mask
    return this->GetMask(nxu, nxv, nru, nrv);
    }

  /** 
   * The complete mask specification is the number of coefficients in X and the
   * number of coefficients in Rho.
   */
  MaskVector GetMask(size_t nxu, size_t nxv, size_t nru, size_t nrv) const 
    {
    MaskVector mask(m * n * 4, 0);

    size_t off = 0;
    for(size_t j = 0; j < n; j++) for(size_t i = 0; i < m; i++) 
      {
      if(i < nxu && j < nxv)
        mask[off] = mask[off+1] = mask[off+2] = 1;
      if(i < nru && j < nrv)
        mask[off+3] = 1;
      off+=4;
      }

    return mask;
    }

private:
  size_t m, n;
  size_t nlevels;

};

class SubdivisionSurfaceCoarseToFineMappingDescriptor : public CoarseToFineMappingDescriptor
{
public:
  virtual ~SubdivisionSurfaceCoarseToFineMappingDescriptor() {}
  typedef vnl_vector<size_t> MaskVector;

  SubdivisionSurfaceCoarseToFineMappingDescriptor(size_t n)
    { this->n = n; }

  /** 
   * Get the number of levels that a model supports. There is only one level.
   * However, you can subdivide the model 'by hand' to get more coefficients
   */
  virtual size_t GetNumberOfLevels() const { return 1; }

  /** 
   * Get the mask corresponding to a given level. There is only one level where 
   * all the control points are selected
   */
  virtual MaskVector GetMaskForLevel(int level) const 
    { return MaskVector(n * 4, 1); }

  /**
   * Get the mask specified by a registry folder. This method always returns a
   * mask of all ones
   */
  virtual MaskVector GetMask(const CoarseToFineSettings *settings) const 
    { return MaskVector(n * 4, 1); }

private:
  size_t n;

};

#endif // __CoarseToFineMappingDescriptor_h_
