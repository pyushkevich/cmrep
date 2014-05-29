#ifndef __OptimizationParameters_h_
#define __OptimizationParameters_h_

#include "Registry.h"
#include <smlmath.h>

class CoarseToFineSettings;

class OptimizationParameters {
public:
  enum Optimizer{ CONJGRAD, GRADIENT, EVOLUTION };
  enum Mapping { 
    AFFINE, COARSE_TO_FINE, IDENTITY, PCA, 
    RADIUS_SUBSET, POSITION_SUBSET, REFLECTION, LAPLACE_BASIS }; 
  enum ImageMatch { VOLUME, BOUNDARY, PROBABILITY_INTEGRAL, RADIUS_VALUES };
  enum PenaltyTerm { 
    BOUNDARY_JACOBIAN = 0, MEDIAL_REGULARITY, ATOM_BADNESS, 
    BOUNDARY_GRAD_R, RADIUS, MEDIAL_CURVATURE, BND_CURVATURE,
    MEDIAL_ANGLES, LOOP_VALIDITY, BOUNDARY_ANGLES, CROSS_CORRELATION, LOCAL_DISTANCE, DIFFEOMORPHIC, 
    BND_JACOBIAN_DISTORTION, MED_JACOBIAN_DISTORTION, 
    BND_ELASTICITY, CLOSEST_POINT, N_PENALTY_TERMS };
  enum CTFSettings { COSINE_BASIS_PDE, LOOP_SUBDIVISION_PDE, NONE };

  /** Optimizer type */
  Optimizer xOptimizer;

  /** Image match type */
  ImageMatch xImageMatch;

  /** Mapping Type */
  Mapping xMapping;

  /** Coarse-to-fine settings */
  CoarseToFineSettings *xCTFSettings;

  /** PCA file name (defined if PCA exists) */
  std::string xPCAFileName;
  
  /** Number of PCA modes to use in PCA-based optimization */
  size_t nPCAModes;

  /** Number of eigenfunctions to use in coefficient mapping */
  size_t xLaplaceBasisSize;

  /** Reflection plane information (only defined for reflection mapping) */
  SMLVec3d xReflectionPlane;
  double xReflectionIntercept;

  /** Weights of the various penalty terms */
  std::map<PenaltyTerm, double> xTermWeights, xTermDefaultWeights;

  /** Parameters associated with the weights */
  std::map<PenaltyTerm, Registry> xTermParameters;

  /** Constructor */
  OptimizationParameters();

  /** Destructor */
  ~OptimizationParameters();

  /** Read the parameters from registry */
  void ReadFromRegistry(Registry &reg);

  /** Write the parameters to y */
  void WriteToRegistry(Registry &reg);

private:

  // Enum-to-string maps
  RegistryEnumMap<Optimizer> xOptimizerRegMap;
  RegistryEnumMap<Mapping> xMappingRegMap;
  RegistryEnumMap<ImageMatch> xImageMatchRegMap;
  RegistryEnumMap<PenaltyTerm> xPenaltyTermRegMap;
  RegistryEnumMap<CTFSettings> xCTFSettingsRegMap;
};


class CoarseToFineSettings {
public:
  virtual void ReadFromRegistry(Registry &folder) = 0;
  virtual void WriteToRegistry(Registry &folder) = 0;
  virtual OptimizationParameters::CTFSettings GetType() const = 0;
  virtual ~CoarseToFineSettings() {}
};

class FourierCoarseToFineSettings : public CoarseToFineSettings {
public:
  size_t nxu, nxv, nru, nrv;

  void ReadFromRegistry(Registry &folder)
    {
    nxu = folder["SpatialU"][0];
    nxv = folder["SpatialV"][0];
    nru = folder["RadialU"][0];
    nrv = folder["RadialV"][0];
    }

  void WriteToRegistry(Registry &folder)
    {
    folder["SpatialU"] << nxu;
    folder["SpatialV"] << nxv;
    folder["RadialU"] << nru;
    folder["RadialV"] << nrv;
    }

  OptimizationParameters::CTFSettings GetType() const
    { return OptimizationParameters::COSINE_BASIS_PDE; }

  virtual ~FourierCoarseToFineSettings() {}
};


#endif // __OptimizationParameters_h_
