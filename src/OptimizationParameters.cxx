#include "OptimizationParameters.h"

OptimizationParameters
::OptimizationParameters()
{
  // Initialize the enum mappings if necessary
  xOptimizerRegMap.AddPair(CONJGRAD, "ConjGrad");
  xOptimizerRegMap.AddPair(GRADIENT, "Gradient");
  xOptimizerRegMap.AddPair(EVOLUTION, "EvolStrat");

  xMappingRegMap.AddPair(AFFINE, "Affine");
  xMappingRegMap.AddPair(COARSE_TO_FINE, "CoarseFine");
  xMappingRegMap.AddPair(IDENTITY, "Identity");
  xMappingRegMap.AddPair(PCA, "PCA");
  xMappingRegMap.AddPair(RADIUS_SUBSET, "RadiusSubset");
  xMappingRegMap.AddPair(POSITION_SUBSET, "PositionSubset");
  xMappingRegMap.AddPair(REFLECTION, "Reflection");
  xMappingRegMap.AddPair(LAPLACE_BASIS, "LaplaceBasis");

  xImageMatchRegMap.AddPair(VOLUME, "VolumeOverlap");
  xImageMatchRegMap.AddPair(PROBABILITY_INTEGRAL, "ProbabilityIntegral");
  xImageMatchRegMap.AddPair(BOUNDARY, "BoundaryIntegral");
  xImageMatchRegMap.AddPair(RADIUS_VALUES, "CurrentRadius");

  xPenaltyTermRegMap.AddPair(BOUNDARY_JACOBIAN, "BoundaryJacobianEnergyTerm");
  xPenaltyTermRegMap.AddPair(BOUNDARY_GRAD_R, "BoundaryGradRPenaltyTerm");
  xPenaltyTermRegMap.AddPair(MEDIAL_REGULARITY, "MedialRegularityTerm");
  xPenaltyTermRegMap.AddPair(MEDIAL_CURVATURE, "MedialCurvaturePenaltyTerm");
  xPenaltyTermRegMap.AddPair(BND_CURVATURE, "BoundaryCurvaturePenaltyTerm");
  xPenaltyTermRegMap.AddPair(ATOM_BADNESS, "AtomBadnessTerm");
  xPenaltyTermRegMap.AddPair(RADIUS, "RadiusPenaltyTerm");
  xPenaltyTermRegMap.AddPair(MEDIAL_ANGLES, "MedialAnglesPenaltyTerm");
  xPenaltyTermRegMap.AddPair(LOOP_VALIDITY, "LoopTangentSchemeValidityPenaltyTerm");
  xPenaltyTermRegMap.AddPair(BOUNDARY_ANGLES, "BoundaryAnglesPenaltyTerm");
  xPenaltyTermRegMap.AddPair(CROSS_CORRELATION, "CrossCorrelation");
  xPenaltyTermRegMap.AddPair(LOCAL_DISTANCE, "LocalDistancePenaltyTerm");
  xPenaltyTermRegMap.AddPair(DIFFEOMORPHIC, "DiffeomorphicEnergyTerm");
  xPenaltyTermRegMap.AddPair(BND_JACOBIAN_DISTORTION, "BoundaryJacobianDistortionPenaltyTerm");
  xPenaltyTermRegMap.AddPair(MED_JACOBIAN_DISTORTION, "MedialJacobianDistortionPenaltyTerm");
  xPenaltyTermRegMap.AddPair(BND_ELASTICITY, "BoundaryElasticityPrior");
  xPenaltyTermRegMap.AddPair(CLOSEST_POINT, "SymmetricClosestPoint");

  xCTFSettingsRegMap.AddPair(COSINE_BASIS_PDE, "CosineBasisPDE");
  xCTFSettingsRegMap.AddPair(LOOP_SUBDIVISION_PDE, "LoopSubdivisionPDE");
  xCTFSettingsRegMap.AddPair(NONE, "NONE");

  // Set the default weights for the terms
  xTermDefaultWeights[BOUNDARY_JACOBIAN] = 0.5;
  xTermDefaultWeights[BOUNDARY_GRAD_R] = 0.0;
  xTermDefaultWeights[MEDIAL_REGULARITY] = 1.0;
  xTermDefaultWeights[MEDIAL_CURVATURE] = 0.0;
  xTermDefaultWeights[BND_CURVATURE] = 0.0;
  xTermDefaultWeights[ATOM_BADNESS] = 0.0;
  xTermDefaultWeights[RADIUS] = 0.1;
  xTermDefaultWeights[MEDIAL_ANGLES] = 0.0;
  xTermDefaultWeights[LOOP_VALIDITY] = 0.0;
  xTermDefaultWeights[BOUNDARY_ANGLES] = 0.0;
  xTermDefaultWeights[CROSS_CORRELATION] = 0.0;
  xTermDefaultWeights[LOCAL_DISTANCE] = 0.0;
  xTermDefaultWeights[DIFFEOMORPHIC] = 0.0;
  xTermDefaultWeights[BND_JACOBIAN_DISTORTION] = 0.0;
  xTermDefaultWeights[MED_JACOBIAN_DISTORTION] = 0.0;
  xTermDefaultWeights[BND_ELASTICITY] = 0.0;
  xTermDefaultWeights[CLOSEST_POINT] = 0.0;

  // Clear the settings pointer
  xCTFSettings = NULL;

  // Set default values
  xOptimizer = CONJGRAD;
  xImageMatch = VOLUME;
  xMapping = IDENTITY;

  xLaplaceBasisSize = 0;
}

OptimizationParameters
::~OptimizationParameters()
{
  if(xCTFSettings != NULL) delete xCTFSettings;
}

void
OptimizationParameters
::ReadFromRegistry(Registry &R)
{
  // Read the enums
  xOptimizer = R["Optimizer"].GetEnum(xOptimizerRegMap, CONJGRAD);
  xImageMatch = R["ImageMatch"].GetEnum(xImageMatchRegMap, VOLUME);
  xMapping = R["Mapping"].GetEnum(xMappingRegMap, IDENTITY);

  // Read the parameters of different methods
  for(size_t i = 0; i < N_PENALTY_TERMS; i++)
    {
    // The current penalty term
    PenaltyTerm pt = static_cast<PenaltyTerm>(BOUNDARY_JACOBIAN + i);

    // Check if the term is a registry folder
    if(R.IsFolder(xPenaltyTermRegMap(pt)))
      {
      // New way of specifying parameters (Folder.Weight for weights, 
      // Folder.ParamName for various parameters)
      xTermWeights[pt] =
        R[xPenaltyTermRegMap(pt)+".Weight"][xTermDefaultWeights[pt]];
      xTermParameters[pt] = R.Folder(xPenaltyTermRegMap(pt));
      }
    else
      {
      // Old way of specifying stuff, where there was only a weight given
      // for each penalty term
      xTermWeights[pt] = 
        R[xPenaltyTermRegMap(pt)][xTermDefaultWeights[pt]];
      }



    }

  // Read the coarse-to-fine specification
  if(xCTFSettings != NULL) delete xCTFSettings;
  CTFSettings c2f = R["CoarseToFine.Mode"].GetEnum(xCTFSettingsRegMap, NONE);
  if(c2f == COSINE_BASIS_PDE)
    xCTFSettings = new FourierCoarseToFineSettings();
  else xCTFSettings = NULL;

  // Read the coarse to fine settings if they are there
  if(xCTFSettings != NULL)
    xCTFSettings->ReadFromRegistry(R.Folder("CoarseToFine"));

  // Read the PCA settings if they are there
  xPCAFileName = R["PCA.FileName"][""];
  nPCAModes = R["PCA.NumberOfModes"][10];

  // Read reflection plane info
  if(xMapping == REFLECTION)
    {
    xReflectionPlane[0] = R["Reflection.Normal[0]"][1.0];
    xReflectionPlane[1] = R["Reflection.Normal[1]"][0.0];
    xReflectionPlane[2] = R["Reflection.Normal[2]"][0.0];
    xReflectionIntercept = R["Reflection.Intercept"][0.0];
    }

  else if(xMapping == LAPLACE_BASIS)
    {
    xLaplaceBasisSize = R["LaplaceBasis.Size"][10];
    }
}
