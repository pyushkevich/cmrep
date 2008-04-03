class ParameterFile {
public:
  enum Optimizer{ CONJGRAD, GRADIENT, EVOLUTION };
  enum Mapping { AFFINE, COARSE_TO_FINE, IDENTITY, PCA }; 
  enum ImageMatch { VOLUME, BOUNDARY };
