#include "ScriptInterface.h"
#include <iostream>
#include <string>

using namespace std;
using namespace medialpde;

int usage()
{
  cout << "usage: mpdeprof input.mpde input.img" << endl;
  return -1;
}
  

int main(int argc, char *argv[])
{
  if(argc < 3) return usage();

  // Load the image
  FloatImage img;
  img.LoadFromFile(argv[2]);

  // Load the model
  MedialPDE mp(4, 6, 32, 80, 0.5, 2, 2);
  mp.LoadFromParameterFile(argv[1]);
  mp.SetNumberOfCoefficients(4, 6);
  mp.SetOptimizerToConjugateGradientDescent(0.1);
  mp.SetMatchToVolumeOverlap();
  mp.SetOptimizationToCoarseToFine(4, 6, 4, 6);
  mp.RunOptimization(&img, 6);
}
