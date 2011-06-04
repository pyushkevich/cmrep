#include <iostream>
#include <string>
#include <sstream>

#include "ScriptInterface.h"

using namespace std;

int main(int argc, char *argv[])
{
  cout << "Usage: cmrep_meanmodel model1.cmrep ... modelN.cmrep output.cmrep" << endl;
  cout << "  Computes the mean of n models " << endl;
  if(argc < 4) return -1;

  // Set up PCA
  MedialPCA mpca;

  // Read all models
  std::vector<MedialPDE *> mpde;
  MedialPDE *target;
  for(int i = 1; i < argc-1;i++)
    mpca.AddSample(target = new MedialPDE(argv[i]));

  // Compute PCA
  mpca.ComputePCA(target);
  mpca.SetFSLocationToMean();
  mpca.GetShapeAtFSLocation(target);
  target->SaveToParameterFile(argv[argc-1]);
}
