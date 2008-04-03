#include "ScriptInterface.h"
#include "BasisFunctions2D.h"
#include "MedialAtom.h"
#include "CartesianMedialModel.h"
#include "OptimizationTerms.h"
#include "CoefficientMapping.h"
#include "MedialAtomGrid.h"
#include "PrincipalComponents.h"
#include "TestSolver.h"

using namespace medialpde;

int usage()
{
  cout << "cmr_fitbin: fit cm-rep to binary image" << endl;
  cout << "usage cmrep_binary template.cmrep image.img result.cmrep param.reg nSteps" << endl;
}

int main(int argc, char *argv[])
{
  if(argc != 6)
    return usage();

  string fnimg = argv[2];
  size_t i = fnimg.rfind('.');
  string fnbase = fnimg.substr(0,i);
  string fnext = fnimg.substr(i+1);

  MedialPDE cmrep(argv[1]);
  FloatImage jet;
  jet.LoadFromPath(fnbase.c_str(), fnext.c_str());
  jet.SetOutsideValue(-1.0);
  
  cmrep.RunOptimization(&jet, atoi(argv[5]), argv[4]);
  cmrep.SaveToParameterFile(argv[3]);
}
