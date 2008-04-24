#include "ScriptInterface.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << "cmrep_subdivide: apply subdivision to cm-rep model" << endl;
  cout << "usage:" << endl;
  cout << "  cmrep_subdivide input.cmrep output.cmrep n_ctl n_atom" << endl;
  cout << "parameters:" << endl;
  cout << "  n_ctl  : number of subdivision levels applied to control pts." << endl;
  cout << "  n_atom : number of levels applied to the sampling grid" << endl;
  return -1;
}

int main(int argc, char *argv[])
  {
  if(argc != 5)
    return usage();
  
  try 
    {
    // Load a cmrep from file  
    SubdivisionMPDE cmrep(argv[1]);

    // Apply subdivision
    size_t nc = atoi(argv[3]), na = atoi(argv[4]);
    cout << "Subdividing coeffs by " << nc << " levels; atoms by " << na << " levels" << endl;
    cmrep.SubdivideMeshes(nc, na);

    // Save the cm-rep to a new file
    cmrep.SaveToParameterFile(argv[2]);
    }
  catch(std::exception &exc)
    {
    cerr << "Exception: " << exc.what() << endl;
    return -1;
    }
  return 0;
  }

