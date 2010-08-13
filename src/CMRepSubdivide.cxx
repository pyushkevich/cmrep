#include "ScriptInterface.h"
#include <iostream>
#include <cstdlib>

using namespace std;

int usage()
{
  cout << "cmrep_subdivide: apply subdivision to cm-rep model" << endl;
  cout << "usage:" << endl;
  cout << "  cmrep_subdivide [options] input.cmrep output.cmrep n_ctl n_atom" << endl;
  cout << "parameters:" << endl;
  cout << "  n_ctl  : number of subdivision levels applied to control pts." << endl;
  cout << "  n_atom : number of levels applied to the sampling grid" << endl;
  cout << "options:" << endl;
  cout << "  -b     : apply subdivision only to boundary triangles" << endl;
  return -1;
}

int main(int argc, char *argv[])
  {
  if(argc < 5)
    return usage();

  bool adaptive = false;
  for(size_t iarg = 1; iarg < argc - 4; iarg++)
    {
    if(!strcmp(argv[iarg], "-b"))
      adaptive = true;
    else
      { 
      cerr << "Unrecognized option " << argv[iarg] << endl;
      return -1;
      }
    }
  
  try 
    {
    // Load a cmrep from file  
    SubdivisionMPDE cmrep(argv[argc - 4]);

    // Apply subdivision
    size_t nc = atoi(argv[argc - 2]), na = atoi(argv[argc - 1]);
    if(adaptive)
      cout << "Subdividing boundary triangles\n "
        << "coeffs by " << nc << " levels; atoms by " << na << " levels" << endl;
    else
      cout << "Subdividing coeffs by " << nc << " levels; atoms by " << na << " levels" << endl;
    cmrep.SubdivideMeshes(nc, na, adaptive);

    // Save the cm-rep to a new file
    cmrep.SaveToParameterFile(argv[argc - 3]);
    }
  catch(std::exception &exc)
    {
    cerr << "Exception: " << exc.what() << endl;
    return -1;
    }
  return 0;
  }

