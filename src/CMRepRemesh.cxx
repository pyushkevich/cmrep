#include "ScriptInterface.h"
#include <string>
#include <iostream>

using namespace std;

int usage()
{
  cout << "cmrep_remesh : remesh medial surface in subdivision surface cm-rep" << endl;
  cout << "Usage:" << endl;
  cout << "  cmrep_remesh [options] input.cmrep output.cmrep " << endl;
  cout << "Options: " << endl;
  cout << "  none so far" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 3) return usage();

  try 
    {
    // Load the cm-rep (should be subdivision type)
    SubdivisionMPDE mpde(argv[1]);

    // Get the subdivision surface
    mpde.Remesh();

    // Save as a new mesh
    mpde.SaveToParameterFile(argv[2]);
    }
  catch(std::exception &exc)
    {
    cerr << "Exception: " << exc.what() << endl;
    return -1;
    }
  catch(string &error)
    {
    cerr << "Exception: " << error << endl;
    return -1;
    }
  catch(...)
    {
    cerr << "Unknown exception" << endl;
    return -1;
    }
  return 0;
}
