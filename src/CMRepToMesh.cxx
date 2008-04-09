#include "ScriptInterface.h"
#include <string>
#include <iostream>

using namespace std;

using namespace medialpde;

int usage()
{
  cout << "cmrep2mesh : create medial and boundary meshes from cm-rep" << endl;
  cout << "Usage:" << endl;
  cout << "  cmrep2mesh model.cmrep outname " << endl;
  cout << "Output: " << endl;
  cout << "  creates VTK mesh files outname.med.vtk and outname.bnd.vtk" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  if(argc < 3) return usage();

  try 
    {
    // Load the cm-rep
    MedialPDE mpde(argv[1]);

    // Create mesh filenames
    string fnmed = string(argv[2]) + ".med.vtk";
    string fnbnd = string(argv[2]) + ".bnd.vtk";

    // Generate a mesh
    mpde.SaveVTKMesh(fnmed.c_str(), fnbnd.c_str());

    // Write what we did
    cout << "Wrote meshes " << fnmed << " and " << fnbnd << endl;
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
