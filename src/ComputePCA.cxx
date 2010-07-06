#include "ScriptInterface.h"
#include <iostream>

using namespace std;

int usage()
{
  cout << 
    "cmrep_pca_generate: generate PCA from a set of models \n"
    "usage:\n"
    "  cmrep_pca_generate [options] -o outfile sample1.cmrep ... sampleN.cmrep\n"
    "options:\n"
    "  -o FILE        : Output filename for the PCA (no extension)\n"
    "  -r FILE.cmrep  : Specify cm-rep model to use as a reference\n";
  return -1;
}

extern void SetupSignalHandlers();

int main(int argc, char *argv[])
{
  int i = 0;

  // Report errors
  SetupSignalHandlers();

  // Check the number of input parameters
  if(argc < 7)
    {
    cerr << "Run cmrep_pca_generate -h for usage" << endl;
    return -1;
    }

  // Required parameters
  string fn_ref, fn_out;

  // Parse command line
  for(i = 1; i < argc; i++)
    {
    string arg = argv[i];
    if(arg == "-r")
      {
      fn_ref = argv[++i];
      }
    else if(arg == "-o")
      {
      fn_out = argv[++i];
      }
    else if(arg == "-h" || arg == "--help")
      {
      usage();
      return 0;
      }
    else if(arg[0] == '-')
      {
      cerr << "Unknown option " << arg << endl;
      return -1;
      }
    else
      {
      break;
      }
    }

  // Check required parameters
  if(!fn_out.length())
    {
    cerr << "Missing output filename" << endl;
    return -1;
    }
  
  // Check optional parameters
  if(!fn_ref.length())
    {
    fn_ref = argv[i];
    cout << "Using first model, " << fn_ref << " as reference model" << endl;
    }

  // Create Medial PCA object
  MedialPCA mpca;
  std::vector<MedialPDE *> pde;
  MedialPDE *ref = NULL;

  // Read the filenames and models
  for(; i < argc; i++)
    {
    MedialPDE *model = NULL;
    try 
      {
      MedialPDE *model = new MedialPDE(argv[i]);
      mpca.AddSample(model);
      pde.push_back(model);
      }
    catch(std::exception &exc)
      {
      cerr << "Exception reading sample " << argv[i] << endl;
      cerr << exc.what() << endl;
      return -1;
      }
    }

  // Read the reference model
  try 
    {
    ref = new MedialPDE(fn_ref.c_str());
    }
  catch(std::exception &exc)
    {
    cerr << "Exception reading reference model " << fn_ref << endl;
    cerr << exc.what() << endl;
    return -1;
    }

  // Compute PCA
  cout << "Starting PCA computation..." << endl;
  mpca.ComputePCA(ref);

  // Export the PCA matrix
  string outmat = fn_out + "_shapemat.mat";
  mpca.ExportShapeMatrix(outmat.c_str());
}
