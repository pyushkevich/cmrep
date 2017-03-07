#include "Procrustes.h"
#include "ReadWriteVTK.h"

#include <vtkSmartPointer.h>
#include <iostream>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_svd.h>

#include <vtkGenericDataObjectReader.h>

using namespace std;

int usage()
{
  cout << "vtkprocrustes: compute procrustes alignment between two meshes" << endl;
  cout << "Usage:" << endl;
  cout << "  vtkprocrustes [options] fixed.vtk moving.vtk output.mat" << endl;
  cout << "Options: " << endl;
  cout << "  -r                 : do a rigid (no scaling) match between meshes" << endl;
  cout << "  -a                 : do a full affine match between meshes" << endl;
  cout << "  -f                 : flip meshes before procrustes" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  bool flip = false;
  bool affine = false;
  bool rigid = false;
  if(argc < 4)
    return usage();

  std::string fnOut = argv[--argc];
  std::string fnMoving = argv[--argc];
  std::string fnFixed = argv[--argc];
  for(int i = 1; i < argc; i++)
    {
    std::string arg = argv[i];
    if(arg == "-a")
      {
      affine = true;
      }
    else if(arg == "-r")
      {
      rigid = true;
      }
    else if(arg == "-f")
      {
      flip = true;
      }
    else
      {
      cerr << "Unknown argument " << arg << endl;
      return usage();
      }
    }

  if(affine && rigid)
    {
    cerr << "Affine and rigid cannot be used together" << endl;
    return -1;
    }

  // Read the two meshes
  vtkSmartPointer<vtkGenericDataObjectReader> rFixed = vtkGenericDataObjectReader::New();
  rFixed->SetFileName(fnFixed.c_str());
  rFixed->Update();
  vtkSmartPointer<vtkPointSet> mFixed = dynamic_cast<vtkPointSet *>(rFixed->GetOutput());
  
  vtkSmartPointer<vtkGenericDataObjectReader> rMoving = vtkGenericDataObjectReader::New();
  rMoving->SetFileName(fnMoving.c_str());
  rMoving->Update();
  vtkSmartPointer<vtkPointSet> mMoving = dynamic_cast<vtkPointSet *>(rMoving->GetOutput());

  if(!mFixed)
    {
    cerr << "Unable to read point set for fixed mesh" << endl;
    return -1;
    }

  if(!mMoving)
    {
    cerr << "Unable to read point set for moving mesh" << endl;
    return -1;
    }

  // Confirm points match
  int k = mFixed->GetNumberOfPoints();
  if(k != mMoving->GetNumberOfPoints())
    {
    cerr << "Number of points does not match between meshes" << endl;
    return -1;
    }

  // Read meshes into a pair of matrices
  vnl_matrix<double> Xfix(k,3), Xmov(k,3), Xfix2mov, A(4,4);
  A.set_identity();

  // Apply flip in x if specified - arbitrary
  for(int i = 0; i < k; i++)
    {
    double x[3];
    mFixed->GetPoint(i, x);
    Xfix(i,0) = (flip) ? -x[0] : x[0];
    Xfix(i,1) = x[1];
    Xfix(i,2) = x[2];

    mMoving->GetPoint(i, x);
    Xmov(i,0) = x[0];
    Xmov(i,1) = x[1];
    Xmov(i,2) = x[2];
    }

  // Set up the results of computation
  if(affine)
    {
    // In case of affine match, this is a 'simple' least square problem
    
    // Add columns of ones
    vnl_matrix<double> X(k, 4); X.fill(1.0); X.update(Xfix, 0, 0);
    vnl_matrix<double> Y(k, 4); Y.fill(1.0); Y.update(Xmov, 0, 0);

    // Construct the LSQ matrix
    vnl_matrix<double> M(12,12); M.fill(0.0);
    M.update(X.transpose() * X, 0, 0);
    M.update(X.transpose() * X, 4, 4);
    M.update(X.transpose() * X, 8, 8);

    // Construct the right hand side
    vnl_matrix<double> YX = Y.transpose() * X;
    vnl_vector<double> b(YX.data_block(), 12);

    // Solve system Mp = b
    vnl_svd<double> svd(M); 
    vnl_vector<double> z = svd.solve(b);

    // Copy parameters into the A matrix
    int iz = 0;
    for(int i = 0; i < 3; i++)
      {
      for(int j = 0; j < 4; j++)
        {
        A(i,j) = z(iz++);
        }
      }
    }
  else
    {
    vnl_matrix<double> R;
    vnl_vector<double> b;
    double s;
    ExtendedProcrustesAnalysis(Xfix, Xmov, Xfix2mov, R, b, s);

    if(rigid)
      s = 1.0;

    A.update(R.transpose() * s, 0, 0);
    A(0,3) = b[0];
    A(1,3) = b[1];
    A(2,3) = b[2];
    }


  // Update the matrix with the flip
  vnl_matrix<double> F(4,4); 
  F.set_identity();
  if(flip)
    F(0,0) = -1.0;

  A = A * F;

  // Save the matrix A to output file
  ofstream matrixFile;
  matrixFile.open(fnOut.c_str());
  matrixFile << A;
  matrixFile.close();
}

