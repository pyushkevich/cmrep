#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkTriangleFilter.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_math.h>
#include <vnl/algo/vnl_qr.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyDataWriter.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>

using namespace std;

int usage()
{
  cout << "cmrep_afftran: " << endl;
  cout << "  compute affine transform between two cm-rep models" << endl;
  cout << "usage:" << endl;
  cout << "  cmrep_to_warp [opts] cmrep_reference cmrep_moving out_transform" << endl;
  cout << "parameters:" << endl;
  cout << "  cmrep_reference, " << endl;
  cout << "  cmrep_moving    : VTK meshes (cmr_fit mesh dir, defX.med.vtk) " << endl;
  cout << "  out_transform   : Where to write the transform to" << endl;
  cout << "options:" << endl;
  cout << "  -x Xmax Xstep   : xi coord. sampling (def: 1 0.1)" << endl;
  cout << "                    specifies whether the warp field extends past " << endl;
  cout << "                    the model's boundary (i.e., 2 0.1). The step " << endl;
  cout << "                    size affects the size of tetrahedra used to " << endl;
  cout << "                    the warp field" << endl;
  cout << "  -f              : flip atom normal, i.e., map spoke1 in cmrep1 to " << endl;
  cout << "                    spoke2 in cmrep2" << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;

  if(argc < 4) return usage();

  // Get the main parameters
  string fn_out       = argv[argc-1];
  string fn_cmrep_mov = argv[argc-2];
  string fn_cmrep_ref = argv[argc-3];

  // Get the options
  double xiMax = 1, xiStep = 0.1;
  bool flip = false;

  // Parse the options
  for(size_t iopt = 1; iopt < argc-3; iopt++)
    {
    if(!strcmp(argv[iopt], "-f"))
      {
      flip = true;
      }
    else if(!strcmp(argv[iopt], "-x"))
      {
      xiMax = atof(argv[++iopt]);
      xiStep = atof(argv[++iopt]);
      }
    else
      {
      cerr << "Unknown option " << argv[iopt] << endl;
      return usage();
      }
    }

  // Load the two cm-rep meshes
  vtkPolyDataReader *r1 = vtkPolyDataReader::New();
  r1->SetFileName(fn_cmrep_ref.c_str());
  r1->Update();
  vtkPolyData *mesh_ref = r1->GetOutput();

  vtkPolyDataReader *r2 = vtkPolyDataReader::New();
  r2->SetFileName(fn_cmrep_mov.c_str());
  r2->Update();
  vtkPolyData *mesh_mov = r2->GetOutput();

  // Get the spoke arrays, radius
  vtkDataArray *spoke_ref[2], *spoke_mov[2];
  spoke_ref[0] = mesh_ref->GetPointData()->GetArray("Spoke1");
  spoke_ref[1] = mesh_ref->GetPointData()->GetArray("Spoke2");
  
  if(!flip)
    {
    spoke_mov[0] = mesh_mov->GetPointData()->GetArray("Spoke1");
    spoke_mov[1] = mesh_mov->GetPointData()->GetArray("Spoke2");
    }
  else
    {
    spoke_mov[0] = mesh_mov->GetPointData()->GetArray("Spoke2");
    spoke_mov[1] = mesh_mov->GetPointData()->GetArray("Spoke1");
    }

  if(spoke_ref[0] == NULL || spoke_ref[1] == NULL)
    {
    cerr << "Spoke arrays missing in reference cmrep" << endl;
    return -1;
    }
  if(spoke_mov[0] == NULL || spoke_mov[1] == NULL)
    {
    cerr << "Spoke arrays missing in moving cm-rep" << endl;
    return -1;
    }
  if(mesh_ref->GetNumberOfPoints() != mesh_mov->GetNumberOfPoints())
    {
    cerr << "CM-rep meshes have different number of points " << endl;
    return -1;
    }

  // Compute the xi array 
  std::vector<double> xivec; 
  xivec.push_back(0.0);
  for(double xi = xiStep; xi < xiMax + 0.5 * xiStep; xi += xiStep)
    { 
    xivec.insert(xivec.end(), xi);
    xivec.insert(xivec.begin(), -xi);
    }

  // Compute least square fit of points in m1 to m2
  size_t np = mesh_ref->GetNumberOfPoints();
  size_t ns = xivec.size(); 
  size_t n = np * ns;
  size_t off = 0;
  Mat X_ref(n, 4), X_mov(n, 4);
  for(size_t i = 0; i < np; i++)
    {
    typedef vnl_vector_fixed<double, 3> Vec;
    Vec m_ref  = Vec(mesh_ref->GetPoint(i));
    Vec s1_ref = Vec(spoke_ref[0]->GetTuple(i));
    Vec s2_ref = Vec(spoke_ref[1]->GetTuple(i));

    Vec m_mov  = Vec(mesh_mov->GetPoint(i));
    Vec s1_mov = Vec(spoke_mov[0]->GetTuple(i));
    Vec s2_mov = Vec(spoke_mov[1]->GetTuple(i));

    for(size_t j = 0; j < ns; j++)
      {
      Vec p_ref = m_ref + ((xivec[j] < 0) ? - xivec[j] * s1_ref : xivec[j] * s2_ref);
      Vec p_mov = m_mov + ((xivec[j] < 0) ? - xivec[j] * s1_mov : xivec[j] * s2_mov);
      for(size_t d = 0; d < 3; d++)
        {
        X_ref(off, d) = p_ref(d);
        X_mov(off, d) = p_mov(d);
        }
      X_ref(off, 3) = 1.0;
      X_mov(off, 3) = 1.0;
      off++;
      }
    }

  Mat M(12,12, 0.0);
  Mat XmXm = X_mov.transpose() * X_mov;
  M.update(XmXm, 0, 0);
  M.update(XmXm, 4, 4);
  M.update(XmXm, 8, 8);
  Mat XrXm = X_ref.transpose() * X_mov;
  Vec b(XrXm.data_block(), 12);
  vnl_qr<double> qr(M);
  Vec z = qr.solve(b);
  Mat tran(4, 4);
  tran.set_identity();
  tran.update(Mat(z.data_block(), 3, 4));

  // Write the transformation matrix
  ofstream fout(fn_out.c_str());
  for(size_t i = 0; i < 4; i++)
    fout << tran[i][0] << " " << tran[i][1] << " " << tran[i][2] << " " << tran[i][3] << endl;
  fout.close();

  /* 
  // Apply transform to the mesh
  vtkMatrixToLinearTransform *m_Transform = vtkMatrixToLinearTransform::New();
  m_Transform->GetMatrix()->DeepCopy(tran.data_block());  
  vtkTransformPolyDataFilter *m_TransformFilter = vtkTransformPolyDataFilter::New();
  m_TransformFilter->SetTransform(m_Transform);
  m_TransformFilter->SetInput(mesh_mov);
  m_TransformFilter->Update();

  // Create a writer
  vtkPolyDataWriter *m_Writer = vtkPolyDataWriter::New();
  m_Writer->SetFileName("dummy.vtk");
  m_Writer->SetInput(m_TransformFilter->GetOutput());
  m_Writer->Update();
  */

  /*

  // Write out a transform
  typedef itk::AffineTransform<double, 3> ATranType;
  ATranType::Pointer atran = ATranType::New();
  itk::Matrix<double, 3, 3> mat;
  mat.GetVnlMatrix().update(tran.extract(3,3));
  atran->SetMatrix(mat);
  itk::Vector<double, 3> vec;
  vec.GetVnlVector().update(tran.get_column(3).extract(3));
  atran->SetOffset(vec);

  // Compute the mean squared difference
  double diff = 0.0;
  for(size_t i = 0; i < n; i++)
    {
    vnl_vector<double> p_ref = X_ref.get_row(i).extract(3);
    vnl_vector<double> p_mov = X_mov.get_row(i).extract(3);
    itk::Point<double, 3> pt_mov; pt_mov.GetVnlVector().update(p_mov);
    vnl_vector<double> p_fit = atran->TransformPoint(pt_mov).GetVnlVector();
    diff += (p_ref - p_fit).squared_magnitude();
    }
  diff /= n;
  cout << "Mean Square Difference: " << diff << endl;

  // Invert the transform, because what we really want is the transform
  // from reference space to moving space
  ATranType::Pointer ainv = ATranType::New();
  atran->GetInverse(ainv);

  // Write the inverse transform
  typedef itk::TransformFileWriter WriterType;
  WriterType::Pointer w = WriterType::New();
  w->SetInput(ainv);
  w->SetFileName(fn_out.c_str());
  w->Update();

  */
}

