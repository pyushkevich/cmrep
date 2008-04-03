#include "vnl/vnl_random.h"
#include "MedialAtom.h"
#include "MedialModelIO.h"
#include "vtkOBJReader.h"
#include "vtkBYUWriter.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "SubdivisionSurface.h"
#include "PDESubdivisionMedialModel.h"
#include "SubdivisionSurfaceMedialIterationContext.h"
#include "ScriptInterface.h"
#include "BranchingSubdivisionSurface.h"

#include <string>
#include <iostream>

void ExportMedialMeshToVTK(GenericMedialModel *model, const char *file);
void ExportBoundaryMeshToVTK(GenericMedialModel *model, const char *file);

using namespace std;

/**
 * Test subdivision surface functionality
 */
int TestSubdivisionSurface(const char *objMesh)
{
  // Read the input mesh
  vtkOBJReader *reader = vtkOBJReader::New();
  reader->SetFileName(objMesh);
  reader->Update();
  vtkPolyData *poly = reader->GetOutput();

  // Create a subdivision surface
  SubdivisionSurface::MeshLevel mesh;
  SubdivisionSurface::ImportLevelFromVTK(poly, mesh);

  // Describe the mesh
  cout << "Input mesh has " << mesh.triangles.size() << " triangles";
  cout << " and " << mesh.nVertices << " vertices" << endl;

  // Save the input surface for comparison
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(poly);
  writer->SetFileName("sub_input.vtk");
  writer->Update();

  // Check the mesh after loading
  if(!SubdivisionSurface::CheckMeshLevel(mesh))
    return 1;

  // Subdivide the mesh once
  SubdivisionSurface::MeshLevel meshsub;
  SubdivisionSurface::Subdivide(&mesh, &meshsub);

  cout << "Subdivided mesh has " << meshsub.triangles.size() << " triangles";
  cout << " and " << meshsub.nVertices << " vertices" << endl;

  // Check the subdivided mesh
  if(!SubdivisionSurface::CheckMeshLevel(meshsub))
    return 1;

  // Compute the subdivision surface
  vtkPolyData *polysub = vtkPolyData::New();
  SubdivisionSurface::ApplySubdivision(poly, polysub, meshsub);

  // Save the subdivision surface
  writer->SetInput(polysub);
  writer->SetFileName("sub_subdivide.vtk");
  writer->Update();

  // Subdivide the mesh one more time
  SubdivisionSurface::MeshLevel meshresub;
  SubdivisionSurface::Subdivide(&meshsub, &meshresub);

  cout << "Subdivided mesh has " << meshresub.triangles.size() << " triangles";
  cout << " and " << meshresub.nVertices << " vertices" << endl;

  // Check the subdivided mesh
  if(!SubdivisionSurface::CheckMeshLevel(meshresub))
    return 1;

  // Save the subdivision surface
  SubdivisionSurface::ApplySubdivision(poly, polysub, meshresub);
  writer->SetInput(polysub);
  writer->SetFileName("sub_resubdivide.vtk");
  writer->Update();



  return 0;
}

/**
 * Test branching subdivision surface
 */
int TestBranchingSubdivisionSurface(const char *objMesh)
{
  // Read the mesh using one of the VTK readers
  vtkPolyData *poly = ReadVTKMesh(objMesh);

  // Create a subdivision surface
  BranchingSubdivisionSurface::MeshLevel mesh;
  BranchingSubdivisionSurface::ImportLevelFromVTK(poly, mesh);

  // Describe the mesh
  cout << "Input mesh has " << mesh.triangles.size() << " triangles";
  cout << " and " << mesh.vertices.size() << " vertices" << endl;

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(poly);
  writer->SetFileName("branch_input.vtk");
  writer->Update();

  // Check the mesh after loading
  if(!BranchingSubdivisionSurface::CheckMeshLevel(mesh))
    return 1;

  // Subdivide the mesh once
  BranchingSubdivisionSurface::MeshLevel meshsub;
  BranchingSubdivisionSurface::Subdivide(&mesh, &meshsub);

  cout << "Subdivided mesh has " << meshsub.triangles.size() << " triangles";
  cout << " and " << meshsub.vertices.size() << " vertices" << endl;

  // Check the subdivided mesh
  if(!BranchingSubdivisionSurface::CheckMeshLevel(meshsub))
    return 1;

  // Compute the subdivision surface
  vtkPolyData *polysub = vtkPolyData::New();
  BranchingSubdivisionSurface::ApplySubdivision(poly, polysub, meshsub);

  // Save the subdivision surface
  writer->SetInput(polysub);
  writer->SetFileName("branch_subdivide.vtk");
  writer->Update();

  // Subdivide the mesh one more time
  BranchingSubdivisionSurface::MeshLevel meshresub;
  BranchingSubdivisionSurface::Subdivide(&meshsub, &meshresub);

  cout << "Subdivided mesh has " << meshresub.triangles.size() << " triangles";
  cout << " and " << meshresub.vertices.size() << " vertices" << endl;

  // Check the subdivided mesh
  if(!BranchingSubdivisionSurface::CheckMeshLevel(meshresub))
    return 1;

  // Save the subdivision surface
  BranchingSubdivisionSurface::ApplySubdivision(poly, polysub, meshresub);
  writer->SetInput(polysub);
  writer->SetFileName("branch_resubdivide.vtk");
  writer->Update();

  return 0;
}

/**
 * Test sparse matrix multiplication
 */
int TestSparseMatrix()
{
  int rc = 0;

  // Create a sparse matrix with random elements
  size_t i, j, k, N1 = 5, N2 = 4, N3 = 6;
  typedef vnl_sparse_matrix<int> Mutable;
  typedef ImmutableSparseMatrix<int> Immutable;

  // Create two mutable matrices
  Mutable A(N1, N2), B(N2, N3), C(N1, N3);

  // Initialize them with some random values (-10 to 10)
  vnl_random rnd;
  for(i = 0; i < N1; i++) for(j = 0; j < N2; j++)
    if(rnd.lrand32(0, 4) == 0)
      A(i,j) = rnd.lrand32(0,18) - 9;

  for(i = 0; i < N2; i++) for(j = 0; j < N3; j++)
    if(rnd.lrand32(0, 4) == 0)
      B(i,j) = rnd.lrand32(0,18) - 9;

  // Compute their product using VNL
  A.mult(B, C);

  // Compute the product using immutable matrices
  Immutable AI, BI, CI, DI;
  AI.SetFromVNL(A);
  BI.SetFromVNL(B);
  CI.SetFromVNL(C);
  Immutable::Multiply(DI, AI, BI);

  // Print the results
  cout << "A is " << AI << endl;
  cout << "B is " << BI << endl;
  cout << "C is " << CI << endl;
  cout << "D is " << DI << endl;

  // Compare CI and DI
  int rcMult = (CI == DI) ? 0 : 1;
  if(rcMult)
    { cout << "Sparse Matrix Multiplication Error" << endl; }

  // Create a matrix for sparse ATA calculation
  size_t M1 = 6, M2 = 7;
  Mutable AA_sparse(M1, M2);
  vnl_matrix<int> AA(M1, M2);
  for(i = 0; i < M1; i++) for(j = 0; j < M2; j++)
    if(rnd.lrand32(0, 4) == 0)
      AA_sparse(i,j) = AA(i,j) = rnd.lrand32(0,18) - 9;

  // Try the computation for ATA
  Immutable AAA;
  AAA.SetFromVNL(AA_sparse);

  Immutable ATA1;
  Immutable::InitializeATA(ATA1, AAA);
  Immutable::ComputeATA(ATA1, AAA);

  // Compute ATA directly
  vnl_matrix<int> ATA2 = AA.transpose() * AA;

  // Print the results
  cout << "ATA1 is " << endl << ATA1 << endl;
  cout << "ATA2 is " << endl << ATA2 << endl;

  // Check for correctness
  bool flagATACorrect = true;
  for(size_t i = 0; i < ATA1.GetNumberOfRows(); i++)
    {
    for(Immutable::RowIterator r = ATA1.Row(i); !r.IsAtEnd(); ++r)
      {
      if(ATA2(i, r.Column()) != r.Value())
        flagATACorrect = false;
      }
    }

  int rcATA = (flagATACorrect) ? 0 : 1;
  if(rcATA)
    { cout << "Sparse Matrix ATA Error" << endl; }

  return rcMult + rcATA;
}

int TestSubdivisionPDE(const char *objMesh)
{
  // Read the input mesh
  vtkPolyData *poly = ReadVTKMesh(objMesh);

  // Create a subdivision surface
  SubdivisionSurface::MeshLevel mesh;
  SubdivisionSurface::ImportLevelFromVTK(poly, mesh);

  // Get the polygon data
  vnl_vector<double> xCtl(mesh.nVertices * 4);

  // Set the u and v values of the data
  vnl_vector<double> u(mesh.nVertices, 0.0), v(mesh.nVertices, 0.0);

  // Get the texture arrays
  vtkDataArray *uv = poly->GetPointData()->GetTCoords();

  // Read the coordinates and u,v from the mesh
  for(size_t i = 0; i < mesh.nVertices; i++)
    {
    // Set the coordinate value
    for(size_t d = 0; d < 3; d++)
      xCtl[i * 4 + d] = poly->GetPoint(i)[d];

    // Set the rho value
    xCtl[i * 4 + 3] = -0.25;

    // Set U and V
    if(uv)
      {
      u[i] = uv->GetTuple(i)[0];
      v[i] = uv->GetTuple(i)[1];
      }
    else
      {
      u[i] = xCtl[i * 4];
      v[i] = xCtl[i * 4 + 1];
      }
    }

  // Create the medial model
  PDESubdivisionMedialModel model;

  // Pass the mesh to the medial model with specified number of levels
  model.SetMesh(mesh, xCtl, u, v, 2, 0);
  model.ComputeAtoms();

  // Temp: test derivatives
  int iTest = model.GetSolver()->TestPartialDerivatives();

  // Save the result
  ExportMedialMeshToVTK(&model, "mesh_medial.vtk");
  ExportBoundaryMeshToVTK(&model, "mesh_boundary.vtk");

  // Save the model itself
  MedialModelIO::WriteModel(&model, "subdivision_model.cmrep");

  return iTest;
}

int TestModelSubdivision(const char *file)
{
  medialpde::SubdivisionMPDE mp1(file), mp2(file);
  mp2.SubdivideMeshes(1,0);

  GenericMedialModel *mm1 = mp1.GetMedialModel();
  GenericMedialModel *mm2 = mp2.GetMedialModel();

  double maxerr = 0.0;
  for(size_t i = 0; i < mm1->GetNumberOfAtoms(); i++)
    {
    double del = (mm1->GetAtomArray()[i].X - mm2->GetAtomArray()[i].X).magnitude();
    if(del > maxerr) maxerr = del;
    }

  return (maxerr > 1.0e-10);
}


int usage()
{
  cout << "testsub: MedialPDE Test Module" << endl;
  cout << "  usage: testpde TEST_ID [parameters] " << endl;
  cout << "  tests: " << endl;
  cout << "    SUBSURF1 XX.obj            Test subdivision with OBJ mesh" << endl;
  cout << "    SUBSURF2 XX.obj            Test subdivision PDE" << endl;
  cout << "    BRANCH1 XX.obj             Test branching subdivision surface" << endl;
  cout << "    SPARSEMAT                  Test sparse matrix routines" << endl;
  cout << "    MODELSUB                   Test model subdivistion" << endl;
  cout << endl;
  return -1;
}

int main(int argc, char *argv[])
{
  // Different tests that can be executed
  if(argc == 1) return usage();

  // Choose a test depending on the parameters
  if(0 == strcmp(argv[1], "SUBSURF1") && argc > 2)
    return TestSubdivisionSurface(argv[2]);
  if(0 == strcmp(argv[1], "BRANCH1") && argc > 2)
    return TestBranchingSubdivisionSurface(argv[2]);
  else if(0 == strcmp(argv[1], "SUBSURF2") && argc > 2)
    return TestSubdivisionPDE(argv[2]);
  else if(0 == strcmp(argv[1], "SPARSEMAT"))
    return TestSparseMatrix();
  else if(0 == strcmp(argv[1], "MODELSUB"))
    return TestModelSubdivision(argv[2]);
  else
    return usage();
}

