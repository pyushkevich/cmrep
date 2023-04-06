#include <iostream>
#include <string>
#include "ReadWriteVTK.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkPointData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkDataReader.h>

using namespace std;

// Templated mesh IO methods
template <class TMeshType>
vtkSmartPointer<TMeshType> ReadMesh(const char *fname)
{ return NULL; }

template <>
vtkSmartPointer<vtkUnstructuredGrid> ReadMesh<>(const char *fname)
{
  vtkNew<vtkUnstructuredGridReader> reader;
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}

template <>
vtkSmartPointer<vtkPolyData> ReadMesh<>(const char *fname)
{
  vtkNew<vtkPolyDataReader> reader;
  reader->SetFileName(fname);
  reader->Update();
  return reader->GetOutput();
}


template <class TMeshType>
void WriteMesh(TMeshType *mesh, const char *fname, bool vtk_binary)
{ }

template <>
void WriteMesh<>(vtkUnstructuredGrid *mesh, const char *fname, bool vtk_binary)
{
  vtkNew<vtkUnstructuredGridWriter> writer;
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(vtk_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

template <>
void WriteMesh<>(vtkPolyData *mesh, const char *fname, bool vtk_binary)
{
  vtkNew<vtkPolyDataWriter> writer;
  writer->SetFileName(fname);
  writer->SetInputData(mesh);
  if(vtk_binary)
    writer->SetFileTypeToBinary();
  writer->Update();
}

int usage()
{
  cout << "mesh_tetra_sample: sample tetrahedral mesh using a surface mesh" << endl;
  cout << "Usage: mesh_tetra_sample [options] mesh.vtk tetra.vtk output.vtk array_name" << endl;
  cout << "Required parameters:" << endl;
  cout << "   mesh.vtk       : Mesh at whose vertices you want to sample from the tetrahedra" << endl;
  cout << "   tetra.vtk      : Tetrahedral mesh with some cell array that you want to sample" << endl;
  cout << "   output.vtk     : Mesh that you want to save" << endl;
  cout << "   mesh.vtk       : Name of the cell array" << endl;
  cout << "Options: " << endl;
  cout << "   -d <value>     : If the vertex does not belong to any cell, assign it the value of the" << endl;
  cout << "                    closest cell that is less than <value> distance away" << endl;
  cout << "   -B             : Write VTK files as binary" << endl;
  cout << "   -b <value>     : Background value (when vertex falls outside of the image)" << endl;
  cout << "                    defaults to NaN" << endl;
  cout << "   -D <array>     : Save the distance to closest tetrahedron into an array with given name" << endl;
  return -1;
}

struct Parameters
{
  std::string fnMesh, fnTetra, fnOutput, array_name, dist_array_name;
  double search_dist = 0., fill_value = std::numeric_limits<double>::quiet_NaN();
  bool flagBinary = false;
};

template <class TMesh>
int MeshTetraSample(Parameters p)
{
  // Read the main mesh
  vtkSmartPointer<TMesh> mesh = ReadMesh<TMesh>(p.fnMesh.c_str());
  cout << "Read sampling mesh from " << p.fnMesh << " with " << mesh->GetNumberOfPoints() << " points." << endl;

  // Read the tetrahedral mesh
  vtkSmartPointer<vtkUnstructuredGrid> tetra = ReadMesh<vtkUnstructuredGrid>(p.fnTetra.c_str());
  cout << "Read tetrahedral mesh from " << p.fnTetra << " with " << mesh->GetNumberOfCells() << " cells." << endl;

  // Get the source array
  vtkSmartPointer<vtkDataArray> a_source = tetra->GetCellData()->GetArray(p.array_name.c_str());
  if(!a_source)
    {
    cerr << "Cell array " << p.array_name << " missing in tetrahedral mesh" << endl;
    return -1;
    }

  // Create the output array
  vtkNew<vtkDoubleArray> a_target;
  a_target->SetNumberOfComponents(a_source->GetNumberOfComponents());
  a_target->SetNumberOfTuples(mesh->GetNumberOfPoints());
  a_target->SetName(p.array_name.c_str());

  // Create the output array
  vtkNew<vtkDoubleArray> a_distance;
  a_distance->SetNumberOfComponents(1);
  a_distance->SetNumberOfTuples(mesh->GetNumberOfPoints());
  a_distance->SetName(p.dist_array_name.c_str());

  // Create a vtk locator for inside/outside checks
  vtkNew<vtkCellLocator> locator;
  locator->SetDataSet(tetra);
  locator->BuildLocator();

  // Iterate over the points in the mesh
  vtkSmartPointer<vtkPoints> mesh_pts = mesh->GetPoints();
  vtkNew<vtkIdList> idlist;
  double x[3], bbox[6];
  double weights[4], x_closest[3], pcoords[3];
  int sub_id;
  double dist2;
  double *tuple = new double[a_source->GetNumberOfComponents()];
  int n_inside = 0, n_close = 0;
  for(unsigned int i = 0; i < mesh_pts->GetNumberOfPoints(); i++)
    {
    // First try to just find the cell
    mesh_pts->GetPoint(i, x);
    vtkIdType id = locator->FindCell(x);

    // If that failed, look for nearby cells
    if(id >= 0)
      {
      a_distance->SetTuple1(i, 0.0);
      n_inside++;
      }
    else if(p.search_dist > 0)
      {
      // Find all the cells inside of the search radius
      for(unsigned int d = 0; d < 3; d++)
        {
        bbox[d * 2] = x[d] - p.search_dist;
        bbox[d * 2 + 1] = x[d] + p.search_dist;
        }
      idlist->Reset();
      locator->FindCellsWithinBounds(bbox, idlist);

      // For each tetrahedron check distance to it
      double mindist2 = p.search_dist * p.search_dist;
      for(unsigned int j = 0; j < idlist->GetNumberOfIds(); j++)
        {
        vtkCell *cell = tetra->GetCell(idlist->GetId(j));
        cell->EvaluatePosition(x, x_closest, sub_id, pcoords, dist2, weights);
        if(dist2 < mindist2)
          {
          mindist2 = dist2;
          id = idlist->GetId(j);
          }
        }

      if(id >= 0)
        {
        a_distance->SetTuple1(i, mindist2);
        n_close++;
        }
      }

    // Now fill output array
    if(id >= 0)
      {
      // Successfully found cell, take its value
      a_source->GetTuple(id, tuple);
      a_target->SetTuple(i, tuple);
      }
    else
      {
      a_distance->SetTuple1(i, numeric_limits<double>::quiet_NaN());
      for(unsigned int k = 0; k < a_source->GetNumberOfComponents(); k++)
        a_target->SetComponent(i, k, p.fill_value);
      }
    }

  // Print statistics
  cout << "Located: " << endl;
  cout << "  " << n_inside << " points inside tetrahedra" << endl;
  cout << "  " << n_close << " within tolerance of tetrahedra" << endl;
  cout << "  " << mesh->GetNumberOfPoints() - (n_inside + n_close) << " outside of tolerance" << endl;

  // Save the output mesh
  mesh->GetPointData()->AddArray(a_target);
  if(p.dist_array_name.length())
    mesh->GetPointData()->AddArray(a_distance);

  WriteMesh<TMesh>(mesh, p.fnOutput.c_str(), p.flagBinary);
  cout << "Wrote output mesh to " << p.fnOutput << endl;

  delete [] tuple;
  return 0;
}


int main(int argc, char *argv[])
{
  Parameters p;

  if(argc < 5)
    return usage();

  for(int ip = 1; ip < argc-4; ip++)
    {
    std::string arg = argv[ip];
    if(arg == "-d")
      {
      p.search_dist = atof(argv[++ip]);
      }
    else if(arg == "-D")
      {
      p.dist_array_name = argv[++ip];
      }
    else if(arg == "-b")
      {
      p.fill_value = atof(argv[++ip]);
      }
    else if(arg == "-B")
      {
      p.flagBinary = true;
      }
    else
      {
      cerr << "error: unrecognized parameter " << arg << endl;
      return usage();
      }
    }

  // Assign the positional parameters
  p.fnMesh = argv[argc-4];
  p.fnTetra = argv[argc-3];
  p.fnOutput = argv[argc-2];
  p.array_name = argv[argc-1];

  // Read the main mesh
  vtkDataReader *reader = vtkDataReader::New();
  reader->SetFileName(argv[argc-4]);
  reader->OpenVTKFile();
  reader->ReadHeader();

  bool isPolyData = true;
  // Is this a polydata?
  if(reader->IsFileUnstructuredGrid())
    {
    reader->Delete();
    return MeshTetraSample<vtkUnstructuredGrid>(p);
    }
  else if(reader->IsFilePolyData())
    {
    reader->Delete();
    return MeshTetraSample<vtkPolyData>(p);
    }
  else
    {
    reader->Delete();
    cerr << "Unsupported VTK data type in input file" << endl;
    return -1;
    }
}
