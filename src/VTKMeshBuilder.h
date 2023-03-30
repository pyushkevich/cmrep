#ifndef __VTKMeshBuilder_h_
#define __VTKMeshBuilder_h_

#include <vtkSmartPointer.h>
#include <vnl_matrix.h>
#include <vnl_vector.h>
#include <cstdarg>

class TriangleMesh;

/**
 * A helper class for making VTK meshes and adding stuff to them using vnl matrices
 * and vectors.
 */
template <class TDataSet>
class VTKMeshBuilder
{
public:

  VTKMeshBuilder();

  VTKMeshBuilder(TDataSet* other);

  TDataSet *GetData() const { return pd; }

  void SetPoints(const vnl_matrix<double> &x);

  void SetTriangles(const TriangleMesh &mesh);
  void SetTriangles(const vnl_matrix<unsigned int> &tri);

  void SetNormals(const vnl_matrix<double> &x);

  template<typename TArray> void AddArray(const TArray &arr, const char *format, ...)
  {
    char buffer[256];
    va_list args;
    va_start (args, format);
    vsnprintf (buffer, 256, format, args);
    AddArrayInternal(arr, buffer);
    va_end (args);
  }

  template<typename TArray> void AddCellArray(const TArray &arr, const char *format, ...)
  {
    char buffer[256];
    va_list args;
    va_start (args, format);
    vsnprintf (buffer, 256, format, args);
    AddArrayInternal(arr, buffer, true);
    va_end (args);
  }

  void Save(const char *fn);

protected:
  vtkSmartPointer<TDataSet> pd;

  void AddArrayInternal(const vnl_matrix<double> &x, const char *name, bool cell=false);
  void AddArrayInternal(const vnl_vector<double> &x, const char *name, bool cell=false);
  void AddArrayInternal(const vnl_matrix<int> &x, const char *name, bool cell=false);
  void AddArrayInternal(const vnl_vector<int> &x, const char *name, bool cell=false);

};



#endif
