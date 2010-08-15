#ifndef __ReadWriteVTK_h_
#define __ReadWriteVTK_h_

#include "vtkPolyData.h"
#include <string>

void WriteVTKData(vtkPolyData *data, std::string fn);
vtkPolyData *ReadVTKData(std::string fn);

#endif
