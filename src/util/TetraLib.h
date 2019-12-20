#ifndef __TetraLib_h_
#define __TetraLib_h_

#include "tetgen.h"


// Write tetgenio object to VTK
void WriteTetgenOutputAsUnstructuredMesh(tetgenio &out, const char *fname);
  

#endif
