#ifndef __MedialPDERenderer_h_
#define __MedialPDERenderer_h_

#include "glengine.h"

class GenericMedialModel;

void glDrawWireframeElements(unsigned short width, unsigned short height);
void glDrawQuadStripElements(unsigned short width,unsigned short height);

class PDESplineRenderer : public GLDisplayListRenderer
{
public:
  // Constructor
  PDESplineRenderer(GenericMedialModel *solver);

  // Build the display list
  virtual void build();

private:
  // INternal methods
  void DrawInternalPoints( size_t nCuts );

  GenericMedialModel *solver;
  GLMaterial *matMedial, *matBoundary;
};

#endif
