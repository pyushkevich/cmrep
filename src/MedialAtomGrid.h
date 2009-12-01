#ifndef __MedialAtomGrid_h_
#define __MedialAtomGrid_h_

#include "MedialAtom.h"
#include "MedialAtomIterators.h"
#include "SubdivisionSurface.h"

using namespace std;

/** 
 * Helper function to access a boundary site in an atom array using a
 * boundary point iterator 
 */
inline BoundaryAtom &
GetBoundaryPoint(MedialBoundaryPointIterator &itBoundary, MedialAtom *xAtoms)
{
  return xAtoms[itBoundary.GetAtomIndex()].xBnd[itBoundary.GetBoundarySide()];
}

/**
 * A helper function that accesses the boundary site for a triangle iterator
 */
inline BoundaryAtom &
GetBoundaryPoint(MedialBoundaryTriangleIterator &itTriangle, MedialAtom *xAtoms, size_t i)
{
  return xAtoms[itTriangle.GetAtomIndex(i)].xBnd[itTriangle.GetBoundarySide()];
}

/** Helper function to compute the coordinate of an internal medial point
 * that is referenced by an iterator */
SMLVec3d GetInternalPoint(MedialInternalPointIterator &it, MedialAtom *xAtoms);

// This method computes the weights for integration over the domain of the medial
// surface. Because the domain may be non-uniform, we must scale all integrals in
// du dv by these weights
double ComputeMedialDomainAreaWeights( MedialIterationContext *xGrid, 
  MedialAtom *xAtoms, double *xWeights);

/**
 * This function computes the area associated with each triangle on the
 * boundary of the object and assigns a third of this area to each of the 
 * vertices of the triangle. Thus for each vertex on the boundary a weight
 * is generated. The total of these weights, equal to the area of the boundary
 * surface is the return value 
 */
double ComputeMedialBoundaryAreaWeights( 
  MedialIterationContext *xGrid, MedialAtom *xAtoms,
  SMLVec3d *xNormals, double *xWeights);

double ComputeMedialInternalVolumePartialDerivativeWeights(
  MedialIterationContext *xGrid, MedialAtom *dAtoms, 
  SMLVec3d *xPoints, SMLVec3d *dPoints, size_t nCuts, 
  double *dWeights, double *dInternalProfileWeights);

/**
 * Compute the derivative of the medial boundary weights */
double ComputeMedialBoundaryAreaPartialDerivative(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, MedialAtom *dAtoms, 
  double *xWeights, double *dWeights);

/** Compute integration weights over the medial surface */
double ComputeMedialSurfaceAreaWeights( 
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  SMLVec3d *xNormals, double *xWeights);

double ComputeMedialSurfaceAreaPartialDerivative(
  MedialIterationContext *xGrid, 
  MedialAtom *xAtoms, MedialAtom *dAtoms, 
  double *xWeights, double *dWeights);

/**
 * This is a fast method for volume element computation. It uses a rough
 * approximation of volume element as
 *   V(u,v,w) = r * ((1-w) aelt[medial] + w aelt[bnd])
 * I've played around in mathematica and it seems that the difference 
 * between this and the true volume element
 *   V(u,v,w) = - Zu x Zv . Zw, Z = (1-w) X[medial] + w * X[bnd]
 * is pretty tiny. Perhaps it's not so tiny near the boundary of the 
 * medial surface, but everything is messed up there anyway :)
 */
double FastComputeMedialInternalVolumeWeights(
  MedialIterationContext *xGrid, 
  MedialAtom *xAtoms,
  double *xMedialWeights, 
  double *xBoundaryWeights, 
  size_t nCuts,
  double *xWeights);

double FastComputeMedialInternalVolumeWeightsVariationalDerivative(
  MedialIterationContext *xGrid, 
  MedialAtom *xAtoms, MedialAtom *dAtoms,
  double *xMedialWeights, double *dMedialWeights,
  double *xBoundaryWeights, double *dBoundaryWeights,
  double *xWeights, size_t nCuts,
  double *dWeights);


/** 
 * Interpolate a list of internal medial points from the m-rep.
 */
void ComputeMedialInternalPoints(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, size_t nCuts, SMLVec3d *xPoints);

void ComputeMedialInternalPointsPartialDerivative(
  MedialIterationContext *xGrid, MedialAtom *dAtoms, size_t nCuts, SMLVec3d *dPoints);

// This is a measure that can be computed over a volume (just a R3 function)
class EuclideanFunction {
public:
  virtual ~EuclideanFunction () {}

  virtual double Evaluate(const SMLVec3d &x) = 0;
  virtual void ComputeGradient(const SMLVec3d &x, SMLVec3d &G)
    { G.fill(0.0); }
  virtual double ComputeFunctionAndGradient(const SMLVec3d &x, SMLVec3d &G)
    {
    this->ComputeGradient(x,G);
    return this->Evaluate(x);
    }
};

inline double ScalarTripleProduct(
  const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C)
{
  return
    A[0] * (B[1] * C[2] - B[2] * C[1]) + 
    A[1] * (B[2] * C[0] - B[0] * C[2]) + 
    A[2] * (B[0] * C[1] - B[1] * C[0]);
}

inline double ScalarTripleProductDerivative(
  const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C,
  const SMLVec3d &DA, const SMLVec3d &DB, const SMLVec3d &DC)
{
  return
    DA[0] * (B[1] * C[2] - B[2] * C[1]) + 
    DA[1] * (B[2] * C[0] - B[0] * C[2]) + 
    DA[2] * (B[0] * C[1] - B[1] * C[0]) +
    A[0] * (DB[1] * C[2] - DB[2] * C[1] + B[1] * DC[2] - B[2] * DC[1]) + 
    A[1] * (DB[2] * C[0] - DB[0] * C[2] + B[2] * DC[0] - B[0] * DC[2]) + 
    A[2] * (DB[0] * C[1] - DB[1] * C[0] + B[0] * DC[1] - B[1] * DC[0]);
}


/**
 * Compute an area of a triangle
 */
inline double TriangleArea(const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C)
{
  return 0.5 * vnl_cross_3d(B - A, C - A).magnitude();
}

/** 
 * Compute the partial derivative of the triangle area with respect to 
 * some parameter p.
 */
inline double TriangleAreaPartialDerivative(
  const SMLVec3d &A, const SMLVec3d &B, const SMLVec3d &C, 
  const SMLVec3d &Ap, const SMLVec3d &Bp, const SMLVec3d &Cp,
  double xArea)
{
  SMLVec3d N = 0.5 * vnl_cross_3d(B - A, C - A);
  return - 0.5 * 
    ( dot_product(Ap, vnl_cross_3d(N, (B - C))) + 
      dot_product(Bp, vnl_cross_3d(N, (C - A))) +
      dot_product(Cp, vnl_cross_3d(N, (A - B))) ) / xArea;
}


void TestTriangleAreaPartialDerivative();


/**
 * Compute the volume of a prism formed by two triangles
 */
inline double PrismVolume( SMLVec3d A1, SMLVec3d B1, SMLVec3d C1, 
  SMLVec3d A2, SMLVec3d B2, SMLVec3d C2) 
{
  // TODO: Get/derive the correct formula!!!
  double xArea1 = TriangleArea(A1, B1, C1);
  double xArea2 = TriangleArea(A2, B2, C2);
  double d = ((A1 + B1 + C1) - (A2 + B2 + C2)).magnitude() / 3.0;
  return d * 0.5 * (xArea1 + xArea2);
}

/** Compute the volume of a convex cell. If the cell is not convex, the result
 * will most likely be nonsense */
double CellVolume(
  const SMLVec3d &X000, const SMLVec3d &X001, 
  const SMLVec3d &X010, const SMLVec3d &X011, 
  const SMLVec3d &X100, const SMLVec3d &X101, 
  const SMLVec3d &X110, const SMLVec3d &X111);


/**
 * Here is a similar formula for a wedge volume. The wedge here is treated as
 * a half of a parallellepiped. The vectors forming its edges are computed by 
 * averaging the contributing points
 */
double WedgeVolume(
  const SMLVec3d &A0, const SMLVec3d &B0, const SMLVec3d &C0,
  const SMLVec3d &A1, const SMLVec3d &B1, const SMLVec3d &C1);

/**
 * This function computes the volume associated with each cell in the medial
 * interior and assigns a portion of this cell to each internal point that is
 * the vertex of the cell. The sum of these resulting internal point weights
 * is equal to the volume of the m-rep
 */
double ComputeMedialInternalVolumeWeights(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  SMLVec3d *xPoints, size_t nCuts, 
  double *xWeights, double *xProfileIntervalWeights = NULL);

void ComputeMedialVolumeIntegrationWeights (
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  double *xVolumeEltQuadTerm, double *xVolumeEltSlope, double *xVolumeEltIntercept);

/** This method integrates any given function over the interior of mrep */
double IntegrateFunctionOverInterior (
  MedialIterationContext *xGrid, SMLVec3d *xPoints, 
  double *xWeights, size_t nCuts, EuclideanFunction *fMatch);

/** Integrate any function over the interior of a medial atom */
double IntegrateFunctionOverBoundary (
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  std::vector<double> &xWeights, EuclideanFunction *fMatch);

/** Compute gradient of a function over the boundary */
double ComputeFunctionJetOverBoundary(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  std::vector<double> &xWeights, EuclideanFunction *fMatch, SMLVec3d *xOutGradient);


/** 
 * Integrate a function over the boundary, making sure that irregular quads are 
 * super-sampled. This results in better quality integration, but is more costly
 * to compute.
 */
double AdaptivelyIntegrateFunctionOverBoundary(
  MedialIterationContext *xGrid, MedialAtom *xAtoms, 
  double xMinQuadArea, EuclideanFunction *fMatch);

/**
 * This is an extension of the loop tangent scheme that works with
 * medial atoms rather than arrays of data. It's a convenience class
 * and this seemed like the best place to stick it in
 */

class MedialAtomLoopScheme : public LoopTangentScheme
{
public:

  // Helper methods to get all the partials given a list of atoms
  SMLVec3d Xu(size_t v, MedialAtom *a) const { return TangentX(0, v, a); }
  SMLVec3d Xv(size_t v, MedialAtom *a) const { return TangentX(1, v, a); }
  SMLVec3d Xuu(size_t v, MedialAtom *a) const { return TangentXu(0, v, a); }
  SMLVec3d Xuv(size_t v, MedialAtom *a) const { return TangentXv(0, v, a); }
  SMLVec3d Xvu(size_t v, MedialAtom *a) const { return TangentXu(1, v, a); }
  SMLVec3d Xvv(size_t v, MedialAtom *a) const { return TangentXv(1, v, a); }

  double Fu(size_t v, MedialAtom *a) const { return PartialF(0, v, a); }
  double Fv(size_t v, MedialAtom *a) const { return PartialF(1, v, a); }
  double Fuu(size_t v, MedialAtom *a) const { return PartialFu(0, v, a); }
  double Fuv(size_t v, MedialAtom *a) const { return PartialFv(0, v, a); }
  double Fvu(size_t v, MedialAtom *a) const { return PartialFu(1, v, a); }
  double Fvv(size_t v, MedialAtom *a) const { return PartialFv(1, v, a); }

  double Ru(size_t v, MedialAtom *a) const { return PartialR(0, v, a); }
  double Rv(size_t v, MedialAtom *a) const { return PartialR(1, v, a); }
  // double Ruu(size_t v, MedialAtom *a) const { return PartialRu(0, v, a); }
  // double Ruv(size_t v, MedialAtom *a) const { return PartialRv(0, v, a); }
  // double Rvu(size_t v, MedialAtom *a) const { return PartialRu(1, v, a); }
  // double Rvv(size_t v, MedialAtom *a) const { return PartialRv(1, v, a); }
  
  double AreaEltDu(size_t v, MedialAtom *a) const { return PartialAreaElt(0, v, a); }
  double AreaEltDv(size_t v, MedialAtom *a) const { return PartialAreaElt(1, v, a); }

  SMLVec3d Nu(size_t v, MedialAtom *a) const { return NormalX(0, v, a); }
  SMLVec3d Nv(size_t v, MedialAtom *a) const { return NormalX(1, v, a); }

  // Derivative of the global parameterization (U, V) with respect to the local
  // parameterization
  double Uu(size_t v, MedialAtom *a) const { return PartialU(0, v, a); }
  double Uv(size_t v, MedialAtom *a) const { return PartialU(1, v, a); }
  double Vu(size_t v, MedialAtom *a) const { return PartialV(0, v, a); }
  double Vv(size_t v, MedialAtom *a) const { return PartialV(1, v, a); }
  
private:

  SMLVec3d TangentX(size_t d, size_t v, MedialAtom *atoms) const
    {    
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    SMLVec3d y(0.0);
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].X;

    return y;
    }

  SMLVec3d NormalX(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    SMLVec3d y(0.0);
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].N;

    return y;    
    }

  /**
   * This method is used to compute second order partials Xuu and Xvu
   */
  SMLVec3d TangentXu(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    SMLVec3d y(0.0);
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].Xu;

    return y;    
    } 

  /**
   * This method is used to compute second order partials Xuv and Xvv
   */
  SMLVec3d TangentXv(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    SMLVec3d y(0.0);
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].Xv;

    return y;  
    } 

  double PartialU(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].u;

    return y;  
    }

  double PartialV(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].v;

    return y;  
    }


  double PartialF(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].F;

    return y;  
    }

  double PartialFv(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].Fv;

    return y;  
    }

  double PartialFu(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].Fu;

    return y;  
    }

  double PartialR(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].R;

    return y;  
    }

  double PartialRv(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].Ru;

    return y;  
    }

  double PartialRu(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].Rv;

    return y;  
    }

  double PartialAreaElt(size_t d, size_t v, MedialAtom *atoms) const
    {
    const size_t *ri = W.GetRowIndex();
    const size_t *ci = W.GetColIndex();

    double y = 0.0;
    for(size_t i = ri[v]; i < ri[v+1]; ++i)
      y += W.GetSparseData()[i].w[d] * atoms[ci[i]].aelt;

    return y;  
    }

};


#endif
