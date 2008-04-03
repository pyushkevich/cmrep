#ifndef __MedialPDESites_h_
#define __MedialPDESites_h_

#include "MedialPDEMasks.h"
#include "GeometryDescriptor.h"
#include "MedialAtom.h"

class FDAbstractSite 
{
public:
  // Matrix and vector typedefs
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;
  
  // Initialize the site
  FDAbstractSite(FiniteDifferenceMask *xMask)
    { this->xMask = xMask; }

  // Compute the value of the encoded equation at the site
  virtual double ComputeEquation(const Mat &Y) = 0;

  // Compute the derivative of the encoded equation with respect to the
  // variables that are involved.
  virtual void ComputeDerivative(
    const Mat &Y, double *A, unsigned int iRow) = 0;
  
  // Specify the differential geometry at the site
  virtual void SetGeometry(GeometryDescriptor *gd, double rho) = 0;

  // Setup the equation for computing variational derivative of Phi for a
  // variation of the surface
  virtual void ComputeVariationalDerivativeMatrix(
    const Mat &Y, double *A, MedialAtom *xAtom) = 0;

  // Compute the right hand side for the derivative PDE
  virtual double ComputeVariationalDerivativeRHS(
    const Mat &Y, MedialAtom *xAtom, MedialAtom *dAtom) = 0;

  // Returns the number of boundaries that a site touches 
  // (0 - int, 1 - border, 2 - corner)
  virtual bool IsBorderSite() = 0;
  
  // virtual bool IsBorderSite(unsigned int dim) = 0;
  virtual double GetInitialValue()
    { return IsBorderSite() ? 0.0 : 1.0; }
    
 	// Print a report describing this site
 	virtual void PrintReport() = 0;
  
protected:

  // The finite differences mask
  FiniteDifferenceMask *xMask;
};

/**
 * This site represents internal grid points where the Laplacial equation holds
 */
class FDInternalSite : public FDAbstractSite
{
public:
  // Constructor: set the site
  FDInternalSite(FiniteDifferenceMask *xMask) :
    FDAbstractSite(xMask), xDerivative(xMask->Size()) {}
  
  // Set the geometry of the site
  void SetGeometry(GeometryDescriptor *gd, double rho);

  // Each site represents an equation F[N_i] = 0. This method computes F[N_i]
  double ComputeEquation(const Mat &Y);

  // This function computes dF[N_i]/d[N_i]
  void ComputeDerivative(const Mat &Y, double *A, unsigned int iRow);

  // This site touches no boundaries
  bool IsBorderSite() 
    { return false; }
  
  // bool IsBorderSite(unsigned int) 
  //  { return false; }

  // Get the column of the newton's matrix entry
  FiniteDifferenceMask *GetMask()
    { return xMask; }

  // Compute the variational derivative 
  void ComputeVariationalDerivativeMatrix(
    const Mat &Y, double *A, MedialAtom *xAtom);
  
  double ComputeVariationalDerivativeRHS(
    const Mat &Y, MedialAtom *xAtom, MedialAtom *dAtom);

  // Compute the variational derivative with respect to Rho-variation 
  // void ComputeVariationalDerivativeRho(const Mat &Y, 
  //  double *A, double *b, MedialAtom *xAtom, MedialAtom *dAtom);

	// Print a report
	void PrintReport();

protected:

  /** Coefficients of the finite differences in the equation */
  double Cuu, Cuv, Cvv, Cu, Cv;

  /** Finite difference derivatives of F, computed for variational differentiation */
  double F, Fu, Fv, Fuu, Fuv, Fvv;

  /** Derivatives with respect to each site */
  vnl_vector<double> xDerivative;
  
  /** The value of rho */
  double rho;

  friend class MedialPDESolver;
};


/**
 * This site represents internal grid points where the boundary condition
 * holds
 */
class FDBorderSite : public FDAbstractSite
{
public:
  // Constructor: set the site
  FDBorderSite(FiniteDifferenceMask *xMask) :
    FDAbstractSite(xMask), xDerivativeF(xMask->Size()),
    xDerivativeFu(xMask->Size()), xDerivativeFv(xMask->Size()) {}
  
  // Set the geometry of the site
  void SetGeometry(GeometryDescriptor *gd, double rho);

  // Each site represents an equation F[N_i] = 0. This method computes F[N_i]
  double ComputeEquation(const Mat &Y);

  // This function computes dF[N_i]/d[N_i]
  void ComputeDerivative(const Mat &Y, double *A, unsigned int iRow);

  // This site touches no boundaries
  bool IsBorderSite() 
    { return true; }
  
  // bool IsBorderSite(unsigned int) 
  //  { return false; }

  // Get the column of the newton's matrix entry
  FiniteDifferenceMask *GetMask()
    { return xMask; }

  // Compute the variational derivative 
  void ComputeVariationalDerivativeMatrix(
    const Mat &Y, double *A, MedialAtom *xAtom);
  
  double ComputeVariationalDerivativeRHS(
    const Mat &Y, MedialAtom *xAtom, MedialAtom *dAtom);

	// Print a report
	void PrintReport();

protected:

  /** Coefficients of the finite differences in the equation */
  double CuCu, CuCv, CvCv;

  /** Finite difference derivatives of F, computed for variational differentiation */
  double F, Fu, Fv;

  /** Derivatives with respect to each site */
  Vec xDerivativeFu, xDerivativeFv, xDerivativeF;
  
  /** The value of rho */
  double rho;

  friend class MedialPDESolver;
};




#endif
