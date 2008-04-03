// Definitions of the callback functions
typedef void (*sparse_f)(const double *x, double *fx, void *param);
typedef void (*sparse_fdf)(const double *x, double *fx, double *j, void *param);

/**
 * Newton solver that does backtracking
 */
int sparse_backtracking_newton_solver(
  size_t n,             // Number of rows
  size_t *rowindex,     // Sparse matrix row index
  size_t *colindex,     // Sparse matric column index
  double *function,     // Data array for the function (NULL to use internal array)
  double *jacobian,     // Data array for the jacobian (NULL to use internal array)
  sparse_f f,           // Callback function
  sparse_fdf fdf,       // Callback Jacobian + function 
  double *xinit,        // Initial solution
  size_t niter,         // Number of iterations
  void *param)
{
  // Get the size of the sparse Jacobian matrix
  size_t m = rowindex[n];

  // Allocate the Jacobian array if necessary
  double *J = (jacobian == NULL) ? new double[m] : jacobian;

  // Allocate the function value array if necessary
  double *f = (function == NULL) ? new double[n] : function;

  // Compute the Jacobian matrix
  fdf(xinit, f, J, param);

  // Compute the value of the solution
  double F = dot_product(f, f, n);

  // Compute the Jacobian update direction
  



}

