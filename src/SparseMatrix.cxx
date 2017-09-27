#include "SparseMatrix.h"
#include "SparseMatrix.txx"
#include <vnl/vnl_sparse_matrix.h>
#include "BruteForceSubdivisionMedialModel.h"

// Need an instantiation of sparse matrix
template class vnl_sparse_matrix<int>;
template class ImmutableSparseArray<double>;
template class ImmutableSparseArray<int>;
template class ImmutableSparseMatrix<double>;
template class ImmutableSparseMatrix<int>;
template class ImmutableSparseArray<std::pair<double, double> >;
template class ImmutableSparseArray<BruteForceSubdivisionMedialModel::NonvaryingAtomTerms>;
