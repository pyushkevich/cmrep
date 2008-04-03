/**
 * This class represents a sparse matrix whose zero elements can not
 * be altered once they have been set. This matrix can only be 
 * initialized using a mutable sparse matrix, such as those provided 
 * by VNL
 */
#ifndef __SparseMatrix_h_
#define __SparseMatrix_h_

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_sparse_matrix.h>
#include <iostream>
#include <list>
#include <vector>
#include <stdexcept>

template<class TVal>
class ImmutableSparseArray
{
public:
  // Typedefs
  typedef ImmutableSparseArray<TVal> Self;

  // Typedef for import from VNL
  typedef vnl_sparse_matrix<TVal> VNLSourceType;

  // Typedefs for import from STL structures
  typedef std::pair<size_t, TVal> STLEntryType;
  typedef std::list<STLEntryType> STLRowType;
  typedef std::vector<STLRowType> STLSourceType;

  // Default constructor
  ImmutableSparseArray();

  // Destructor
  ~ImmutableSparseArray();

  // Assignment operator for VNL
  void SetFromVNL(VNLSourceType &src);

  // Assignment operator that takes an array of lists
  void SetFromSTL(STLSourceType &src, size_t nColumns);

  // Set all the arrays in the matrix to external pointers. The caller must
  // relinquish the control of the pointers to the sparse matrix, which will
  // delete the data at some point
  void SetArrays(size_t rows, size_t cols, size_t *xRowIndex, size_t *xColIndex, TVal *data);

  // Pointers to the data stored inside the matrix
  size_t *GetRowIndex() { return xRowIndex; }
  size_t *GetColIndex() { return xColIndex; }
  TVal *GetSparseData() { return xSparseValues; }

  const size_t *GetRowIndex() const { return xRowIndex; }
  const size_t *GetColIndex() const { return xColIndex; }
  const TVal *GetSparseData() const { return xSparseValues; }

  class ConstRowIterator {
  public:
    ConstRowIterator(const Self *p, size_t row)
      { this->p = p; iStart = iPos = p->xRowIndex[row]; iEnd = p->xRowIndex[row+1]; }
    
    ConstRowIterator()
      { p = NULL; iPos = iEnd = iStart = 0; }

    bool IsAtEnd()
      { return iPos == iEnd; }

    const TVal &Value()
      { return p->xSparseValues[iPos]; }

    size_t Column()
      { return p->xColIndex[iPos]; }

    size_t Size()
      { return iEnd - iStart; }

    ConstRowIterator &operator++()
      { ++iPos; return *this; }
    
  private:
    const Self *p;
    size_t iPos, iEnd, iStart;
  };

  // Row iterator goes through nonzero elements of rows
  class RowIterator {
  public:
    RowIterator(Self *p, size_t row)
      { this->p = p; iStart = iPos = p->xRowIndex[row]; iEnd = p->xRowIndex[row+1]; }

    RowIterator()
      { p = NULL; iPos = iEnd = iStart = 0; }
    
    bool IsAtEnd()
      { return iPos == iEnd; }

    TVal &Value()
      { return p->xSparseValues[iPos]; }

    size_t Size()
      { return iEnd - iStart; }

    size_t Column()
      { return p->xColIndex[iPos]; }

    RowIterator &operator++()
      { ++iPos; return *this; }
    
  private:
    Self *p;
    size_t iPos, iEnd, iStart;
  };
    
  // Get the row iterator
  RowIterator Row(size_t iRow)
    { return RowIterator(this, iRow); } 

  // Get the row iterator
  ConstRowIterator Row(size_t iRow) const
    { return ConstRowIterator(this, iRow); } 

  // Set i-th non-zero value in a row
  TVal &GetValueBySparseIndex(size_t iRow, size_t iNZInRow)
    { return xSparseValues[xRowIndex[iRow] + iNZInRow]; }

  // Get the number of sparse values
  size_t GetNumberOfSparseValues() const { return nSparseEntries; }

  // Get the number of rows and columns
  size_t GetNumberOfRows() const { return nRows; }
  size_t GetNumberOfColumns() const { return nColumns; }  

  // Reset the matrix (clear all data, revert to initialized state)
  void Reset();

  // Set all the values in the matrix to some value
  void Fill(const TVal &value)
    { for(size_t i = 0; i < nSparseEntries; i++) xSparseValues[i] = value; }

  // A copy operator that actually copies the data
  Self & operator = (const Self &src);

  // A copy constructor
  ImmutableSparseArray(const ImmutableSparseArray<TVal> &src);

protected:
  
  /** Representation for the sparse matrix */
  TVal *xSparseValues;
  size_t *xRowIndex, *xColIndex, nRows, nColumns, nSparseEntries;

  // Iterator is our friend
  friend class RowIterator;
};


template<class TVal>
class ImmutableSparseMatrix : public ImmutableSparseArray<TVal>
{
public:
  // Typedefs
  typedef ImmutableSparseMatrix<TVal> Self;
  typedef vnl_vector<TVal> Vec;

  // Compare two matrices
  bool operator == (const Self &B);

  // Get a regular VNL matrix
  vnl_matrix<TVal> GetDenseMatrix() const;

  // Set the matrix to identity
  void SetIdentity(size_t n);

  // Compute the matrix product C = A * B
  static void Multiply(Self &C, const Self &A, const Self &B);

  // Compute the matrix product c = A^t * b
  Vec MultiplyByVector(const Vec &b) const;

  // Compute the matrix product c = A^t * b
  Vec MultiplyTransposeByVector(const Vec &b) const;

  // This method initializes the matrix A^t A (it simply creates the 
  // structure of the matrix for future fast computations)
  static void InitializeATA(Self &ATA, const Self &A);

  // This matrix computes the product ATA. The matrix ATA should be
  // initialized before calling this method
  static void ComputeATA(Self &ATA, const Self &A);

  // Print to the standard stream
  void PrintSelf(std::ostream &out) const;
};

// Print the matrix to an output stream
template<class TVal>
std::ostream& operator << (std::ostream &out, const ImmutableSparseMatrix<TVal> &A)
  { A.PrintSelf(out); }


  
#endif
