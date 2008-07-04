# Find the directory where MKL lives
FIND_LIBRARY(LAPACK_LIB lapack PATH "/lib" "/usr/lib" "/usr/local/lib" 
  DOC "Find the Lapack library")
FIND_LIBRARY(BLAS_LIB blas PATH "/lib" "/usr/lib" "/usr/local/lib" 
  DOC "Find the BLAS library")
FIND_LIBRARY(G2C_LIB g2c PATH "/lib" "/usr/lib" "/usr/local/lib" 
  DOC "Find the Fortran-to-C library")

# Set the flag that LAPACK is present
IF(LAPACK_LIB AND BLAS_LIB AND G2C_LIB)
  SET(LAPACK_LIBS ${LAPACK_LIB} ${BLAS_LIB} ${G2C_LIB})
  SET(LAPACK_FOUND 1 CACHE INTERNAL "Lapack is available")
ENDIF(LAPACK_LIB AND BLAS_LIB AND G2C_LIB)


