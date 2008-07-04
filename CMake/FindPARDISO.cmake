# This only makes sense on unix platforms
IF(UNIX)
  
  # Set the library name on GCC or INTEL
  IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(PARDISO_LIB_NAME "pardiso_GNU_IA32")
  ELSE(CMAKE_COMPILER_IS_GNUCXX)
    SET(PARDISO_LIB_NAME "pardiso_INTEL_80_IA32")
  ENDIF(CMAKE_COMPILER_IS_GNUCXX)

  # Find the library
  FIND_LIBRARY(PARDISO_LIB ${PARDISO_LIB_NAME} 
    PATH "/lib" "/usr/lib" "/usr/local/lib" "/opt/pardiso/lib")

  # Set the variable
  IF(PARDISO_LIB)
    SET(PARDISO_FOUND 1 CACHE INTERNAL "PARDISO is available")
  ENDIF(PARDISO_LIB)

ELSE(UNIX)
  MESSAGE(SEND_ERROR
    "PARDISO is only available on the UNIX environment!")
ENDIF(UNIX)
