FIND_LIBRARY(IPOPT_LIBRARY NAMES ipopt)

IF(IPOPT_LIBRARY)

  OPTION(USE_MUMPS "Use MUMPS Solver" ON)
  OPTION(USE_HSL "Use HSL Solver" OFF)

  GET_FILENAME_COMPONENT(IPOPT_LIB_DIR ${IPOPT_LIBRARY} PATH)

  FIND_LIBRARY(IPOPT_GFORTRAN_LIB gfortran)
  FIND_LIBRARY(IPOPT_BLAS_LIB blas)
  FIND_LIBRARY(IPOPT_LAPACK_LIB lapack)

  IF(IPOPT_BLAS_LIB AND IPOPT_LAPACK_LIB)
    SET(IPOPT_LAPACK_LIBS ${IPOPT_LAPACK_LIB} ${IPOPT_BLAS_LIB})
    IF(IPOPT_GFORTRAN_LIB)
      SET(IPOPT_LAPACK_LIBS ${IPOPT_LAPACK_LIBS} ${IPOPT_GFORTRAN_LIB})
    ENDIF()
    SET(IPOPT_LAPACK_FOUND 1 CACHE INTERNAL "Lapack is available from IpOpt")
  ELSE()
    MESSAGE(SEND_ERROR "One of the IPOPT BLAS, LAPACK libs is not found")
  ENDIF()

  SET(IPOPT_LIBRARIES ${IPOPT_LIBRARY} ${IPOPT_LAPACK_LIBS})

  IF(USE_MUMPS)
    FIND_LIBRARY(IPOPT_METIS_LIBRARY coinmetis PATHS ${IPOPT_LIB_DIR})
    FIND_LIBRARY(IPOPT_MUMPS_LIBRARY coinmumps PATHS ${IPOPT_LIB_DIR})
    SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${IPOPT_METIS_LIBRARY} ${IPOPT_MUMPS_LIBRARY})
  ENDIF(USE_MUMPS)

  IF(USE_HSL)
    FIND_LIBRARY(IPOPT_METIS_LIBRARY coinmetis PATHS ${IPOPT_LIB_DIR})
    FIND_LIBRARY(IPOPT_HSL_LIBRARY coinhsl PATHS ${IPOPT_LIB_DIR})
    SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${IPOPT_METIS_LIBRARY} ${IPOPT_HSL_LIBRARY})
  ENDIF(USE_HSL)

  GET_FILENAME_COMPONENT(IPOPT_ROOT_DIR ${IPOPT_LIB_DIR} PATH)

  FIND_PATH(IPOPT_INCLUDE_DIR NAMES IpTNLP.hpp PATHS ${IPOPT_ROOT_DIR}/include/coin)

ENDIF(IPOPT_LIBRARY)


