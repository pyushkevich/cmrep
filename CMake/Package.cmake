# ----------------------------------------------------------------
# INSTALLATION AND PACKAGING with CPack
# ----------------------------------------------------------------

# On Win32, we must include the redistributable
IF(MSVC80 OR MSVC90)
  FIND_PROGRAM(VCREDIST_X86 vcredist_x86.exe)
  IF(VCREDIST_X86)
    INSTALL(FILES ${VCREDIST_X86} DESTINATION bin)
    SET(CPACK_NSIS_EXTRA_INSTALL_COMMANDS 
      "ExecWait '\\\"$INSTDIR\\\\bin\\\\vcredist_x86.exe\\\" /q:a'")
  ENDIF(VCREDIST_X86)
ENDIF(MSVC80 OR MSVC90)

# Allow package generation
SET(CPACK_PACKAGE_NAME "cmrep")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Continuous Medial Representation")
SET(CPACK_PACKAGE_VENDOR "picsl.upenn.edu")
SET(CPACK_PACKAGE_VERSION_MAJOR "${CMREP_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${CMREP_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${CMREP_VERSION_PATCH}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "cmrep-${CMREP_VERSION_FULL}")
SET(CPACK_NSIS_MODIFY_PATH ON)


# Shamelessly stolen from ParaView_
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "cmrep-${CMREP_VERSION_FULL}")
IF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
  EXEC_PROGRAM(uname ARGS "-m" OUTPUT_VARIABLE CMAKE_SYSTEM_PROCESSOR)
ENDIF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
IF(NOT DEFINED CPACK_SYSTEM_NAME)
  SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)
IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
  IF(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win64-${CMAKE_SYSTEM_PROCESSOR})
  ELSE(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win32-${CMAKE_SYSTEM_PROCESSOR})
  ENDIF(CMAKE_CL_64)
ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)
IF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
  SET(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")
ENDIF(NOT DEFINED CPACK_PACKAGE_FILE_NAME)

# Show GPL license
SET(CPACK_RESOURCE_FILE_LICENSE "${CMREP_SOURCE_DIR}/COPYING")

IF(WIN32 AND NOT UNIX)

  SET(CPACK_GENERATOR "ZIP")
  SET(CPACK_EXTENSION "zip")

ELSE(WIN32 AND NOT UNIX)

  # Set the generator to either STGZ or Apple
  IF(NOT APPLE)
    SET(CPACK_GENERATOR "TGZ")
    SET(CPACK_EXTENSION "tar.gz")
  ELSE(NOT APPLE)
    SET(CPACK_GENERATOR "ZIP")
    SET(CPACK_EXTENSION "zip")
  ENDIF(NOT APPLE)

ENDIF(WIN32 AND NOT UNIX)

#--------------------------------------------------------------------------------
# Uploading code to SourceForge
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Configure SCP

FIND_PROGRAM(SCP_PROGRAM NAMES scp DOC "Location of the scp program (optional)")
MARK_AS_ADVANCED(SCP_PROGRAM)

SET(SCP_ARGUMENTS "-v" CACHE STRING "Optional arguments to the scp command for uploads to SourceForge")
MARK_AS_ADVANCED(SCP_ARGUMENTS)

SET(SCP_USERNAME "" CACHE STRING "SourceForge.net account id for uploads")
MARK_AS_ADVANCED(SCP_USERNAME)

SET(NIGHTLY_TARGET "cmrep-nightly-${CPACK_SYSTEM_NAME}.${CPACK_EXTENSION}")
SET(EXPERIMENTAL_TARGET "cmrep-experimental-${CPACK_SYSTEM_NAME}.${CPACK_EXTENSION}")

SET(SCP_ROOT "frs.sourceforge.net:/home/frs/project/c/cm/cmrep/cmrep")

#--------------------------------------------------------------------------------
# Create targets

SET(CPACK_PACKAGE_FILE_NAME_WEXT "${CPACK_PACKAGE_FILE_NAME}.${CPACK_EXTENSION}")

ADD_CUSTOM_TARGET(cmrep_upload_nightly
  VERBATIM COMMAND "${SCP_PROGRAM}" ${SCP_ARGUMENTS}
  ${CPACK_PACKAGE_FILE_NAME_WEXT} ${SCP_USERNAME},itk-snap@${SCP_ROOT}/Nightly/${NIGHTLY_TARGET}
  DEPENDS ${CPACK_TARGET}
  WORKING_DIRECTORY ${CMREP_BINARY_DIR}
  COMMENT "Uploading package ${CPACK_PACKAGE_FILE_NAME_WEXT} to SourceForge.net as ${NIGHTLY_TARGET}")

ADD_CUSTOM_TARGET(cmrep_upload_experimental
  VERBATIM COMMAND "${SCP_PROGRAM}" ${SCP_ARGUMENTS} 
    ${CPACK_PACKAGE_FILE_NAME_WEXT} ${SCP_USERNAME},itk-snap@${SCP_ROOT}/Experimental/${EXPERIMENTAL_TARGET}
  DEPENDS ${CPACK_TARGET}
  WORKING_DIRECTORY ${CMREP_BINARY_DIR}
  COMMENT "Uploading package ${CPACK_PACKAGE_FILE_NAME_WEXT} to SourceForge.net as ${EXPERIMENTAL_TARGET}")


# This has to be at the end
INCLUDE(CPack)
