PROJECT (PSFEstimator)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

# Set version information. Bump this up with every release.
SET (PSFEstimator_MAJOR_NUMBER "1")
SET (PSFEstimator_MINOR_NUMBER "0")
SET (PSFEstimator_REVISION_NUMBER "0")
SET (PSFEstimator_VERSION_STRING "${PSFEstimator_MAJOR_NUMBER}.${PSFEstimator_MINOR_NUMBER}.${PSFEstimator_REVISION_NUMBER}")

SET (EXECUTABLE_OUTPUT_PATH "${PSFEstimator_BINARY_DIR}/bin")
SET (LIBRARY_OUTPUT_PATH "${PSFEstimator_BINARY_DIR}/lib")

OPTION (CMAKE_VERBOSE_MAKEFILE  "Enable/Disable verbose compiler output" ON)
OPTION (CMAKE_COLOR_MAKEFILE "Enable/Disable color cues when building" ON)
MARK_AS_ADVANCED (CLEAR CMAKE_VERBOSE_MAKEFILE CMAKE_COLOR_MAKEFILE)

#######################################
# Include ITK
#######################################

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
     "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

#######################################
# Include VTK
#######################################

FIND_PACKAGE(VTK)
IF (VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR
     "Cannot build without VTK.  Please set VTK_DIR.")
ENDIF (VTK_FOUND)

#######################################
# Project library and support classes.
#
# PSFEstimator is split up into several logical groups.
#######################################
INCLUDE_DIRECTORIES (
    ${PSFEstimator_SOURCE_DIR}/NanoITK/Code/Algorithms
    ${PSFEstimator_SOURCE_DIR}/NanoITK/Code/BasicFilters
    ${PSFEstimator_SOURCE_DIR}/NanoITK/Code/Common
    ${PSFEstimator_SOURCE_DIR}/lib
    ${PSFEstimator_SOURCE_DIR}/lib/IO
    ${PSFEstimator_SOURCE_DIR}/lib/Model
    ${PSFEstimator_SOURCE_DIR}/lib/Visualization
    ${PSFEstimator_BINARY_DIR}/lib
)

ADD_SUBDIRECTORY(lib)

#######################################
# Applications
#######################################
ADD_SUBDIRECTORY(apps)

#######################################
# Documentation
#######################################
ADD_SUBDIRECTORY(docs)
