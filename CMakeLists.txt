PROJECT (VisualPSFOptimizer)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

# Set version information. Bump this up with every release.
SET (VisualPSFOptimizer_MAJOR_NUMBER "1")
SET (VisualPSFOptimizer_MINOR_NUMBER "0")
SET (VisualPSFOptimizer_REVISION_NUMBER "0")
SET (VisualPSFOptimizer_VERSION_STRING "${VisualPSFOptimizer_MAJOR_NUMBER}.${VisualPSFOptimizer_MINOR_NUMBER}.${VisualPSFOptimizer_REVISION_NUMBER}")

SET (EXECUTABLE_OUTPUT_PATH "${VisualPSFOptimizer_BINARY_DIR}/bin")
SET (LIBRARY_OUTPUT_PATH "${VisualPSFOptimizer_BINARY_DIR}/lib")

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
# Personal filter library and support classes.
#
# FibrinAnalysis is split up into several logical groups.
# apps      - a place for various applications related to the main application
#######################################
INCLUDE_DIRECTORIES (
    ${VisualPSFOptimizer_SOURCE_DIR}/lib
    ${VisualPSFOptimizer_SOURCE_DIR}/lib/filters
    ${VisualPSFOptimizer_SOURCE_DIR}/lib/filters/itk
    ${VisualPSFOptimizer_SOURCE_DIR}/lib/io
    ${VisualPSFOptimizer_SOURCE_DIR}/lib/model
    ${VisualPSFOptimizer_SOURCE_DIR}/lib/visualization
    ${VisualPSFOptimizer_BINARY_DIR}/lib)

ADD_SUBDIRECTORY(lib)

#######################################
# Applications
#######################################
ADD_SUBDIRECTORY(apps)


#######################################
# Testing components
#######################################
#OPTION (BUILD_IT_TESTS "Include VisualPSFOptimizer unit tests" ON)
#IF (BUILD_IT_TESTS)
#    SUBDIRS(TestSuite Tests)
#ENDIF (BUILD_IT_TESTS)
