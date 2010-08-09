IF ( WIN32 )
  INCLUDE( "${VisualPSFOptimizer_SOURCE_DIR}/apps/VisualPSFOptimizer/installer/Win32/CMake/Win32Package.txt" )
ENDIF ( WIN32 )

IF ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  INCLUDE( "${VisualPSFOptimizer_SOURCE_DIR}/apps/VisualPSFOptimizer/installer/MacOSX/CMake/MacOSXPackage.txt" )
ENDIF ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
