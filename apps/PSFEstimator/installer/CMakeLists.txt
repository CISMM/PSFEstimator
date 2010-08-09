IF ( WIN32 )
  INCLUDE( "${PSFEstimator_SOURCE_DIR}/apps/PSFEstimator/installer/Win32/CMake/Win32Package.txt" )
ENDIF ( WIN32 )

IF ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  INCLUDE( "${PSFEstimator_SOURCE_DIR}/apps/PSFEstimator/installer/MacOSX/CMake/MacOSXPackage.txt" )
ENDIF ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
