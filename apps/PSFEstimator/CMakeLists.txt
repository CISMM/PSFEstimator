INCLUDE(${VTK_USE_FILE})

# use what QVTK built with
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
SET(DESIRED_QT_VERSION ${VTK_DESIRED_QT_VERSION} CACHE FILEPATH "")
FIND_PACKAGE(Qt4)
IF(QT_USE_FILE)
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
  SET(QT_LIBRARIES   ${QT_QT_LIBRARY})
ENDIF(QT_USE_FILE)

SET (QtGuiSrcs
  Main.cxx
  QBeadSpreadFunctionPropertyTableModel.cxx
  PSFEstimator.cxx
)

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
  ${QT_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

SET(UIS
  PSFEstimator.ui
  NewImageDialog.ui
)
SET (ApplicationHeaders
  PSFEstimator.h
  QBeadSpreadFunctionPropertyTableModel.h
)
QT4_WRAP_UI(UIHeaders ${UIS})
QT4_WRAP_CPP(MOCSrcs ${ApplicationHeaders} )
ADD_DEFINITIONS(-DQT_GUI_LIBS -DQT_CORE_LIB -DQT3_SUPPORT)
SET_SOURCE_FILES_PROPERTIES(${QtGuiSrcs} PROPERTIES
                          OBJECT_DEPENDS "${UIHeaders}")


# It's nice to have the ui in the windows project file...just double click on it
# and designer comes up on that ui file :)
IF (${CMAKE_BUILD_TOOL} MATCHES "msdev")
  SET (QtGuiSrcs ${QtGuiSrcs} ${UIS})
ENDIF (${CMAKE_BUILD_TOOL} MATCHES "msdev")
IF (${CMAKE_BUILD_TOOL} MATCHES "devenv")
  SET (QtGuiSrcs ${QtGuiSrcs} ${UIS})
ENDIF (${CMAKE_BUILD_TOOL} MATCHES "devenv")


# Workaround for rc.exe in Visual Studio 2005 being unable to handle many /I arguments.
If (WIN32)
  SET( WinIcon icon.res )
  ADD_CUSTOM_COMMAND( OUTPUT icon.res COMMAND rc.exe /fo icon.res ${CMAKE_CURRENT_SOURCE_DIR}/installer/Win32/icon.rc )
ENDIF (WIN32)

ADD_EXECUTABLE( PSFEstimator
  ${WinIcon}
  ${QtGuiSrcs}
  ${UISrcs}
  ${MOCSrcs}
)

TARGET_LINK_LIBRARIES( PSFEstimator
  ${QT_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  psfeModel
  psfeIO
  psfeVisualization
)

# Set package properties
SET( CPACK_PACKAGE_NAME "PSFEstimator" )
SET( CPACK_PACKAGE_VENDOR "University of North Carolina CISMM" )
SET( CPACK_PACKAGE_EXECUTABLES "PSFEstimator;PSF Estimator ${PSFEstimator_VERSION_STRING}" )
SET( CPACK_PACKAGE_DESCRIPTION_SUMMARY "PSF Estimator" )
SET( CPACK_PACKAGE_VERSION "${PSFEstimator_VERSION_STRING}" )
SET( CPACK_PACKAGE_VERSION_MAJOR "${PSFEstimator_MAJOR_NUMBER}" )
SET( CPACK_PACKAGE_VERSION_MINOR "${PSFEstimator_MINOR_NUMBER}" )
SET( CPACK_PACKAGE_VERSION_PATCH "${PSFEstimator_REVISION_NUMBER}" )
SET( CPACK_PACKAGE_INSTALL_DIRECTORY "PSF Estimator ${PSFEstimator_VERSION_STRING}" )
SET( CPACK_PACKAGE_INSTALL_REGISTRY_KEY "PSF Estimator ${PSFEstimator_VERSION_STRING}" )
SET( CPACK_RESOURCE_FILE_LICENSE "${PSFEstimator_SOURCE_DIR}/apps/PSFEstimator/installer/licenses/PSFEstimator-license.txt" )
SET( CPACK_RESOURCE_FILE_README  "${PSFEstimator_SOURCE_DIR}/apps/PSFEstimator/installer/licenses/README.txt" )

INCLUDE( installer/CMakeLists.txt )

INSTALL( TARGETS PSFEstimator
  RUNTIME DESTINATION bin
)

# This must always be last!
INCLUDE( CPack )

# This includes Microsoft runtime libraries when run in Windows dev environment.
INCLUDE( InstallRequiredSystemLibraries )
