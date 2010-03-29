INCLUDE(${VTK_USE_FILE})

# use what QVTK built with
SET(QT_MOC_EXECUTABLE ${VTK_QT_MOC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_UIC_EXECUTABLE ${VTK_QT_UIC_EXECUTABLE} CACHE FILEPATH "")
SET(QT_QMAKE_EXECUTABLE ${VTK_QT_QMAKE_EXECUTABLE} CACHE FILEPATH "")
SET(DESIRED_QT_VERSION ${VTK_DESIRED_QT_VERSION} CACHE FILEPATH "")
FIND_PACKAGE(Qt)
IF(QT_USE_FILE)
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
  SET(QT_LIBRARIES   ${QT_QT_LIBRARY})
ENDIF(QT_USE_FILE)

SET (QtGuiSrcs
  Main.cxx
  QPointSpreadFunctionPropertyTableModel.cxx
  VisualPSFOptimizer.cxx
)

# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
  ${QT_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

SET(UIS
  VisualPSFOptimizer.ui
  NewImageDialog.ui
)
SET (ApplicationHeaders
  VisualPSFOptimizer.h
  QPointSpreadFunctionPropertyTableModel.h
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

IF ( WIN32 )
  SET( WinResources "${VisualPSFOptimizer_SOURCE_DIR}/apps/VisualPSFOptimizer/installer/Win32/icon.rc" )
ENDIF ( WIN32 )

ADD_EXECUTABLE( VisualPSFOptimizer ${QtGuiSrcs} ${UISrcs} ${MOCSrcs} ${WinResources} )

TARGET_LINK_LIBRARIES( VisualPSFOptimizer
  QVTK
  ${QT_LIBRARIES}
  vtksys
  vtkCommon
  vtkFiltering
  vtkGraphics
  vtkHybrid
  vtkImaging
  vtkRendering
  vtkVolumeRendering
  itksys
  ITKCommon
  ITKBasicFilters
  ITKIO
  ITKNumerics
  vpsfoFilter
  vpsfoIO
  vpsfoModel
  vpsfoVisualization
)

# Set package properties
SET( CPACK_PACKAGE_NAME "VisualPSFOptimizer" )
SET( CPACK_PACKAGE_VENDOR "University of North Carolina CISMM" )
SET( CPACK_PACKAGE_EXECUTABLES "VisualPSFOptimizer;Visual PSF Optimizer ${VisualPSFOptimizer_VERSION_STRING}" )
SET( CPACK_PACKAGE_DESCRIPTION_SUMMARY "Visual PSF Optimizer" )
SET( CPACK_PACKAGE_VERSION "${VisualPSFOptimizer_VERSION_STRING}" )
SET( CPACK_PACKAGE_VERSION_MAJOR "${VisualPSFOptimizer_MAJOR_NUMBER}" )
SET( CPACK_PACKAGE_VERSION_MINOR "${VisualPSFOptimizer_MINOR_NUMBER}" )
SET( CPACK_PACKAGE_VERSION_PATCH "${VisualPSFOptimizer_REVISION_NUMBER}" )
SET( CPACK_PACKAGE_INSTALL_DIRECTORY "Visual PSF Optimizer ${VisualPSFOptimizer_VERSION_STRING}" )
SET( CPACK_PACKAGE_INSTALL_REGISTRY_KEY "Visual PSF Optimizer ${VisualPSFOptimizer_VERSION_STRING}" )
SET( CPACK_RESOURCE_FILE_LICENSE "${VisualPSFOptimizer_SOURCE_DIR}/apps/VisualPSFOptimizer/installer/licenses/VisualPSFOptimizer-license.txt" )
SET( CPACK_RESOURCE_FILE_README  "${VisualPSFOptimizer_SOURCE_DIR}/apps/VisualPSFOptimizer/installer/licenses/README.txt" )

INCLUDE( installer/CMakeLists.txt )

INSTALL( TARGETS VisualPSFOptimizer
  RUNTIME DESTINATION bin
)

# This must always be last!
INCLUDE( CPack )

# This includes Microsoft runtime libraries when run in Windows dev environment.
INCLUDE( InstallRequiredSystemLibraries )
