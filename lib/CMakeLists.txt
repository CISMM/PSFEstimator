ADD_SUBDIRECTORY(Filters)
ADD_SUBDIRECTORY(IO)
ADD_SUBDIRECTORY(Model)
ADD_SUBDIRECTORY(Visualization)

CONFIGURE_FILE( "${VisualPSFOptimizer_SOURCE_DIR}/lib/Version.h.in"
  "${VisualPSFOptimizer_BINARY_DIR}/lib/Version.h" )
  
