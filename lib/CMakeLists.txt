ADD_SUBDIRECTORY(IO)
ADD_SUBDIRECTORY(Model)
ADD_SUBDIRECTORY(Visualization)

CONFIGURE_FILE( "${PSFEstimator_SOURCE_DIR}/lib/Version.h.in"
  "${PSFEstimator_BINARY_DIR}/lib/Version.h" )

