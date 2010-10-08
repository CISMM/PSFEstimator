# Configure the documentation guide with the library/program version number
CONFIGURE_FILE( "${PSFEstimator_SOURCE_DIR}/docs/Configuration.in.tex"
  "${PSFEstimator_BINARY_DIR}/docs/Configuration.tex" @ONLY )

# This is a simple file copy
CONFIGURE_FILE( "${PSFEstimator_SOURCE_DIR}/docs/PSFEstimator-UserGuide.tex"
  "${PSFEstimator_BINARY_DIR}/docs/PSFEstimator-${PSFEstimator_VERSION_STRING}-UserGuide.tex"
  COPYONLY )
