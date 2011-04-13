ADD_EXECUTABLE( BatchPSFEstimator BatchPSFEstimator.cxx )

TARGET_LINK_LIBRARIES( BatchPSFEstimator
  ${ITK_LIBRARIES}
  psfeModel
  psfeIO
)
