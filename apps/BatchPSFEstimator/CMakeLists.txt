ADD_EXECUTABLE( BatchPSFEstimator BatchPSFEstimator.cxx )

TARGET_LINK_LIBRARIES( BatchPSFEstimator
  itksys
  ITKCommon
  ITKBasicFilters
  ITKIO
  ITKNumerics
  psfeModel
  psfeIO
)
