ADD_EXECUTABLE( BatchPSFOptimizer BatchPSFOptimizer.cxx )

TARGET_LINK_LIBRARIES( BatchPSFOptimizer
  itksys
  ITKCommon
  ITKBasicFilters
  ITKIO
  ITKNumerics
  vpsfoModel
  vpsfoFilter
  vpsfoIO
)
