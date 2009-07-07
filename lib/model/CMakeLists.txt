SET(modelSrc
  DataModel.cxx
)

ADD_LIBRARY(vpsfoModel ${modelSrc})

TARGET_LINK_LIBRARIES( vpsfoModel
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
)
