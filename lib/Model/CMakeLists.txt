SET(modelSrc
  DataModel.cxx
)

ADD_LIBRARY(psfeModel ${modelSrc})

TARGET_LINK_LIBRARIES( psfeModel
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
)
