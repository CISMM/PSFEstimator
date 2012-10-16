SET(modelSrc
  DataModel.cxx
  StringUtils.cxx
)

ADD_LIBRARY(psfeModel ${modelSrc})

TARGET_LINK_LIBRARIES( psfeModel
  ${VTK_LIBRARIES}
  ${ITK_LIBRARIES}
  ITKMicroscopyPSFToolkit
)
