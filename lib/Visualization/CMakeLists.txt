SET(visualizationSrc
  ImagePlaneVisualizationPipeline.cxx
  IsosurfaceVisualizationPipeline.cxx
  OutlineVisualizationPipeline.cxx
  VisualizationPipeline.cxx
  Visualization.cxx
)

ADD_LIBRARY(psfeVisualization ${visualizationSrc})
