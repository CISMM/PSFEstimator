SET(filterSrc
  PassThroughFilter.cxx
  itk/itkGibsonLanniPSFImageSource.cxx
  itk/itkParameterizedImageSource.cxx
  itk/itkImageToParameterizedImageSourceMetric.cxx
)

ADD_LIBRARY(vpsfoFilter ${filterSrc})
