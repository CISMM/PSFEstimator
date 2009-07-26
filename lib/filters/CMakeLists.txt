SET(filterSrc
  PassThroughFilter.cxx
  itk/itkGibsonLanniPSFImageSource.cxx
  itk/itkImageToParameterizedImageSourceMetric.cxx
  itk/itkParameterizedImageSource.cxx
  itk/itkScanImageFilter.cxx
)

ADD_LIBRARY(vpsfoFilter ${filterSrc})
