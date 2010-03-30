SET(filterSrc
  PassThroughFilter.cxx
  itk/itkGibsonLanniBSFImageSource.cxx
  itk/itkGibsonLanniPSFImageSource.cxx
  itk/itkImageToParameterizedImageSourceMetric.cxx
  itk/itkParameterizedImageSource.cxx
  itk/itkPoissonNoiseImageToImageMetric.cxx
  itk/itkScanImageFilter.cxx
  itk/itkSphereConvolutionFilter.cxx
)

ADD_LIBRARY(vpsfoFilter ${filterSrc})
