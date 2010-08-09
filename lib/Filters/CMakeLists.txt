SET(filterSrc
  PassThroughFilter.cxx
  itk/itkGibsonLanniBSFImageSource.cxx
  itk/itkGibsonLanniPSFImageSource.cxx
  itk/itkImageToParametricImageSourceMetric.txx
  itk/itkParametricImageSource.txx
  itk/itkPoissonNoiseImageToImageMetric.cxx
  itk/itkScanImageFilter.cxx
  itk/itkSphereConvolutionFilter.cxx
)

ADD_LIBRARY(vpsfoFilter ${filterSrc})
