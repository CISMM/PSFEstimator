SET(filterSrc
  PassThroughFilter.cxx
  itk/itkGibsonLanniBSFImageSource.txx
  itk/itkGibsonLanniPSFImageSource.txx
  itk/itkImageToParametricImageSourceMetric.txx
  itk/itkParametricImageSource.txx
  itk/itkPoissonNoiseImageToImageMetric.txx
  itk/itkScanImageFilter.txx
  itk/itkSphereConvolutionFilter.txx
)

ADD_LIBRARY(psfeFilter ${filterSrc})
