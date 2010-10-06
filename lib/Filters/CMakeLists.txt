SET(filterSrc
  PassThroughFilter.cxx
  itkGibsonLanniBSFImageSource.txx
  itkGibsonLanniPSFImageSource.txx
  itkImageToParametricImageSourceMetric.txx
  itkParametricImageSource.txx
  itkPoissonNoiseImageToImageMetric.txx
  itkScanImageFilter.txx
  itkSphereConvolutionFilter.txx
)

ADD_LIBRARY(psfeFilter ${filterSrc})
