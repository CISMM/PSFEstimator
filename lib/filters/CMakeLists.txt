SET(filterSrc
  PassThroughFilter.cxx
  itk/itkGibsonLanniPSFImageSource.cxx
)

ADD_LIBRARY(vpsfoFilter ${filterSrc})
