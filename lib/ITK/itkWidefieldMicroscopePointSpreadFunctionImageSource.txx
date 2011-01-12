#ifndef __itkWidefieldMicroscopePointSpreadFunctionImageSource_txx
#define __itkWidefieldMicroscopePointSpreadFunctionImageSource_txx

#include "itkWidefieldMicroscopePointSpreadFunctionImageSource.h"


namespace itk
{


template< class TOutputImage >
WidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::WidefieldMicroscopePointSpreadFunctionImageSource()
{
  this->m_EmissionWavelength = 550.0;
}


template< class TOutputImage >
WidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::~WidefieldMicroscopePointSpreadFunctionImageSource()
{
}


template< class TOutputImage >
void
WidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  std::cout << indent << "EmissionWavelength: " << m_EmissionWavelength << std::endl;
}


} // end namespace itk

#endif // __itkWidefieldMicroscopePointSpreadFunctionImageSource_txx
