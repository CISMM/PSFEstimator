#ifndef __itkMicroscopePointSpreadFunctionImageSource_txx
#define __itkMicroscopePointSpreadFunctionImageSource_txx

#include "itkMicroscopePointSpreadFunctionImageSource.h"


namespace itk
{


template< class TOutputImage >
MicroscopePointSpreadFunctionImageSource< TOutputImage >
::MicroscopePointSpreadFunctionImageSource()
{
  this->m_Magnification = 60.0;
  this->m_NumericalAperture = 1.4;
}


template< class TOutputImage >
MicroscopePointSpreadFunctionImageSource< TOutputImage >
::~MicroscopePointSpreadFunctionImageSource()
{
}


template< class TOutputImage >
void
MicroscopePointSpreadFunctionImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Magnification: " << m_Magnification << std::endl;
  os << indent << "NumericalAperture: " << m_NumericalAperture << std::endl;
}

} // end namespace itk

#endif // __itkMicroscopePointSpreadFunctionImageSource_txx
