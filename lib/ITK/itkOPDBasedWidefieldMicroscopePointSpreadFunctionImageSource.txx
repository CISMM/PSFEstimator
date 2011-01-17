#ifndef __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_txx
#define __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_txx

#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource.h"


namespace itk
{


template< class TOutputImage >
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource()
{
  this->m_DesignCoverSlipRefractiveIndex    = 1.522; // unitless
  this->m_ActualCoverSlipRefractiveIndex    = 1.522; // unitless
  this->m_DesignCoverSlipThickness          = 170.0; // in micrometers
  this->m_ActualCoverSlipThickness          = 170.0; // in micrometers
  this->m_DesignImmersionOilRefractiveIndex = 1.515; // unitless
  this->m_ActualImmersionOilRefractiveIndex = 1.515; // unitless
  this->m_DesignImmersionOilThickness       = 100.0; // in micrometers

  this->m_DesignSpecimenLayerRefractiveIndex         =  1.33; // unitless
  this->m_ActualSpecimenLayerRefractiveIndex         =  1.33; // unitless
  this->m_ActualPointSourceDepthInSpecimenLayer      =   0.0; // in micrometers
}


template< class TOutputImage >
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::~OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource()
{
}


template< class TOutputImage >
void
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  std::cout << indent << "DesignCoverSlipRefractiveIndex: "
            << m_DesignCoverSlipRefractiveIndex << std::endl;
  std::cout << indent << "ActualCoverSlipRefractiveIndex: "
            << m_ActualCoverSlipRefractiveIndex << std::endl;
  std::cout << indent << "DesignCoverSlipThickness: "
            << m_DesignCoverSlipThickness << std::endl;
  std::cout << indent << "ActualCoverSlipThickness: "
            << m_ActualCoverSlipThickness << std::endl;
  std::cout << indent << "DesignImmersionOilRefractiveIndex: "
            << m_DesignImmersionOilRefractiveIndex << std::endl;
  std::cout << indent << "ActualImmersionOilRefractiveIndex: "
            << m_ActualImmersionOilRefractiveIndex << std::endl;
  std::cout << indent << "DesignImmersionOilThickness: "
            << m_DesignImmersionOilThickness << std::endl;
  std::cout << indent << "DesignSpecimenLayerRefractiveIndex: "
            << m_DesignSpecimenLayerRefractiveIndex << std::endl;
  std::cout << indent << "ActualSpecimenLayerRefractiveIndex: "
            << m_ActualSpecimenLayerRefractiveIndex << std::endl;
  std::cout << indent << "ActualPointSourceDepthInSpecimenLayer: "
            << m_ActualPointSourceDepthInSpecimenLayer << std::endl;
}

} // end namespace itk


#endif // __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_txx
