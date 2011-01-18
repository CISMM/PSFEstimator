/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniPSFImageSource.cxx,v $
  Language:  C++
  Date:      $Date: 2010/05/17 15:41:35 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGibsonLanniPointSpreadFunctionImageSource_txx
#define __itkGibsonLanniPointSpreadFunctionImageSource_txx

#include "itkGibsonLanniPointSpreadFunctionImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMath.h"
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"

namespace itk
{

//----------------------------------------------------------------------------
template< class TOutputImage >
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::GibsonLanniPointSpreadFunctionImageSource()
{
  this->m_Size.Fill(32);
  this->m_Spacing.Fill(65.0);
  this->m_Origin.Fill(0.0);

  // Set default PSF model parameters.
  this->m_EmissionWavelength   = 550.0; // in nanometers
  this->m_NumericalAperture    = 1.4;   // unitless
  this->m_Magnification        = 60.0;  // unitless

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


//----------------------------------------------------------------------------
template< class TOutputImage >
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::~GibsonLanniPointSpreadFunctionImageSource()
{
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::SetParameter(unsigned int index, ParametersValueType value)
{
  switch (index)
    {
    case 0:
      this->SetEmissionWavelength(value);
      break;

    case 1:
      this->SetNumericalAperture(value);
      break;

    case 2:
      this->SetMagnification(value);
      break;

    case 3:
      this->SetDesignCoverSlipRefractiveIndex(value);
      break;

    case 4:
      this->SetActualCoverSlipRefractiveIndex(value);
      break;

    case 5:
      this->SetDesignCoverSlipThickness(value);
      break;

    case 6:
      this->SetActualCoverSlipThickness(value);
      break;

    case 7:
      this->SetDesignImmersionOilRefractiveIndex(value);
      break;

    case 8:
      this->SetActualImmersionOilRefractiveIndex(value);
      break;

    case 9:
      this->SetDesignImmersionOilThickness(value);
      break;

    case 10:
      this->SetDesignSpecimenLayerRefractiveIndex(value);
      break;

    case 11:
      this->SetActualSpecimenLayerRefractiveIndex(value);
      break;

    case 12:
      this->SetActualPointSourceDepthInSpecimenLayer(value);
      break;
    }
}


//----------------------------------------------------------------------------
template< class TOutputImage >
typename GibsonLanniPointSpreadFunctionImageSource<TOutputImage>::ParametersValueType
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::GetParameter(unsigned int index) const
{
  switch (index)
    {
    case 0:
      return this->GetEmissionWavelength();
      break;

    case 1:
      return this->GetNumericalAperture();
      break;

    case 2:
      return this->GetMagnification();
      break;

    case 3:
      return this->GetDesignCoverSlipRefractiveIndex();
      break;

    case 4:
      return this->GetActualCoverSlipRefractiveIndex();
      break;

    case 5:
      return this->GetDesignCoverSlipThickness();
      break;

    case 6:
      return this->GetActualCoverSlipThickness();
      break;

    case 7:
      return this->GetDesignImmersionOilRefractiveIndex();
      break;

    case 8:
      return this->GetActualImmersionOilRefractiveIndex();
      break;

    case 9:
      return this->GetDesignImmersionOilThickness();
      break;

    case 10:
      return this->GetDesignSpecimenLayerRefractiveIndex();
      break;

    case 11:
      return this->GetActualSpecimenLayerRefractiveIndex();
      break;

    case 12:
      return this->GetActualPointSourceDepthInSpecimenLayer();
      break;

    default:
      return 0.0;
      break;
    }
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::SetParameters(const ParametersType& parameters)
{
  int index = 0;

  SetEmissionWavelength(parameters[index++]);
  SetNumericalAperture(parameters[index++]);
  SetMagnification(parameters[index++]);

  SetDesignCoverSlipRefractiveIndex(parameters[index++]);
  SetActualCoverSlipRefractiveIndex(parameters[index++]);
  SetDesignCoverSlipThickness(parameters[index++]);
  SetActualCoverSlipThickness(parameters[index++]);
  SetDesignImmersionOilRefractiveIndex(parameters[index++]);
  SetActualImmersionOilRefractiveIndex(parameters[index++]);
  SetDesignImmersionOilThickness(parameters[index++]);

  SetDesignSpecimenLayerRefractiveIndex(parameters[index++]);
  SetActualSpecimenLayerRefractiveIndex(parameters[index++]);
  SetActualPointSourceDepthInSpecimenLayer(parameters[index++]);
}


template< class TOutputImage >
typename GibsonLanniPointSpreadFunctionImageSource<TOutputImage>::ParametersType
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::GetParameters() const
{
  ParametersType parameters(GetNumberOfParameters());

  int index = 0;
  parameters[index++] = this->GetEmissionWavelength();
  parameters[index++] = this->GetNumericalAperture();
  parameters[index++] = this->GetMagnification();

  parameters[index++] = this->GetDesignCoverSlipRefractiveIndex();
  parameters[index++] = this->GetActualCoverSlipRefractiveIndex();
  parameters[index++] = this->GetDesignCoverSlipThickness();
  parameters[index++] = this->GetActualCoverSlipThickness();
  parameters[index++] = this->GetDesignImmersionOilRefractiveIndex();
  parameters[index++] = this->GetActualImmersionOilRefractiveIndex();
  parameters[index++] = this->GetDesignImmersionOilThickness();

  parameters[index++] = this->GetDesignSpecimenLayerRefractiveIndex();
  parameters[index++] = this->GetActualSpecimenLayerRefractiveIndex();
  parameters[index++] = this->GetActualPointSourceDepthInSpecimenLayer();

  return parameters;
}


//----------------------------------------------------------------------------
template< class TOutputImage >
unsigned int
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return 13;
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::BeforeThreadedGenerateData()
{
  this->m_IntegrandFunctor.CopySettings(this);
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::ThreadedGenerateData(const RegionType& outputRegionForThread, int threadId )
{
  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  typename TOutputImage::Pointer image = this->GetOutput(0);

  ImageRegionIteratorWithIndex<OutputImageType> it(image, outputRegionForThread);

  for (; !it.IsAtEnd(); ++it)
    {
    IndexType index = it.GetIndex();
    PointType point;
    image->TransformIndexToPhysicalPoint(index, point);

    it.Set( ComputeSampleValue( point ));
    progress.CompletedPixel();
    }
}


//----------------------------------------------------------------------------
template< class TOutputImage >
double
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::ComputeSampleValue(typename TOutputImage::PointType& point)
{
  PixelType px = point[0] * 1e-9;
  PixelType py = point[1] * 1e-9;
  PixelType pz = point[2] * 1e-9;

  double mag = this->m_Magnification;

  // We have to convert to coordinates of the detector points
  double x_o = px * mag;
  double y_o = py * mag;
  double z_o = pz; // No conversion needed

  // Return squared magnitude of the integrated value
  return static_cast<PixelType>(
    norm( Integrate( this->m_IntegrandFunctor, 0.0, 1.0, 20, x_o, y_o, z_o)));
}

} // end namespace itk

#endif
