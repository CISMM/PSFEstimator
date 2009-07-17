/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToParameterizedImageSourceMetric.cxx,v $
  Language:  C++
  Date:      $Date: 2009/07/17 16:10:19 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToParameterizedImageSourceMetric_txx
#define __itkImageToParameterizedImageSourceMetric_txx

// First, make sure that we include the configuration file.
// This line may be removed once the ThreadSafeTransform gets
// integrated into ITK.
#include "itkConfigure.h"

#include "itkImageToParameterizedImageSourceMetric.h"

namespace itk
{

/**
 * Constructor
 */
template <class TFixedImage, class TMovingImageSource> 
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::ImageToParameterizedImageSourceMetric()
{
  m_FixedImage            = 0; // has to be provided by the user.
  m_MovingImageSource     = 0; // has to be provided by the user.
  m_ImageToImageMetric    = 0; // has to be provided by the user.
  m_Transform             = TransformType::New();
  m_Interpolator          = InterpolatorType::New();
  m_NumberOfPixelsCounted = 0; // initialize to zero
}

/**
 * Destructor
 */
template <class TFixedImage, class TMovingImageSource> 
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::~ImageToParameterizedImageSourceMetric()
{

}


/**
 * Get the derivative of the cost function.
 */
template <class TFixedImage, class TMovingImageSource>
void
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::GetDerivative(const ParametersType& parameters, DerivativeType& derivative) const {
  // Do nothing.
}


/**
 * Get the value of the cost function. 
 */
template <class TFixedImage, class TMovingImageSource>
typename ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>::MeasureType
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::GetValue(const ParametersType& parameters) const {
  // Send the parameters to the parametric image source.
  m_MovingImageSource->SetParameters(parameters);
  m_MovingImageSource->Update();

  m_ImageToImageMetric->SetFixedImage(m_FixedImage);
  m_ImageToImageMetric->SetFixedImageRegion(m_FixedImage->
					    GetLargestPossibleRegion());

  // Have to set the new moving image in the interpolator manually because
  // the delegate image to image metric does this only at initialization.
  MovingImageSourceOutputImagePointerType movingImage = 
    m_MovingImageSource->GetOutput();
  m_Interpolator->SetInputImage(movingImage);

  // Now we can set the moving image in the image to image metric.
  m_ImageToImageMetric->SetMovingImage(movingImage);

  return m_ImageToImageMetric->GetValue(parameters);
}


/**
 * Set the parameters of the moving image source.
 */
template <class TFixedImage, class TMovingImageSource> 
void
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::SetParameters( const ParametersType & parameters ) const
{
  if( !m_MovingImageSource )
    {
    itkExceptionMacro(<<"Moving image source has not been assigned");
    }
  m_MovingImageSource->SetParameters( parameters );
}


/**
 * Initialize
 */
template <class TFixedImage, class TMovingImageSource> 
void
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::Initialize(void) throw ( ExceptionObject )
{

  if( !m_MovingImageSource )
    {
    itkExceptionMacro(<<"MovingImageSource is not present");
    }

  if( !m_FixedImage )
    {
    itkExceptionMacro(<<"FixedImage is not present");
    }

  if( !m_ImageToImageMetric )
    {
    itkExceptionMacro(<<"ImageToImageMetric is not present");
    }

  // Make sure the FixedImageRegion is within the FixedImage buffered region
  if ( !m_FixedImageRegion.Crop( m_FixedImage->GetBufferedRegion() ) )
    {
    itkExceptionMacro(
      <<"FixedImageRegion does not overlap the fixed image buffered region" );
    }

  // If there are any observers on the metric, call them to give the
  // user code a chance to set parameters on the metric
  this->InvokeEvent( InitializeEvent() );
}


/**
 * PrintSelf
 */
template <class TFixedImage, class TMovingImageSource> 
void
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Moving Image Source: " << m_MovingImageSource.GetPointer()  << std::endl;
  os << indent << "Fixed  Image: " << m_FixedImage.GetPointer()   << std::endl;
  os << indent << "FixedImageRegion: " << m_FixedImageRegion << std::endl;
  os << indent << "Number of Pixels Counted: " << m_NumberOfPixelsCounted 
     << std::endl;

}


} // end namespace itk


#endif
