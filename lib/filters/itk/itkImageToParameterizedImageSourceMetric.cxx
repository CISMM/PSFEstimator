/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToParameterizedImageSourceMetric.cxx,v $
  Language:  C++
  Date:      $Date: 2009/09/17 20:30:15 $
  Version:   $Revision: 1.5 $

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
  m_ParametersMask        = ParametersMaskType(0);
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
 * Set the moving image source.
 */
template <class TFixedImage, class TMovingImageSource>
void
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::SetMovingImageSource(MovingImageSourceType* source) {
  itkDebugMacro("setting MovingImageSource to " << source );
  if (this->m_MovingImageSource != source) {
    this->m_MovingImageSource = source;
    this->Modified();
    
    // Now reinitialize the parameter mask array to match the number
    // of the parameters in the moving image source.
    m_ParametersMask = ParametersMaskType(source->GetNumberOfParameters());
    
    // Initialize to have no active parameters.
    for (unsigned int i = 0; i < m_ParametersMask.Size(); i++) {
      m_ParametersMask[i] = 0;
    }
  }
}


/**
 * Set the ImageToImageMetric.
 */
template <class TFixedImage, class TMovingImageSource>
void
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::SetImageToImageMetric(ImageToImageMetricType* source) {
  itkDebugMacro("setting ImageToImageMetric to " << source );
  if (this->m_ImageToImageMetric != source) {
    this->m_ImageToImageMetric = source;
    this->Modified();
    
    m_ImageToImageMetric->SetTransform(m_Transform);
    m_ImageToImageMetric->SetInterpolator(m_Interpolator);
  }
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
 * Get the value of the cost function. The parameters passed into this method
 * are assumed to be the active parameters which is a subset of all parameters.
 */
template <class TFixedImage, class TMovingImageSource>
typename ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>::MeasureType
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::GetValue(const ParametersType& parameters) const {
  // Send the parameters to the parametric image source.
  //m_MovingImageSource->SetParameters(parameters);
  SetParameters(parameters);
  std::cout << "Parameters: " << parameters << " - ";
  m_MovingImageSource->GetOutput()->SetRequestedRegionToLargestPossibleRegion();
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

  MeasureType value = m_ImageToImageMetric->GetValue(parameters);
  std::cout << "Value: " << value << std::endl;
  return value;
}


/**
 * Set the active parameters of the moving image source.
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

  // Iterate through the parameters mask and set only the active parameters
  ParametersType allParameters = m_MovingImageSource->GetParameters();
  int activeIndex = 0;
  for (unsigned int i = 0; i < allParameters.Size(); i++) {
    if (m_ParametersMask[i])
      allParameters[i] = parameters[activeIndex++];
  }

  m_MovingImageSource->SetParameters( allParameters );
}


/**
 * Get the number of parameters.
 */
template <class TFixedImage, class TMovingImageSource>
unsigned int
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::GetNumberOfParameters(void) const {
  ParametersMaskType* mask = const_cast< ImageToParameterizedImageSourceMetric< TFixedImage,TMovingImageSource >* >(this)->GetParametersMask();

  // Count up the active parameters
  unsigned int activeCount = 0;
  for (unsigned int i = 0; i < m_MovingImageSource->GetNumberOfParameters(); i++) {
    if (mask->GetElement(i))
      activeCount++;
  }
  
  return activeCount;
}


/**
 * Get parameters mask array.
 */
template <class TFixedImage, class TMovingImageSource> 
typename ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>::ParametersMaskType*
ImageToParameterizedImageSourceMetric<TFixedImage,TMovingImageSource>
::GetParametersMask() throw ( ExceptionObject ) {
  if (m_ParametersMask.Size() == 0) {
    itkExceptionMacro(<<"MovingImageSource is not present so the parameters mask has not been initialized");
  }

  return &m_ParametersMask;
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
}


} // end namespace itk


#endif
