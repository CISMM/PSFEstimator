/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToParametricImageSourceMetric.txx,v $
  Language:  C++
  Date:      $Date: 2009/09/17 20:30:15 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToParametricImageSourceMetric_txx
#define __itkImageToParametricImageSourceMetric_txx

// First, make sure that we include the configuration file.
// This line may be removed once the ThreadSafeTransform gets
// integrated into ITK.
#include "itkConfigure.h"

#include "itkImageToParametricImageSourceMetric.h"

namespace itk
{

/**
 * Constructor
 */
template <class TFixedImage, class TMovingImageSource>
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::ImageToParametricImageSourceMetric()
{
  m_FixedImage        = 0; // has to be provided by the user.
  m_MovingImageSource = 0; // has to be provided by the user.
  m_DelegateMetric    = 0; // has to be provided by the user.
  m_Transform         = TransformType::New(); // immutable
  m_Interpolator      = 0; // has to be provided by the user.
  m_ParametersMask    = ParametersMaskType(0);
}


/**
 * Destructor
 */
template <class TFixedImage, class TMovingImageSource>
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::~ImageToParametricImageSourceMetric()
{

}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetMovingImageSource(MovingImageSourceType* source)
{
  itkDebugMacro("setting MovingImageSource to " << source );
  if (m_MovingImageSource != source ||
      m_MovingImageSource->GetNumberOfParameters() != m_ParametersMask.Size())
    {
    this->m_MovingImageSource = source;
    this->Modified();

    // Now reinitialize the parameter mask array to match the number
    // of the parameters in the moving image source.
    m_ParametersMask = ParametersMaskType(source->GetNumberOfParameters());

    // Initialize to have no active parameters.
    for ( unsigned int i = 0; i < m_ParametersMask.Size(); i++ )
      {
      m_ParametersMask[i] = 0;
      }
    }
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetTransform(TransformType* transform)
{
  itkWarningMacro(<< "Setting the Transform on an "
                  << "ImageToParametrizedImageSourceMetric has no effect. "
                  << "Default IdentityTransform still in use.");
}


template <class TFixedImage, class TMovingImageSource>
const unsigned long &
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetNumberOfPixelsCounted() const
{
  if ( m_DelegateMetric )
    {
    return m_DelegateMetric->GetNumberOfPixelsCounted();
    }
  return 0L;
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetFixedImageRegion(FixedImageRegionType region) {
  if ( m_DelegateMetric )
    {
    m_DelegateMetric->SetFixedImageRegion(region);
    }
}


template <class TFixedImage, class TMovingImageSource>
const typename ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>::FixedImageRegionType & 
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetFixedImageRegion()
{
  return m_DelegateMetric->GetFixedImageRegion();
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetDelegateMetric(DelegateMetricType* source)
{
  itkDebugMacro("setting DelegateMetric to " << source );
  if ( this->m_DelegateMetric != source )
    {
    this->m_DelegateMetric = source;
    this->Modified();

    m_DelegateMetric->SetTransform(m_Transform);
    m_DelegateMetric->SetInterpolator(m_Interpolator);
    }
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetDerivative(const ParametersType& parameters, DerivativeType& derivative) const
{
  // Do nothing.
}


template <class TFixedImage, class TMovingImageSource>
typename ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>::MeasureType
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetValue(const ParametersType& parameters) const
{
  // Send the parameters to the parametric image source.
  std::cout << "Parameters: " << parameters << std::endl;
  SetParameters(parameters);

  // Now update the parametric image source.
  m_MovingImageSource->GetOutput()->SetRequestedRegionToLargestPossibleRegion();
  m_MovingImageSource->Update();

  m_DelegateMetric->SetFixedImage(m_FixedImage);
  m_DelegateMetric->
    SetFixedImageRegion(m_FixedImage->GetLargestPossibleRegion());

  // Have to set the new moving image in the interpolator manually because
  // the delegate image to image metric does this only at initialization.
  MovingImageSourceOutputImagePointerType movingImage =
    m_MovingImageSource->GetOutput();
  m_Interpolator->SetInputImage(movingImage);

  // Now we can set the moving image in the image to image metric.
  m_DelegateMetric->SetMovingImage(movingImage);

  // We have to initialize the delegate metric here to avoid an exception
  m_DelegateMetric->Initialize();
  MeasureType value = m_DelegateMetric->GetValue(parameters);
  std::cout << "Value: " << value << std::endl;

  return value;
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::SetParameters( const ParametersType & parameters ) const
{
  if( !m_MovingImageSource )
    {
    itkExceptionMacro(<<"Moving image source has not been assigned");
    }

  // Iterate through the parameters mask and set only the active parameters
  ParametersType allParameters = m_MovingImageSource->GetParameters();
  int activeIndex = 0;
  for (unsigned int i = 0; i < allParameters.Size(); i++)
    {
    if ( m_ParametersMask[i] )
      {
      allParameters[i] = parameters[activeIndex++];
      }
    }

  m_MovingImageSource->SetParameters( allParameters );
}


template <class TFixedImage, class TMovingImageSource>
unsigned int
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetNumberOfParameters(void) const
{
  ParametersMaskType* mask = const_cast< ImageToParametricImageSourceMetric< TFixedImage,TMovingImageSource >* >(this)->GetParametersMask();

  // Count up the active parameters
  unsigned int activeCount = 0;
  for (unsigned int i = 0; i < m_MovingImageSource->GetNumberOfParameters(); i++)
    {
    if ( mask->GetElement(i) )
      {
      activeCount++;
      }
    }

  return activeCount;
}


template <class TFixedImage, class TMovingImageSource>
typename ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>::ParametersMaskType*
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::GetParametersMask() throw ( ExceptionObject )
{
  if ( m_ParametersMask.Size() == 0 )
    {
    itkExceptionMacro(<<"MovingImageSource is not present so the parameters mask has not been initialized");
    }

  return &m_ParametersMask;
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
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

  if( !m_DelegateMetric )
    {
    itkExceptionMacro(<<"ImageToImageMetric is not present");
    }

  // If there are any observers on the metric, call them to give the
  // user code a chance to set parameters on the metric
  this->InvokeEvent( InitializeEvent() );
}


template <class TFixedImage, class TMovingImageSource>
void
ImageToParametricImageSourceMetric<TFixedImage,TMovingImageSource>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Fixed  Image: " << m_FixedImage.GetPointer() << std::endl;
  os << indent << "Moving Image Source: " << m_MovingImageSource.GetPointer() << std::endl;
  os << indent << "Image to Image Metric: " << m_DelegateMetric.GetPointer() << std::endl;
  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "ParametersMask: " << m_ParametersMask << std::endl;
}

} // end namespace itk


#endif
