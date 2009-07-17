/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageToParameterizedImageSourceMetric.h,v $
  Language:  C++
  Date:      $Date: 2009/07/17 16:10:18 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToParameterizedImageSourceMetric_h
#define __itkImageToParameterizedImageSourceMetric_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "itkCostFunction.h"
#include "itkExceptionObject.h"
#include "itkIdentityTransform.h"
#include "itkImageBase.h"
#include "itkImageToImageMetric.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkSingleValuedCostFunction.h"

namespace itk
{
  
/** \class ImageToParameterizedImageSourceMetric
 * \brief Computes similarity between two images, one of which is fixed and
 * the other generated from a parameterized image source.
 *
 * This class computes a value that measures the similarity 
 * between the Fixed image and the parameterized Moving image. "Moving"
 * in this metric and subclasses refers to changes in the parameters used
 * to generate the image, not movement induced by a spatial transformation.
 *
 * This class is parameterized over two types. The first template class
 * is the fixed image data and the second template class is the source of
 * the moving image source, typically a subclass of 
 * itkParameterizedImageSource.
 *
 * \ingroup RegistrationMetrics
 *
 */

template <class TFixedImage,  class TMovingImageSource> 
class ITK_EXPORT ImageToParameterizedImageSourceMetric : public SingleValuedCostFunction 
{
public:
  /** Standard class typedefs. */
  typedef ImageToParameterizedImageSourceMetric           Self;
  typedef SingleValuedCostFunction     Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToParameterizedImageSourceMetric, SingleValuedCostFunction);

  /**  Type of the moving Image. */
  typedef TMovingImageSource                         MovingImageSourceType;
  typedef typename TMovingImageSource::PixelType     MovingImageSourcePixelType;
  typedef typename MovingImageSourceType::Pointer    MovingImageSourcePointer;
  typedef typename MovingImageSourceType::OutputImageType 
    MovingImageSourceOutputImageType;
  typedef typename MovingImageSourceOutputImageType::Pointer
    MovingImageSourceOutputImagePointerType;

  /**  Type of the fixed Image. */
  typedef TFixedImage                                FixedImageType;
  typedef typename FixedImageType::ConstPointer      FixedImageConstPointer;
  typedef typename FixedImageType::RegionType        FixedImageRegionType;

  /**  Type of the delegate image comparison metric. */
  typedef ImageToImageMetric<FixedImageType, MovingImageSourceOutputImageType> ImageToImageMetricType;
  typedef typename ImageToImageMetricType::Pointer ImageToImageMetricTypePointer;

  /** Constants for the image dimensions */
  itkStaticConstMacro(MovingImageSourceDimension, 
                      unsigned int,
                      TMovingImageSource::ImageDimension);
  itkStaticConstMacro(FixedImageDimension, 
                      unsigned int,
                      TFixedImage::ImageDimension);
  
  typedef typename NumericTraits<MovingImageSourcePixelType>::RealType RealType;

  /**  Type of the measure. */
  typedef typename Superclass::MeasureType         MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType      DerivativeType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::ParametersType      ParametersType;

  /** Transform and interpolator to pass to the delegate image to image metric. */
  typedef IdentityTransform<double>        TransformType;
  typedef typename TransformType::Pointer  TransformTypePointer;

  typedef NearestNeighborInterpolateImageFunction<FixedImageType, double>
    InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorTypePointer;

  /** Connect the Fixed Image.  */
  itkSetConstObjectMacro( FixedImage, FixedImageType );

  /** Get the Fixed Image. */
  itkGetConstObjectMacro( FixedImage, FixedImageType );

  /** Connect the Moving Image Source.  */
  itkSetObjectMacro( MovingImageSource, MovingImageSourceType );

  /** Get the Moving Image Source. */
  itkGetObjectMacro( MovingImageSource, MovingImageSourceType );

  /** Get the number of pixels considered in the computation. */
  itkGetConstReferenceMacro( NumberOfPixelsCounted, unsigned long );

  /** Set the region over which the metric will be computed */
  itkSetMacro( FixedImageRegion, FixedImageRegionType );

  /** Get the region over which the metric will be computed */
  itkGetConstReferenceMacro( FixedImageRegion, FixedImageRegionType );

  /** Set the delegate ImageToImageMetric. */
  //itkSetObjectMacro( ImageToImageMetric, ImageToImageMetricType );
  virtual void SetImageToImageMetric(ImageToImageMetricType* _arg) {
    itkDebugMacro("setting ImageToImageMetric to " << _arg );
    if (this->m_ImageToImageMetric != _arg) {
	this->m_ImageToImageMetric = _arg;
	this->Modified();

	m_ImageToImageMetric->SetTransform(m_Transform);
	m_ImageToImageMetric->SetInterpolator(m_Interpolator);
    }
  }

  /** Get the delegate ImageToImageMetric. */
  itkGetConstObjectMacro( ImageToImageMetric, ImageToImageMetricType );

  /** Get the derivative of the cost function (calling this returns an
      undefined derivative); */
  virtual void GetDerivative(const ParametersType& parameters, DerivativeType& derivative) const;

  /** Get the value of the cost function. */
  virtual MeasureType GetValue(const ParametersType& parameters) const;
 
  /** Set the parameters for the Moving Image Source. */
  void SetParameters( const ParametersType & parameters ) const;

  /** Return the number of parameters required by the Transform */
  virtual unsigned int GetNumberOfParameters(void) const 
    { return m_MovingImageSource->GetNumberOfParameters(); }

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject );

protected:
  ImageToParameterizedImageSourceMetric();
  virtual ~ImageToParameterizedImageSourceMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

  mutable unsigned long         m_NumberOfPixelsCounted;

  FixedImageConstPointer        m_FixedImage;
  MovingImageSourcePointer      m_MovingImageSource;
  
  ImageToImageMetricTypePointer m_ImageToImageMetric;

  /* Disable spatial registration with an identity transform. */
  TransformTypePointer    m_Transform;

  /* Use nearest neighbor interpolation because we will just be looking up
     pixel values at grid positions. */
  InterpolatorTypePointer m_Interpolator;

private:
  ImageToParameterizedImageSourceMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  FixedImageRegionType        m_FixedImageRegion;  

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToParameterizedImageSourceMetric.cxx"
#endif

#endif
