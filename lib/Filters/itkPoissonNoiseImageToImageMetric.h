/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPoissonNoiseImageToImageMetric.h,v $
  Language:  C++
  Date:      $Date: 2010/03/30 04:17:56 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPoissonNoiseImageToImageMetric_h
#define __itkPoissonNoiseImageToImageMetric_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "itkImageToImageMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"


namespace itk
{
/** \class PoissonNoiseImageToImageMetric
 * \brief Computes similarity between two objects to be registered
 *
 * This Class is templated over the type of the fixed and moving
 * images to be compared.
 *
 * This metrics computes the negative log likelihood of the fixed image by
 * treating voxels in the moving image as the mean (and standard deviation)
 * of a Poisson distribution and calculating the probability of the intensity
 * of the corresponding voxel in the fixed image according to that
 * distribution. Each voxel distribution is treated as independent.
 *
 * \ingroup RegistrationMetrics
 */
template < class TFixedImage, class TMovingImage >
class ITK_EXPORT PoissonNoiseImageToImageMetric :
    public ImageToImageMetric< TFixedImage, TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef PoissonNoiseImageToImageMetric                 Self;
  typedef ImageToImageMetric<TFixedImage, TMovingImage > Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PoissonNoiseImageToImageMetric, ImageToImageMetric);


  /** Types transferred from the base class */
  typedef typename Superclass::RealType                 RealType;
  typedef typename Superclass::TransformType            TransformType;
  typedef typename Superclass::TransformPointer         TransformPointer;
  typedef typename Superclass::TransformParametersType  TransformParametersType;
  typedef typename Superclass::TransformJacobianType    TransformJacobianType;
  typedef typename Superclass::GradientPixelType        GradientPixelType;
  typedef typename Superclass::GradientImageType        GradientImageType;
  typedef typename Superclass::InputPointType           InputPointType;
  typedef typename Superclass::OutputPointType          OutputPointType;

  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::FixedImageConstPointer   FixedImageConstPointer;
  typedef typename Superclass::MovingImageConstPointer  MovingImageConstPointer;


  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
                      DerivativeType & derivative ) const;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
                              MeasureType& Value, DerivativeType& Derivative ) const;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(MovingPixelTypeHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TMovingImage::PixelType>));
  itkConceptMacro(MovingRealTypeAdditivieOperatorsCheck,
    (Concept::AdditiveOperators<
     typename NumericTraits<typename TMovingImage::PixelType>::RealType>));
  itkConceptMacro(MovingRealTypeMultiplyOperatorCheck,
    (Concept::MultiplyOperator<
     typename NumericTraits<typename TMovingImage::PixelType>::RealType>));

  /** End concept checking */
#endif

protected:
  PoissonNoiseImageToImageMetric();
  virtual ~PoissonNoiseImageToImageMetric() {};

private:
  PoissonNoiseImageToImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPoissonNoiseImageToImageMetric.txx"
#endif

#endif
