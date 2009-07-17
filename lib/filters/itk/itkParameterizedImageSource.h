/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParameterizedImageSource.h,v $
  Language:  C++
  Date:      $Date: 2009/07/17 16:10:19 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkParameterizedImageSource_h
#define __itkParameterizedImageSource_h

#include "itkArray.h"
#include "itkImageSource.h"

namespace itk
{

/** \class ParameterizedParameterizedImageSource
 *  \brief Extends ParameterizedImageSource so that parameters can be passed 
 *  to it as a vector of values using the SetParameters() method. This method
 *  enables parameterized image sources to be used within ITK's optimization/
 *  registration framework.
 *
 * \ingroup DataSources
 */
template <class TOutputImage>
class ITK_EXPORT ParameterizedImageSource : public ImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef ParameterizedImageSource               Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Smart Pointer type to a DataObject. */
  typedef DataObject::Pointer DataObjectPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ParameterizedImageSource,ImageSource);

  /** Some convenient typedefs. */
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;

  typedef double                       ParametersValueType;
  typedef Array< ParametersValueType > ParametersType;

   /** ImageDimension constant */
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Set the parameters for this source. */
  virtual void SetParameters(const ParametersType& parameters) = 0;
  virtual ParametersType GetParameters() const = 0;
  virtual unsigned int GetNumberOfParameters() const = 0;
  
protected:
  ParameterizedImageSource();
  virtual ~ParameterizedImageSource() {}
    
private:
  ParameterizedImageSource(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_ParameterizedImageSource(_, EXPORT, x, y) namespace itk { \
  _(1(class EXPORT ParameterizedImageSource< ITK_TEMPLATE_1 x >)) \
  namespace Templates { typedef ParameterizedImageSource< ITK_TEMPLATE_1 x > ParameterizedImageSource##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkParameterizedImageSource+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkParameterizedImageSource.cxx"
#endif

#endif
