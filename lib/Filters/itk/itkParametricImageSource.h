/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkParametricImageSource.h,v $
  Language:  C++
  Date:      $Date: 2010/04/19 18:50:02 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkParametricImageSource_h
#define __itkParametricImageSource_h

#include "itkArray.h"
#include "itkImageSource.h"

namespace itk
{

/** \class ParametricImageSource
 *  \brief Extends ImageSource so that parameters can be passed
 *  to it as a vector of values using the SetParameters() method. This method
 *  enables parameterized image sources to be used within ITK's optimization/
 *  registration framework.
 *
 * \ingroup DataSources
 */
template <class TOutputImage>
class ITK_EXPORT ParametricImageSource : public ImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef ParametricImageSource     Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Smart Pointer type to a DataObject. */
  typedef DataObject::Pointer DataObjectPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ParametricImageSource,ImageSource);

  /** Some convenient typedefs. */
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;

  typedef double                               ParametersValueType;
  typedef Array< ParametersValueType >         ParametersType;

   /** ImageDimension constant */
  itkStaticConstMacro(OutputImageDimension,
                      unsigned int,
                      TOutputImage::ImageDimension);

  /** Set the parameters for this source. */
  virtual void           SetParameters(const ParametersType& parameters) = 0;
  virtual ParametersType GetParameters() const = 0;
  virtual unsigned int   GetNumberOfParameters() const = 0;

protected:
  ParametricImageSource();
  virtual ~ParametricImageSource() {}

private:
  ParametricImageSource(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParametricImageSource.txx"
#endif

#endif
