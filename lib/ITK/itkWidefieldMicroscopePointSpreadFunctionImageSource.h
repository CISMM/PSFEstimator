/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkWidefieldMicroscopePointSpreadFunctionImageSource_h
#define __itkWidefieldMicroscopePointSpreadFunctionImageSource_h

#include "itkMicroscopePointSpreadFunctionImageSource.h"

namespace itk
{

/** \class WidefieldMicroscopePointSpreadFunctionImageSource
 * \brief Base class for point-spread function image sources from
 * widefield microscopes. This class has a parameter for the emission
 * wavelength of light.
 *
 * \ingroup DataSources Multithreaded
 */
template< class TOutputImage >
class ITK_EXPORT WidefieldMicroscopePointSpreadFunctionImageSource :
  public MicroscopePointSpreadFunctionImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef WidefieldMicroscopePointSpreadFunctionImageSource        Self;
  typedef MicroscopePointSpreadFunctionImageSource< TOutputImage > Superclass;
  typedef SmartPointer< Self >                                     Pointer;
  typedef SmartPointer< const Self >                               ConstPointer;

  /** Typedef for the output image. */
  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::PixelType      PixelType;
  typedef typename OutputImageType::IndexType      IndexType;
  typedef typename OutputImageType::RegionType     RegionType;
  typedef typename OutputImageType::PointType      PointType;
  typedef typename OutputImageType::PointValueType PointValueType;
  typedef typename OutputImageType::SpacingType    SpacingType;
  typedef typename OutputImageType::SizeType       SizeType;
  typedef typename OutputImageType::SizeValueType  SizeValueType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(WidefieldMicroscopePointSpreadFunctionImageSource,
               MicroscopePointSpreadFunctionImageSource);

  /** Set/get the emission wavelength. */
  itkSetMacro(EmissionWavelength, double);
  itkGetConstMacro(EmissionWavelength, double);

protected:
  WidefieldMicroscopePointSpreadFunctionImageSource();
  ~WidefieldMicroscopePointSpreadFunctionImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  double m_EmissionWavelength;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWidefieldMicroscopePointSpreadFunctionImageSource.txx"
#endif

#endif // __itkWidefieldMicroscopePointSpreadFunctionImageSource_h
