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
#ifndef __itkMicroscopePointSpreadFunctionImageSource_h
#define __itkMicroscopePointSpreadFunctionImageSource_h

#include "itkParametricImageSource.h"
#include "itkMacro.h"

namespace itk
{

/** \class MicroscopePointSpreadFunctionImageSource
 * \brief Base class for point-spread function image sources specific
 * to microscopy. This class defines terms common to all optical
 * microscopes: magnification, numerical aperture of the objective
 * lens, and the spatial position of the point source of light.
 *
 * \ingroup DataSources Multithreaded
 */
template< class TOutputImage >
class ITK_EXPORT MicroscopePointSpreadFunctionImageSource :
  public ParametricImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef MicroscopePointSpreadFunctionImageSource Self;
  typedef ParametricImageSource< TOutputImage >    Superclass;
  typedef SmartPointer< Self >                     Pointer;
  typedef SmartPointer< const Self >               ConstPointer;

  /** Typedef for the output image. */
  typedef TOutputImage                             OutputImageType;
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
  itkTypeMacro(MicroscopePointSpreadFunctionImageSource,
               ParametricImageSource);

  /** Set/get  the magnification. */
  itkSetMacro(Magnification, double);
  itkGetConstMacro(Magnification, double);

  /** Set/get the numerical aperture (NA). */
  itkSetMacro(NumericalAperture, double);
  itkGetConstMacro(NumericalAperture, double);

protected:
  MicroscopePointSpreadFunctionImageSource();
  ~MicroscopePointSpreadFunctionImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  double m_Magnification;
  double m_NumericalAperture;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMicroscopePointSpreadFunctionImageSource.txx"
#endif

#endif // __itkMicroscopePointSpreadFunctionImageSource_h
