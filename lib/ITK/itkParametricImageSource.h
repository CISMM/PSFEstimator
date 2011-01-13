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
#ifndef __itkParametricImageSource_h
#define __itkParametricImageSource_h

#include "itkArray.h"
#include "itkAffineTransform.h"
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
  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename OutputImageType::RegionType     OutputImageRegionType;
  typedef typename OutputImageType::PointType      PointType;
  typedef typename OutputImageType::PointValueType PointValueType;
  typedef typename OutputImageType::SpacingType    SpacingType;
  typedef typename OutputImageType::IndexType      IndexType;
  typedef typename OutputImageType::RegionType     RegionType;
  typedef typename OutputImageType::SizeType       SizeType;
  typedef typename OutputImageType::SizeValueType  SizeValueType;

  typedef double                                   ParametersValueType;
  typedef Array< ParametersValueType >             ParametersType;

   /** ImageDimension constant */
  itkStaticConstMacro(OutputImageDimension,
                      unsigned int,
                      TOutputImage::ImageDimension);

  typedef AffineTransform< PointValueType, Self::OutputImageDimension >
    TransformType;
  typedef typename TransformType::Pointer TransformPointer;

  /** Set/get the image origin. */
  itkSetMacro(Origin, PointType);
  itkGetConstReferenceMacro(Origin, PointType);

  /** Set/get the image spacing. */
  itkSetMacro(Spacing, SpacingType);
  itkGetConstReferenceMacro(Spacing, SpacingType);

  /** Set/get the image size. */
  itkSetMacro(Size, SizeType);
  itkGetConstReferenceMacro(Size, SizeType);

  /** Set/get the transform. */
  itkSetObjectMacro(Transform, TransformType);
  itkGetObjectMacro(Transform, TransformType);

  /** Set the parameters for this source. */
  virtual void SetParameter(unsigned int index, ParametersValueType value) = 0;
  virtual ParametersValueType GetParameter(unsigned int index) const = 0;

  virtual void SetParameters(const ParametersType& parameters) = 0;
  virtual ParametersType GetParameters() const = 0;
  virtual unsigned int   GetNumberOfParameters() const = 0;

protected:
  ParametricImageSource();
  virtual ~ParametricImageSource() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  PointType        m_Origin;    // origin
  SpacingType      m_Spacing;   // spacing
  SizeType         m_Size;      // size of the output image
  TransformPointer m_Transform; // transform to apply to image sample coordinates

  virtual void GenerateOutputInformation();

private:
  ParametricImageSource(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParametricImageSource.txx"
#endif

#endif
