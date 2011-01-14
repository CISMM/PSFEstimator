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
#ifndef __itkGaussianPointSpreadFunctionImageSource_h
#define __itkGaussianPointSpreadFunctionImageSource_h

#include "itkParametricImageSource.h"
#include "itkFixedArray.h"
#include "itkGaussianImageSource.h"

namespace itk
{

/** \class GaussianPointSpreadFunctionImageSource
 * \brief Adapter class that addes the interface for the ParametricImageSource
 * to the GaussianImageSource.
 *
 * The adapted parameters are the Mean, standard deviation (Sigma), and Scale.
 *
 * This class was developed by Cory Quammen at the Center for Computer
 * Integrated Systems for Microscopy and Manipulation (http:://www.cismm.org)
 * at the University of North Carolina at Chapel Hill and was supported by
 * NIH grant P41-EB002025-25A1.
 *
 * \ingroup DataSources
 */
template <typename TOutputImage>
class ITK_EXPORT GaussianPointSpreadFunctionImageSource :
    public ParametricImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GaussianPointSpreadFunctionImageSource       Self;
  typedef ParametricImageSource<TOutputImage> Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  typedef typename Superclass::OutputImageType       OutputImageType;
  typedef typename Superclass::OutputImagePixelType  OutputImagePixelType;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename Superclass::SpacingType           SpacingType;
  typedef typename Superclass::PointType             PointType;
  typedef GaussianImageSource< OutputImageType >     GaussianImageSourceType;
  typedef typename GaussianImageSourceType::Pointer  GaussianImageSourcePointer;

  /** Dimensionality of the output image */
  itkStaticConstMacro(NDimensions, unsigned int, TOutputImage::ImageDimension);

  typedef typename GaussianImageSourceType::ArrayType ArrayType;
  typedef typename Superclass::SizeType               SizeType;
  typedef typename Superclass::SizeValueType          SizeValueType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GaussianPointSpreadFunctionImageSource,ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef double                       ParametersValueType;
  typedef Array< ParametersValueType > ParametersType;

  void SetSigma(const ArrayType& sigma);
  const ArrayType& GetSigma() const;

  void SetMean(const ArrayType& mean);
  const ArrayType& GetMean() const;

  /** Define methods defined by ParametricImageSource. */
  void                SetParameter(unsigned int index, ParametersValueType value);
  ParametersValueType GetParameter(unsigned int index) const;
  void           SetParameters(const ParametersType& parameters);
  ParametersType GetParameters() const;
  unsigned int   GetNumberOfParameters() const;

protected:
  GaussianPointSpreadFunctionImageSource();
  ~GaussianPointSpreadFunctionImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

  GaussianImageSourcePointer m_GaussianImageSource;

private:
  GaussianPointSpreadFunctionImageSource(const GaussianPointSpreadFunctionImageSource&); //purposely not implemented
  void operator=(const GaussianPointSpreadFunctionImageSource&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianPointSpreadFunctionImageSource.txx"
#endif

#endif
