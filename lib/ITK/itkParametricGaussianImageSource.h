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
#ifndef __itkParametricGaussianImageSource_h
#define __itkParametricGaussianImageSource_h

#include "itkParametricImageSource.h"
#include "itkFixedArray.h"
#include "itkGaussianImageSource.h"

namespace itk
{

/** \class ParametricGaussianImageSource
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
class ITK_EXPORT ParametricGaussianImageSource :
    public ParametricImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ParametricGaussianImageSource       Self;
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
  itkTypeMacro(ParametricGaussianImageSource,ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef double                       ParametersValueType;
  typedef Array< ParametersValueType > ParametersType;

  void SetSigma(const ArrayType& sigma)
  {
    if (sigma != this->m_GaussianImageSource->GetSigma())
      {
      this->m_GaussianImageSource->SetSigma(sigma);
      this->Modified();
      }
  }

  const ArrayType& GetSigma() const
  {
    return this->m_GaussianImageSource->GetSigma();
  }

  void SetMean(const ArrayType& mean)
  {
    if (mean != this->m_GaussianImageSource->GetMean())
      {
      this->m_GaussianImageSource->SetMean(mean);
      this->Modified();
      }
  }

  const ArrayType& GetMean() const
  {
    return this->m_GaussianImageSource->GetMean();
  }

  void SetScale(double scale)
  {
    if (scale != this->m_GaussianImageSource->GetScale())
      {
      this->m_GaussianImageSource->SetScale(scale);
      this->Modified();
      }
  }

  double GetScale() const
  {
    return this->m_GaussianImageSource->GetScale();
  }

  /** Define methods defined by ParametricImageSource. */
  void                SetParameter(unsigned int index, ParametersValueType value);
  ParametersValueType GetParameter(unsigned int index) const;
  void           SetParameters(const ParametersType& parameters);
  ParametersType GetParameters() const;
  unsigned int   GetNumberOfParameters() const;

protected:
  ParametricGaussianImageSource();
  ~ParametricGaussianImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

  GaussianImageSourcePointer m_GaussianImageSource;

private:
  ParametricGaussianImageSource(const ParametricGaussianImageSource&); //purposely not implemented
  void operator=(const ParametricGaussianImageSource&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParametricGaussianImageSource.txx"
#endif

#endif
