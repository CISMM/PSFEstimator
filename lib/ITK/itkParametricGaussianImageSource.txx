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
#ifndef __itkParametricGaussianImageSource_txx
#define __itkParametricGaussianImageSource_txx

#include "itkParametricGaussianImageSource.h"

namespace itk
{

template <class TOutputImage>
ParametricGaussianImageSource<TOutputImage>
::ParametricGaussianImageSource()
{
  this->m_GaussianImageSource = GaussianImageSourceType::New();
}

template <class TOutputImage>
ParametricGaussianImageSource<TOutputImage>
::~ParametricGaussianImageSource()
{
}


template<typename TOutputImage>
void
ParametricGaussianImageSource<TOutputImage>
::SetParameter(unsigned int index, ParametersValueType value)
{
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  if (index >= 0 && index < dimensions)
    {
    ArrayType sigma = this->GetSigma();
    sigma[index] = value;
    this->SetSigma(sigma);
    }
  else if (index >= dimensions && index < 2*dimensions)
    {
    ArrayType mean = this->GetMean();
    mean[index-dimensions] = value;
    this->SetMean(mean);
    }
  else
    {
    this->SetScale(value);
    }
}


template<typename TOutputImage>
typename ParametricGaussianImageSource<TOutputImage>::ParametersValueType
ParametricGaussianImageSource<TOutputImage>
::GetParameter(unsigned int index) const
{
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  if (index >= 0 && index < dimensions)
    {
    return this->GetSigma()[index];
    }
  else if (index >= dimensions && index < 2*dimensions)
    {
    return this->GetMean()[index-dimensions];
    }
  else
    {
    return this->GetScale();
    }
}


template<typename TOutputImage>
void
ParametricGaussianImageSource<TOutputImage>
::SetParameters(const ParametersType& parameters)
{
  unsigned int i;
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  ArrayType sigma, mean;
  for (i=0; i<dimensions; i++)
    {
    sigma[i] = parameters[i];
    mean[i]  = parameters[dimensions + i];
    }
  this->SetSigma(sigma);
  this->SetMean(mean);
  this->SetScale(parameters[2*dimensions]);
}


template<typename TOutputImage>
typename ParametricGaussianImageSource<TOutputImage>::ParametersType
ParametricGaussianImageSource<TOutputImage>
::GetParameters() const
{
  ParametersType parameters(GetNumberOfParameters());
  ArrayType sigma = this->GetSigma();
  ArrayType mean  = this->GetMean();

  unsigned int i;
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  for (i=0; i<dimensions; i++)
    {
    parameters[i] = sigma[i];
    parameters[dimensions + i] = mean[i];
    }
  parameters[2*dimensions] = this->GetScale();

  return parameters;
}


template<typename TOutputImage>
unsigned int
ParametricGaussianImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return 2*itkGetStaticConstMacro(OutputImageDimension) + 1;
}


template<typename TOutputImage>
void
ParametricGaussianImageSource<TOutputImage>
::GenerateData()
{
  m_GaussianImageSource->GraftOutput(this->GetOutput());
  this->GraftOutput(m_GaussianImageSource->GetOutput());
}


template <class TOutputImage>
void
ParametricGaussianImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  std::cout << this->m_GaussianImageSource << std::endl;
}


} // end namespace itk

#endif
