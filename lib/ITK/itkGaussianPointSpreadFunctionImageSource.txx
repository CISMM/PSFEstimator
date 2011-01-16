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
#ifndef __itkGaussianPointSpreadFunctionImageSource_txx
#define __itkGaussianPointSpreadFunctionImageSource_txx

#include "itkGaussianPointSpreadFunctionImageSource.h"

namespace itk
{

template <class TOutputImage>
GaussianPointSpreadFunctionImageSource<TOutputImage>
::GaussianPointSpreadFunctionImageSource()
{
  this->m_GaussianImageSource = GaussianImageSourceType::New();

  ArrayType sigma(200.0);
  this->SetSigma(sigma);
  ArrayType mean(0.0);
  this->SetMean(mean);
}

template <class TOutputImage>
GaussianPointSpreadFunctionImageSource<TOutputImage>
::~GaussianPointSpreadFunctionImageSource()
{
}


template<typename TOutputImage>
void
GaussianPointSpreadFunctionImageSource<TOutputImage>
::SetSigma(const ArrayType& sigma)
{
  if (sigma != this->m_GaussianImageSource->GetSigma())
    {
    this->m_GaussianImageSource->SetSigma(sigma);
    this->Modified();
    }
}


template< typename TOutputImage >
const typename GaussianPointSpreadFunctionImageSource< TOutputImage >::ArrayType&
GaussianPointSpreadFunctionImageSource< TOutputImage >
::GetSigma() const
{
  return this->m_GaussianImageSource->GetSigma();
}


template< typename TOutputImage >
void
GaussianPointSpreadFunctionImageSource< TOutputImage >
::SetMean(const ArrayType& mean)
{
  if (mean != this->m_GaussianImageSource->GetMean())
    {
    this->m_GaussianImageSource->SetMean(mean);
    this->Modified();
    }
}


template< typename TOutputImage >
const typename GaussianPointSpreadFunctionImageSource< TOutputImage >::ArrayType&
GaussianPointSpreadFunctionImageSource< TOutputImage >
::GetMean() const
{
  return this->m_GaussianImageSource->GetMean();
}


template<typename TOutputImage>
void
GaussianPointSpreadFunctionImageSource<TOutputImage>
::SetParameter(unsigned int index, ParametersValueType value)
{
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  if (index >= 0 && index < dimensions)
    {
    ArrayType sigma = this->GetSigma();
    sigma[index] = value;
    this->SetSigma(sigma);
    }
}


template<typename TOutputImage>
typename GaussianPointSpreadFunctionImageSource<TOutputImage>::ParametersValueType
GaussianPointSpreadFunctionImageSource<TOutputImage>
::GetParameter(unsigned int index) const
{
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  if (index >= 0 && index < dimensions)
    {
    return this->GetSigma()[index];
    }

  return 0.0;
}


template<typename TOutputImage>
void
GaussianPointSpreadFunctionImageSource<TOutputImage>
::SetParameters(const ParametersType& parameters)
{
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  ArrayType sigma;
  for ( unsigned int i=0; i<dimensions; i++)
    {
    sigma[i] = parameters[i];
    }
  this->SetSigma(sigma);
}


template<typename TOutputImage>
typename GaussianPointSpreadFunctionImageSource<TOutputImage>::ParametersType
GaussianPointSpreadFunctionImageSource<TOutputImage>
::GetParameters() const
{
  ParametersType parameters(GetNumberOfParameters());
  ArrayType sigma = this->GetSigma();

  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  for (unsigned int i=0; i<dimensions; i++)
    {
    parameters[i] = sigma[i];
    }

  return parameters;
}


template<typename TOutputImage>
unsigned int
GaussianPointSpreadFunctionImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return itkGetStaticConstMacro(OutputImageDimension);
}


template<typename TOutputImage>
void
GaussianPointSpreadFunctionImageSource<TOutputImage>
::GenerateData()
{
  m_GaussianImageSource->SetSize(this->m_Size);
  m_GaussianImageSource->SetSpacing(this->m_Spacing);
  m_GaussianImageSource->SetOrigin(this->m_Origin);

  m_GaussianImageSource->GraftOutput(this->GetOutput());

  m_GaussianImageSource->UpdateLargestPossibleRegion();

  this->GraftOutput(m_GaussianImageSource->GetOutput());
}


template <class TOutputImage>
void
GaussianPointSpreadFunctionImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  std::cout << this->m_GaussianImageSource << std::endl;
}


} // end namespace itk

#endif
