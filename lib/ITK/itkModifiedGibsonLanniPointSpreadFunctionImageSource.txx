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
#ifndef __itkModifiedGibsonLanniPointSpreadFunctionImageSource_txx
#define __itkModifiedGibsonLanniPointSpreadFunctionImageSource_txx

#include "itkModifiedGibsonLanniPointSpreadFunctionImageSource.h"

#include "itkGibsonLanniPointSpreadFunctionImageSource.txx"
#include "itkGaussianImageSource.txx"

namespace itk
{

template< class TOutputImage >
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::ModifiedGibsonLanniPointSpreadFunctionImageSource()
{
  m_GibsonLanniSource = GibsonLanniSourceType::New();
  m_GaussianSource = GaussianSourceType::New();
  m_GaussianSource->SetNormalized(false);

  m_SubtractFilter = SubtractImageFilterType::New();
  m_SubtractFilter->SetInput1( m_GibsonLanniSource->GetOutput() );
  m_SubtractFilter->SetInput2( m_GaussianSource->GetOutput() );

  typename GaussianSourceType::ArrayType mean;
  mean.Fill(0.0);
  m_GaussianSource->SetMean(mean);

  typename GaussianSourceType::ArrayType sigma;
  sigma.Fill(10.0);
  m_GaussianSource->SetSigma(sigma);
  m_GaussianSource->SetScale(0.0);

  m_ModifiedEventCommand = MemberCommandType::New();
  m_ModifiedEventCommand->SetCallbackFunction(this, &Self::ModifiedCallback);

  m_GibsonLanniSourceObserverTag = m_GibsonLanniSource->
    AddObserver(ModifiedEvent(), m_ModifiedEventCommand);
  m_GaussianSourceObserverTag = m_GaussianSource->
    AddObserver(ModifiedEvent(), m_ModifiedEventCommand);
}


template< class TOutputImage >
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::~ModifiedGibsonLanniPointSpreadFunctionImageSource()
{
  m_GibsonLanniSource->RemoveObserver(m_GibsonLanniSourceObserverTag);
  m_GaussianSource->RemoveObserver(m_GaussianSourceObserverTag);
}


template< class TOutputImage >
void
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::ModifiedCallback()
{
  this->Modified();
}


template< class TOutputImage >
void
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::SetParameter(unsigned int index, ParametersValueType value)
{
  unsigned int numGibsonLanniParameters = m_GibsonLanniSource->GetNumberOfParameters();
  if ( index < numGibsonLanniParameters )
    {
    m_GibsonLanniSource->SetParameter(index, value);
    }
  else
    {
    // Set Gaussian parameter
    typename GaussianSourceType::ArrayType mean  = m_GaussianSource->GetMean();
    typename GaussianSourceType::ArrayType sigma = m_GaussianSource->GetSigma();
    index -= numGibsonLanniParameters;
    switch ( index )
      {
      case 0:
      case 1:
      case 2:
        mean[index] = value;
        m_GaussianSource->SetMean(mean);
        break;

      case 3:
      case 4:
      case 5:
        sigma[index - 3] = value;
        m_GaussianSource->SetSigma(sigma);
        break;

      case 6:
        m_GaussianSource->SetScale(value);
        break;

      default:
        break;
      }
    }
}


template< class TOutputImage >
typename ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >::ParametersValueType
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::GetParameter(unsigned int index) const
{
  unsigned int numGibsonLanniParameters = m_GibsonLanniSource->GetNumberOfParameters();
  if ( index < numGibsonLanniParameters )
    {
    return m_GibsonLanniSource->GetParameter(index);
    }
  else
    {
    // Get Gaussian parameter
    typename GaussianSourceType::ArrayType mean  = m_GaussianSource->GetMean();
    typename GaussianSourceType::ArrayType sigma = m_GaussianSource->GetSigma();
    index -= numGibsonLanniParameters;
    switch ( index )
      {
      case 0:
      case 1:
      case 2:
        return mean[index];
        break;

      case 3:
      case 4:
      case 5:
        return sigma[index - 3];
        break;

      case 6:
        return m_GaussianSource->GetScale();
        break;

      default:
        break;
      }
    }

  return 0.0;
}


template< class TOutputImage >
void
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::SetParameters(const ParametersType& parameters)
{
  m_GibsonLanniSource->SetParameters(parameters);

  // Set the Gaussian parameters
  typename GaussianSourceType::ArrayType mean  = m_GaussianSource->GetMean();
  typename GaussianSourceType::ArrayType sigma = m_GaussianSource->GetSigma();

  int index = m_GibsonLanniSource->GetNumberOfParameters();
  mean[0] = parameters[index++];
  mean[1] = parameters[index++];
  mean[2] = parameters[index++];
  m_GaussianSource->SetMean(mean);

  sigma[0] = parameters[index++];
  sigma[1] = parameters[index++];
  sigma[2] = parameters[index++];
  m_GaussianSource->SetSigma(sigma);

  m_GaussianSource->SetScale(parameters[index++]);
}


template< class TOutputImage >
typename ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >::ParametersType
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::GetParameters() const
{
  // Concatenate the Gibson-Lanni parameters and the Gaussian
  // parameters.
  ParametersType gibsonLanniParams = m_GibsonLanniSource->GetParameters();
  ParametersType parameters( this->GetNumberOfParameters() );
  unsigned int index;
  for ( index = 0; index < m_GibsonLanniSource->GetNumberOfParameters(); index++ )
    {
    parameters[index] = gibsonLanniParams[index];
    }

  // Get the Gaussian parameters
  typename GaussianSourceType::ArrayType mean  = m_GaussianSource->GetMean();
  typename GaussianSourceType::ArrayType sigma = m_GaussianSource->GetSigma();

  parameters[index++] = mean[0];
  parameters[index++] = mean[1];
  parameters[index++] = mean[2];

  parameters[index++] = sigma[0];
  parameters[index++] = sigma[1];
  parameters[index++] = sigma[2];

  parameters[index++] = m_GaussianSource->GetScale();

  return parameters;
}


template< class TOutputImage >
unsigned int
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::GetNumberOfParameters() const
{
  return m_GibsonLanniSource->GetNumberOfParameters() + 7;
}


template< class TOutputImage >
void
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::GenerateData()
{
  m_GibsonLanniSource->SetSize( this->GetSize() );
  m_GibsonLanniSource->SetOrigin( this->GetOrigin() );
  m_GibsonLanniSource->SetSpacing( this->GetSpacing() );

  m_GaussianSource->SetSize( this->GetSize() );
  m_GaussianSource->SetOrigin( this->GetOrigin() );
  m_GaussianSource->SetSpacing( this->GetSpacing() );

  m_SubtractFilter->GraftOutput( this->GetOutput() );
  m_SubtractFilter->Update();
  this->GraftOutput( m_SubtractFilter->GetOutput() );
}


template< class TOutputImage >
void
ModifiedGibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  std::cout << indent << m_GibsonLanniSource << std::endl;
  std::cout << indent << m_GaussianSource << std::endl;
}


} // end namespace itk

#endif // __itkModifiedGibsonLanniPointSpreadFunctionImageSource_txx
