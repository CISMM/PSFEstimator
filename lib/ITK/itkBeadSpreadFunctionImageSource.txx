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
#ifndef _itkBeadSpreadFunctionImageSource_txx
#define _itkBeadSpreadFunctionImageSource_txx

#include "itkBeadSpreadFunctionImageSource.h"


namespace itk
{

template< class TOutputImage >
BeadSpreadFunctionImageSource< TOutputImage >
::BeadSpreadFunctionImageSource()
{
  m_IntensityShift = 0.0;
  m_IntensityScale = 1.0;

  m_KernelSource = NULL;

  m_Convolver = ConvolverType::New();

  m_RescaleFilter = RescaleImageFilterType::New();
  m_RescaleFilter->SetInput(m_Convolver->GetOutput());
}


template< class TOutputImage >
BeadSpreadFunctionImageSource< TOutputImage >
::~BeadSpreadFunctionImageSource()
{
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetParameter(unsigned int index, ParametersValueType value)
{
  unsigned int numberOfBSFParameters = 2 * ImageDimension + 5;
  if (index < numberOfBSFParameters)
    {
    SpacingType spacing = this->GetSpacing();
    PointType   center  = this->GetBeadCenter();

    switch (index)
      {
      case 0:
      case 1:
      case 2:
        spacing[index] = value;
        this->SetSpacing(spacing);
        break;

      case 3:
        this->SetBeadRadius(value);
        break;

      case 4:
      case 5:
      case 6:
        center[index - 4] = value;
        this->SetBeadCenter(center);
        break;

      case 7:
        this->SetShearX(value);
        break;

      case 8:
        this->SetShearY(value);
        break;

      case 9:
        this->SetIntensityScale(value);
        break;

      case 10:
        this->SetIntensityShift(value);
        break;
      }
    }
  else
    {
    this->m_KernelSource->SetParameter(index - numberOfBSFParameters, value);
    }
}


template< class TOutputImage >
typename BeadSpreadFunctionImageSource< TOutputImage >::ParametersValueType
BeadSpreadFunctionImageSource< TOutputImage >
::GetParameter(unsigned int index) const
{
  unsigned int numberOfBSFParameters = 2 * ImageDimension + 5;
  if (index < numberOfBSFParameters)
    {
    SpacingType spacing = this->GetSpacing();
    PointType   center  = this->GetBeadCenter();

    switch (index)
      {
      case 0:
      case 1:
      case 2:
        return spacing[index];
        break;

      case 3:
        return this->GetBeadRadius();
        break;

      case 4:
      case 5:
      case 6:
        return center[index - 4];
        break;

      case 7:
        return this->GetShearX();
        break;

      case 8:
        return this->GetShearY();
        break;

      case 9:
        return this->GetIntensityScale();
        break;

      case 10:
        return this->GetIntensityShift();
        break;

      default:
        return 999.0;
      }
    }
  else
    {
    return this->m_KernelSource->GetParameter(index - numberOfBSFParameters);
    }
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetParameters(const ParametersType& parameters)
{
  int index = 0;

  // The first parameters are bead-spread function parameters
  SpacingType spacing;
  for (int i = 0; i < ImageDimension; i++)
    {
    spacing[i] = parameters[index++];
    }

  this->SetBeadRadius(parameters[index++]);

  PointType center;
  for (int i = 0; i < ImageDimension; i++)
    {
    center[i] = parameters[index++];
    }

  this->SetShearX(parameters[index++]);
  this->SetShearY(parameters[index++]);

  this->SetIntensityShift(parameters[index++]);
  this->SetIntensityScale(parameters[index++]);

  // The last parameters go to the kernel source
  ParametersType kernelParameters(this->m_KernelSource->GetNumberOfParameters());
  for (unsigned int i = 0; i < kernelParameters.GetSize(); i++)
    {
    kernelParameters[i] = parameters[index++];
    }

  this->m_KernelSource->SetParameters(kernelParameters);
}


template< class TOutputImage >
typename BeadSpreadFunctionImageSource< TOutputImage >::ParametersType
BeadSpreadFunctionImageSource< TOutputImage >
::GetParameters() const
{
  ParametersType parameters(GetNumberOfParameters());
  int index = 0;

  // The first parameters come from the bead-spread function
  const SpacingType spacing = this->GetSpacing();
  for (int i = 0; i < ImageDimension; i++)
    {
    parameters[index++] = spacing[i];
    }

  parameters[index++] = this->GetBeadRadius();

  const PointType beadCenter = this->GetBeadCenter();
  for (int i = 0; i < ImageDimension; i++)
    {
    parameters[index++] = beadCenter[i];
    }

  parameters[index++] = this->GetShearX();
  parameters[index++] = this->GetShearY();

  // The last parameters come from the kernel source
  ParametersType kernelParameters = this->m_KernelSource->GetParameters();
  for (unsigned int i = 0; i < kernelParameters.GetSize(); i++)
    {
    parameters[index++] = kernelParameters[i];
    }

  return parameters;
}


template< class TOutputImage >
unsigned int
BeadSpreadFunctionImageSource< TOutputImage >
::GetNumberOfParameters() const
{
  return this->m_KernelSource->GetNumberOfParameters() + 2*ImageDimension + 5;
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetZCoordinate(unsigned int index, double coordinate)
{
  m_Convolver->SetZCoordinate(index, coordinate);
}


template< class TOutputImage >
double
BeadSpreadFunctionImageSource< TOutputImage >
::GetZCoordinate(unsigned int index)
{
  return m_Convolver->GetZCoordinate(index);
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetUseCustomZCoordinates(bool use)
{
  m_Convolver->SetUseCustomZCoordinates(use);
  this->Modified();
}


template< class TOutputImage >
bool
BeadSpreadFunctionImageSource< TOutputImage >
::GetUseCustomZCoordinates()
{
  return m_Convolver->GetUseCustomZCoordinates();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::GenerateData()
{
  // Set the PSF sampling spacing and size parameters, and update.
  PointType   psfTableOrigin;
  SpacingType psfTableSpacing;
  SizeType    psfTableSize;

  psfTableSpacing[0] = 40.0; // A somewhat arbitrary spacing.
  psfTableSpacing[1] = 40.0;
  psfTableSpacing[2] = 100.0;

  // Determine necessary spatial extent of PSF table.
  PointType minExtent;
  PointType maxExtent;
  for ( int i = 0; i < 3; i++ )
    {
    // First calculate extent of BSF in this dimension.
    minExtent[i] = GetOrigin()[i];
    maxExtent[i] = static_cast<PointValueType>(GetSize()[i]-1)*GetSpacing()[i] +
      minExtent[i];

    // Now modify calculated PSF dimension to account for bead shift and radius
    minExtent[i] += -GetBeadCenter()[i] - GetBeadRadius();
    maxExtent[i] += -GetBeadCenter()[i] + GetBeadRadius();

    // Determine logical extent of the PSF table for the min and max extents.
    long iDimMin = Math::Floor<long>(minExtent[i] / psfTableSpacing[i]);
    psfTableOrigin[i] = static_cast<double>(iDimMin) * psfTableSpacing[i];
    long iDimMax = Math::Ceil<long>(maxExtent[i] / psfTableSpacing[i]);

    // Determine the logical extent of the PSF table in this dimension.
    psfTableSize[i] = iDimMax - iDimMin + 1;
    }

  // We need only half the radial profile image here.
  PointType   profileOrigin(psfTableOrigin);
  SpacingType profileSpacing(psfTableSpacing);
  SizeType    profileSize(psfTableSize);

  // Calculate distance from image corners to bead center, projected
  // to the xy-plane.
  minExtent[2] = 0.0;
  maxExtent[2] = 0.0;
  PointType beadCenter(GetBeadCenter());
  beadCenter[2] = 0.0;

  PointType pt[4];
  pt[0][0] = minExtent[0];  pt[0][1] = minExtent[1];  pt[0][2] = 0.0;
  pt[1][0] = minExtent[0];  pt[1][1] = maxExtent[1];  pt[1][2] = 0.0;
  pt[2][0] = maxExtent[0];  pt[2][1] = minExtent[1];  pt[2][2] = 0.0;
  pt[3][0] = maxExtent[0];  pt[3][1] = maxExtent[1];  pt[3][2] = 0.0;

  double maxRadialDistance = NumericTraits<double>::min();
  for ( unsigned int i = 0; i < 4; i++)
    {
    VectorType v = pt[i] - beadCenter;
    double distance = v.GetNorm();
    if (distance > maxRadialDistance)
      {
      maxRadialDistance = distance;
      }
    }

  // Need to change some values here
  profileSpacing[0] = 0.5 * psfTableSpacing[0];
  profileSpacing[1] = 0.5 * psfTableSpacing[1];
  profileOrigin[0] = 0.0;
  profileOrigin[1] = 0.0;
  long maxRadialSize = Math::Ceil<long>(maxRadialDistance / profileSpacing[0]);
  profileSize[0] = maxRadialSize;
  profileSize[1] = 1;

  m_KernelSource->SetSize(profileSize);
  m_KernelSource->SetSpacing(profileSpacing);
  m_KernelSource->SetOrigin(profileOrigin);
  m_KernelSource->UpdateLargestPossibleRegion();

  m_Convolver->SetInput(m_KernelSource->GetOutput());
  m_Convolver->UpdateLargestPossibleRegion();

  m_RescaleFilter->GraftOutput(this->GetOutput());
  m_RescaleFilter->SetShift(m_IntensityShift);
  m_RescaleFilter->SetScale(m_IntensityScale);
  m_RescaleFilter->UpdateLargestPossibleRegion();
  this->GraftOutput(m_RescaleFilter->GetOutput());
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::GenerateOutputInformation()
{
  OutputImageType *output;
  IndexType index = {{0}};
  SizeType size( m_Convolver->GetSize() );

  output = this->GetOutput(0);

  RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( size );
  largestPossibleRegion.SetIndex( index );
  output->SetLargestPossibleRegion( largestPossibleRegion );

  output->SetSpacing( m_Convolver->GetSpacing() );
  output->SetOrigin( m_Convolver->GetOrigin() );
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << m_KernelSource << std::endl;
  os << indent << m_Convolver << std::endl;
  os << indent << m_RescaleFilter << std::endl;
}


} // end namespace itk

#endif // _itkBeadSpreadFunctionImageSource_txx
