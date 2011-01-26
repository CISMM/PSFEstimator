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
  m_KernelIsRadiallySymmetric = false;

  m_Convolver = ConvolverType::New();

  // Specify multiple integration samples in x and y but not z.
  typename ConvolverType::SizeType voxelSamples = {{1, 1, 1}};
  m_Convolver->SetNumberOfIntegrationSamples(voxelSamples);
  m_Convolver->WeightIntegrationByAreaOn();

  m_RescaleFilter = RescaleImageFilterType::New();
  m_RescaleFilter->SetInput(m_Convolver->GetOutput());

  m_ModifiedEventCommand = MemberCommandType::New();
  m_ModifiedEventCommand->SetCallbackFunction(this, &Self::KernelModified);
  m_ObserverTag = 0;
}


template< class TOutputImage >
BeadSpreadFunctionImageSource< TOutputImage >
::~BeadSpreadFunctionImageSource()
{
  this->m_KernelSource->RemoveObserver(this->m_ObserverTag);
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::KernelModified()
{
  this->Modified();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetSize(const SizeType& size)
{
  if (size != m_Convolver->GetSize())
    {
    this->Modified();
    }
  m_Convolver->SetSize(size);
}


template< class TOutputImage >
const typename BeadSpreadFunctionImageSource< TOutputImage >::SizeType&
BeadSpreadFunctionImageSource< TOutputImage >
::GetSize() const
{
  return m_Convolver->GetSize();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetSpacing(const SpacingType& spacing)
{
  if (spacing != m_Convolver->GetSpacing())
    {
    this->Modified();
    }
  m_Convolver->SetSpacing(spacing);
}


template< class TOutputImage >
const typename BeadSpreadFunctionImageSource< TOutputImage >::SpacingType&
BeadSpreadFunctionImageSource< TOutputImage >
::GetSpacing() const
{
  return m_Convolver->GetSpacing();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetOrigin(const PointType& origin)
{
  if (origin != m_Convolver->GetOrigin())
    {
    this->Modified();
    }
  m_Convolver->SetOrigin(origin);
}


template< class TOutputImage >
const typename BeadSpreadFunctionImageSource< TOutputImage >::PointType&
BeadSpreadFunctionImageSource< TOutputImage >
::GetOrigin() const
{
  return m_Convolver->GetOrigin();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetBeadCenter(const PointType& center)
{
  if (center != m_Convolver->GetSphereCenter())
    {
    m_Convolver->SetSphereCenter(center);
    this->Modified();
    }
}


template< class TOutputImage >
const typename BeadSpreadFunctionImageSource< TOutputImage >::PointType&
BeadSpreadFunctionImageSource< TOutputImage >
::GetBeadCenter() const
{
  return m_Convolver->GetSphereCenter();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetBeadRadius(double radius)
{
  if (radius != m_Convolver->GetSphereRadius())
    {
    m_Convolver->SetSphereRadius(radius);
    this->Modified();
    }
}


template< class TOutputImage >
double
BeadSpreadFunctionImageSource< TOutputImage >
::GetBeadRadius() const
{
  return m_Convolver->GetSphereRadius();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetShearX(double shear)
{
  if (shear != m_Convolver->GetShearX())
    {
    m_Convolver->SetShearX(shear);
    this->Modified();
    }
}


template< class TOutputImage >
double
BeadSpreadFunctionImageSource< TOutputImage >
::GetShearX() const
{
  return m_Convolver->GetShearX();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetShearY(double shear)
{
  if (shear != m_Convolver->GetShearY())
    {
    m_Convolver->SetShearY(shear);
    this->Modified();
    }
}


template< class TOutputImage >
double
BeadSpreadFunctionImageSource< TOutputImage >
::GetShearY() const
{
  return m_Convolver->GetShearY();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource< TOutputImage >
::SetKernelSource( KernelImageSourceType* source )
{
  if ( this->m_KernelSource != source )
    {
    if ( this->m_KernelSource )
      {
      this->m_KernelSource->RemoveObserver(this->m_ObserverTag);
      }
    this->m_KernelSource = source;
    this->m_ObserverTag = this->m_KernelSource->
      AddObserver(ModifiedEvent() , m_ModifiedEventCommand);
    this->Modified();
    }

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
        this->SetIntensityShift(value);
        break;

      case 10:
        this->SetIntensityScale(value);
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
  unsigned int numberOfBSFParameters = this->GetNumberOfBeadSpreadFunctionParameters();
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
        return this->GetIntensityShift();
        break;

      case 10:
        return this->GetIntensityScale();
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
  this->SetSpacing(spacing);

  this->SetBeadRadius(parameters[index++]);

  PointType center;
  for (int i = 0; i < ImageDimension; i++)
    {
    center[i] = parameters[index++];
    }
  this->SetBeadCenter(center);

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
  parameters[index++] = this->GetIntensityShift();
  parameters[index++] = this->GetIntensityScale();

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
  return this->m_KernelSource->GetNumberOfParameters() +
    this->GetNumberOfBeadSpreadFunctionParameters();
}


template< class TOutputImage >
unsigned int
BeadSpreadFunctionImageSource< TOutputImage >
::GetNumberOfBeadSpreadFunctionParameters() const
{
  return 2*ImageDimension + 5;
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
  SizeType    psfTableSize;
  SpacingType psfTableSpacing(50.0); // An arbitrary spacing

  // Determine necessary spatial extent of PSF table.
  PointType minExtent(this->GetOrigin());
  PointType maxExtent;
  const unsigned int dimensions = itkGetStaticConstMacro(OutputImageDimension);
  for ( unsigned int i = 0; i < dimensions; i++ )
    {
    // First calculate extent of BSF in this dimension.
    maxExtent[i] = static_cast<PointValueType>
      (this->GetSize()[i]-1) * this->GetSpacing()[i] + minExtent[i];

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

  // Generate just a radial profile of the PSF if it is radially symmetric
  if (this->m_KernelIsRadiallySymmetric)
    {

    // Calculate distance from image corners to bead center, projected
    // to the xy-plane.
    typedef Point< double, 2> Point2DType;
    Point2DType beadCenter(GetBeadCenter().GetDataPointer());

    Point2DType pt[4];
    pt[0][0] = minExtent[0];  pt[0][1] = minExtent[1];
    pt[1][0] = minExtent[0];  pt[1][1] = maxExtent[1];
    pt[2][0] = maxExtent[0];  pt[2][1] = minExtent[1];
    pt[3][0] = maxExtent[0];  pt[3][1] = maxExtent[1];

    double maxRadialDistance = NumericTraits<double>::min();
    for ( unsigned int i = 0; i < 4; i++)
      {
      double distance = pt[i].EuclideanDistanceTo(beadCenter);
      if (distance > maxRadialDistance) maxRadialDistance = distance;
      }

    // Need to change some values here
    psfTableOrigin[0] = 0.0;
    psfTableOrigin[1] = 0.0;
    long maxRadialSize = Math::Ceil<long>(maxRadialDistance / psfTableSpacing[0]);
    psfTableSize[0] = maxRadialSize;
    psfTableSize[1] = 1;
    }

  m_KernelSource->SetSize(psfTableSize);
  m_KernelSource->SetSpacing(psfTableSpacing);
  m_KernelSource->SetOrigin(psfTableOrigin);
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
  m_KernelSource->Print(os,indent);
  m_Convolver->Print(os,indent);
  m_RescaleFilter->Print(os,indent);
}


} // end namespace itk

#endif // _itkBeadSpreadFunctionImageSource_txx
