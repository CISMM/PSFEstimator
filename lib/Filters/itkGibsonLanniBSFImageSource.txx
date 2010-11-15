/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniBSFImageSource.cxx,v $
  Language:  C++
  Date:      $Date: 2010/05/24 19:01:23 $
  Version:   $Revision: 1.13 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGibsonLanniBSFImageSource_txx
#define __itkGibsonLanniBSFImageSource_txx

#include "itkGibsonLanniBSFImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"
#include "itkMath.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRotationalExtrusionTransform.h"
#include "itkScanImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkSumProjectionImageFilter.h"

namespace itk
{

/**
 *
 */
template <class TOutputImage>
GibsonLanniBSFImageSource<TOutputImage>
::GibsonLanniBSFImageSource()
{
  m_Size.Fill(32);
  m_Spacing.Fill(65.0);
  m_Origin.Fill(0.0);
  m_IntensityShift = 0.0;
  m_IntensityScale = 1.0;

  m_PSFSource = PSFSourceType::New();

  m_ExtrusionFilter = ExtrusionFilterType::New();
  RotationalExtrusionTransform< double >::Pointer extrusionTransform =
    RotationalExtrusionTransform< double >::New();
  m_ExtrusionFilter->SetTransform(extrusionTransform);
  m_ExtrusionFilter->SetInput(m_PSFSource->GetOutput());

  m_Convolver = ConvolverType::New();
  m_Convolver->SetInput(m_ExtrusionFilter->GetOutput());

  m_RescaleFilter = RescaleImageFilterType::New();
  m_RescaleFilter->SetInput(m_Convolver->GetOutput());
}


template <class TOutputImage>
GibsonLanniBSFImageSource<TOutputImage>
::~GibsonLanniBSFImageSource()
{
}


template <class TOutputImage>
void
GibsonLanniBSFImageSource<TOutputImage>
::SetParameters(const ParametersType& parameters)
{
  int index = 0;
  SpacingType spacing;
  spacing[0] = parameters[index++];
  spacing[1] = parameters[index++];
  spacing[2] = parameters[index++];
  SetSpacing(spacing);

  SetBeadRadius(parameters[index++]);

  PointType center;
  center[0] = parameters[index++];
  center[1] = parameters[index++];
  center[2] = parameters[index++];
  SetBeadCenter(center);

  SetShearX(parameters[index++]);
  SetShearY(parameters[index++]);

  SetEmissionWavelength(parameters[index++]);
  SetNumericalAperture(parameters[index++]);
  SetMagnification(parameters[index++]);

  SetDesignCoverSlipRefractiveIndex(parameters[index++]);
  SetActualCoverSlipRefractiveIndex(parameters[index++]);
  SetDesignCoverSlipThickness(parameters[index++]);
  SetActualCoverSlipThickness(parameters[index++]);
  SetDesignImmersionOilRefractiveIndex(parameters[index++]);
  SetActualImmersionOilRefractiveIndex(parameters[index++]);
  SetDesignImmersionOilThickness(parameters[index++]);

  SetDesignSpecimenLayerRefractiveIndex(parameters[index++]);
  SetActualSpecimenLayerRefractiveIndex(parameters[index++]);
  SetActualPointSourceDepthInSpecimenLayer(parameters[index++]);
  SetDesignDistanceFromBackFocalPlaneToDetector(parameters[index++]);
  SetActualDistanceFromBackFocalPlaneToDetector(parameters[index++]);

  SetIntensityShift(parameters[index++]);
  SetIntensityScale(parameters[index++]);
}


template <class TOutputImage>
typename GibsonLanniBSFImageSource<TOutputImage>::ParametersType
GibsonLanniBSFImageSource<TOutputImage>
::GetParameters() const {
  ParametersType parameters(GetNumberOfParameters());

  int index = 0;

  const SpacingType spacing = GetSpacing();
  parameters[index++] = spacing[0];
  parameters[index++] = spacing[1];
  parameters[index++] = spacing[2];

  parameters[index++] = GetBeadRadius();

  const PointType beadCenter = GetBeadCenter();
  parameters[index++] = beadCenter[0];
  parameters[index++] = beadCenter[1];
  parameters[index++] = beadCenter[2];

  // Shear goes here
  parameters[index++] = GetShearX();
  parameters[index++] = GetShearY();

  parameters[index++] = GetEmissionWavelength();
  parameters[index++] = GetNumericalAperture();
  parameters[index++] = GetMagnification();

  parameters[index++] = GetDesignCoverSlipRefractiveIndex();
  parameters[index++] = GetActualCoverSlipRefractiveIndex();
  parameters[index++] = GetDesignCoverSlipThickness();
  parameters[index++] = GetActualCoverSlipThickness();
  parameters[index++] = GetDesignImmersionOilRefractiveIndex();
  parameters[index++] = GetActualImmersionOilRefractiveIndex();
  parameters[index++] = GetDesignImmersionOilThickness();

  parameters[index++] = GetDesignSpecimenLayerRefractiveIndex();
  parameters[index++] = GetActualSpecimenLayerRefractiveIndex();
  parameters[index++] = GetActualPointSourceDepthInSpecimenLayer();
  parameters[index++] = GetDesignDistanceFromBackFocalPlaneToDetector();
  parameters[index++] = GetActualDistanceFromBackFocalPlaneToDetector();

  parameters[index++] = GetIntensityShift();
  parameters[index++] = GetIntensityScale();

  return parameters;
}


template <class TOutputImage>
unsigned int
GibsonLanniBSFImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return 26;
}


template <class TOutputImage>
void
GibsonLanniBSFImageSource<TOutputImage>
::SetZCoordinate(unsigned int index, double coordinate)
{
  m_Convolver->SetZCoordinate(index, coordinate);
}


template <class TOutputImage>
double
GibsonLanniBSFImageSource<TOutputImage>
::GetZCoordinate(unsigned int index)
{
  return m_Convolver->GetZCoordinate(index);
}


template <class TOutputImage>
void
GibsonLanniBSFImageSource<TOutputImage>
::SetUseCustomZCoordinates(bool use)
{
  m_Convolver->SetUseCustomZCoordinates(use);
  this->Modified();
}


template <class TOutputImage>
bool
GibsonLanniBSFImageSource<TOutputImage>
::GetUseCustomZCoordinates()
{
  return m_Convolver->GetUseCustomZCoordinates();
}


template <class TOutputImage>
void
GibsonLanniBSFImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  unsigned int i;
  os << indent << "Origin: [";
  for ( i=0; i < m_Origin.Size() - 1; i++ )
    {
    os << m_Origin[i] << ", ";
    }
  os << m_Origin[i] << "]" << std::endl;

  os << indent << "Spacing: [";
  for ( i=0; i < m_Spacing.Size() - 1; i++ )
    {
    os << m_Spacing[i] << ", ";
    }
  os << m_Spacing[i] << "] (nanometers)" << std::endl;

  os << indent << "Size: [";
  for ( i=0; i < m_Size.GetSizeDimension() - 1; i++ )
    {
    os << m_Size[i] << ", ";
    }
  os << m_Size[i] << "]" << std::endl;

  os << m_PSFSource << std::endl;
}

//----------------------------------------------------------------------------
template <typename TOutputImage>
void
GibsonLanniBSFImageSource<TOutputImage>
::GenerateOutputInformation()
{
  OutputImageType *output;
  IndexType index = {{0}};
  SizeType size( m_Size );

  output = this->GetOutput(0);

  RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( size );
  largestPossibleRegion.SetIndex( index );
  output->SetLargestPossibleRegion( largestPossibleRegion );

  output->SetSpacing(m_Spacing);
  output->SetOrigin(m_Origin);
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
void
GibsonLanniBSFImageSource<TOutputImage>
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

  m_ExtrusionFilter->SetOutputSpacing(psfTableSpacing);
  m_ExtrusionFilter->SetOutputOrigin(psfTableOrigin);
  m_ExtrusionFilter->SetSize(psfTableSize);

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

  m_PSFSource->SetSize(profileSize);
  m_PSFSource->SetSpacing(profileSpacing);
  m_PSFSource->SetOrigin(profileOrigin);
  m_PSFSource->UpdateLargestPossibleRegion();

  m_ExtrusionFilter->UpdateLargestPossibleRegion();

  m_Convolver->UpdateLargestPossibleRegion();

  m_RescaleFilter->GraftOutput(this->GetOutput());
  m_RescaleFilter->SetShift(m_IntensityShift);
  m_RescaleFilter->SetScale(m_IntensityScale);
  m_RescaleFilter->UpdateLargestPossibleRegion();
  this->GraftOutput(m_RescaleFilter->GetOutput());
}


} // end namespace itk

#endif
