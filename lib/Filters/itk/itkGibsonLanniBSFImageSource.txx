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
#include "itkGaussianImageSource.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScanImageFilter.h"
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
  m_Size       = new unsigned long [TOutputImage::GetImageDimension()];
  m_Spacing    = new float [TOutputImage::GetImageDimension()];
  m_Origin     = new float [TOutputImage::GetImageDimension()];
  m_BeadCenter = new float[TOutputImage::GetImageDimension()];

  //Initial image is 63 wide in each direction.
  for ( unsigned int i = 0; i<TOutputImage::GetImageDimension(); i++ )
    {
    m_Size[i] = 63;
    m_Spacing[i] = 65.0f;
    m_Origin[i] = 0.0f;
    m_BeadCenter[i] = 0.0f;
    }

  m_BeadRadius = 200.0;

  m_BackgroundIntensity = 0.0;
  m_MaximumIntensity = 1.0;

  m_PSFSource = PSFSourceType::New();
  m_Convolver = ConvolverType::New();
  m_Convolver->SetInput(m_PSFSource->GetOutput());

  m_RescaleFilter = RescaleImageFilterType::New();
  m_RescaleFilter->SetInput(m_Convolver->GetOutput());
}


template <class TOutputImage>
GibsonLanniBSFImageSource<TOutputImage>
::~GibsonLanniBSFImageSource()
{
  delete [] m_Size;
  delete [] m_Spacing;
  delete [] m_Origin;
  delete [] m_BeadCenter;
}


template <class TOutputImage>
void
GibsonLanniBSFImageSource<TOutputImage>
::SetParameters(const ParametersType& parameters)
{
  Array<float> floatParams(GetNumberOfParameters());
  for ( unsigned int i = 0; i < GetNumberOfParameters(); i++ )
    {
    floatParams[i] = static_cast<float>(parameters[i]);
    }

  int index = 0;
  float spacing[3];
  spacing[0] = floatParams[index++];
  spacing[1] = floatParams[index++];
  spacing[2] = floatParams[index++];
  SetSpacing(spacing);

  // CCD border goes here
  index += 2;

  SetBeadRadius(floatParams[index++]);

  float center[3];
  center[0] = floatParams[index++];
  center[1] = floatParams[index++];
  center[2] = floatParams[index++];
  SetBeadCenter(center);

  SetShearX(floatParams[index++]);
  SetShearY(floatParams[index++]);

  SetEmissionWavelength(floatParams[index++]);
  SetNumericalAperture(floatParams[index++]);
  SetMagnification(floatParams[index++]);

  SetDesignCoverSlipRefractiveIndex(floatParams[index++]);
  SetActualCoverSlipRefractiveIndex(floatParams[index++]);
  SetDesignCoverSlipThickness(floatParams[index++]);
  SetActualCoverSlipThickness(floatParams[index++]);
  SetDesignImmersionOilRefractiveIndex(floatParams[index++]);
  SetActualImmersionOilRefractiveIndex(floatParams[index++]);
  SetDesignImmersionOilThickness(floatParams[index++]);

  SetDesignSpecimenLayerRefractiveIndex(floatParams[index++]);
  SetActualSpecimenLayerRefractiveIndex(floatParams[index++]);
  SetActualPointSourceDepthInSpecimenLayer(floatParams[index++]);
  SetDesignDistanceFromBackFocalPlaneToDetector(floatParams[index++]);
  SetActualDistanceFromBackFocalPlaneToDetector(floatParams[index++]);

  SetBackgroundIntensity(floatParams[index++]);
  SetMaximumIntensity(floatParams[index++]);
}


template <class TOutputImage>
typename GibsonLanniBSFImageSource<TOutputImage>::ParametersType
GibsonLanniBSFImageSource<TOutputImage>
::GetParameters() const {
  Array<float> floatParams(GetNumberOfParameters());

  int index = 0;
  floatParams[index++] = GetSpacing()[0];
  floatParams[index++] = GetSpacing()[1];
  floatParams[index++] = GetSpacing()[2];

  // CCD border goes here
  floatParams[index++] = 0.0;
  floatParams[index++] = 0.0;

  floatParams[index++] = GetBeadRadius();

  float* beadCenter = GetBeadCenter();
  floatParams[index++] = beadCenter[0];
  floatParams[index++] = beadCenter[1];
  floatParams[index++] = beadCenter[2];

  // Shear goes here
  floatParams[index++] = GetShearX();
  floatParams[index++] = GetShearY();

  floatParams[index++] = GetEmissionWavelength();
  floatParams[index++] = GetNumericalAperture();
  floatParams[index++] = GetMagnification();

  floatParams[index++] = GetDesignCoverSlipRefractiveIndex();
  floatParams[index++] = GetActualCoverSlipRefractiveIndex();
  floatParams[index++] = GetDesignCoverSlipThickness();
  floatParams[index++] = GetActualCoverSlipThickness();
  floatParams[index++] = GetDesignImmersionOilRefractiveIndex();
  floatParams[index++] = GetActualImmersionOilRefractiveIndex();
  floatParams[index++] = GetDesignImmersionOilThickness();

  floatParams[index++] = GetDesignSpecimenLayerRefractiveIndex();
  floatParams[index++] = GetActualSpecimenLayerRefractiveIndex();
  floatParams[index++] = GetActualPointSourceDepthInSpecimenLayer();
  floatParams[index++] = GetDesignDistanceFromBackFocalPlaneToDetector();
  floatParams[index++] = GetActualDistanceFromBackFocalPlaneToDetector();

  floatParams[index++] = GetBackgroundIntensity();
  floatParams[index++] = GetMaximumIntensity();

  ParametersType parameters(GetNumberOfParameters());
  for (unsigned int i = 0; i < GetNumberOfParameters(); i++)
    {
    parameters[i] = static_cast<double>(floatParams[i]);
    }

  return parameters;
}


template <class TOutputImage>
unsigned int
GibsonLanniBSFImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return 28;
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
  for ( i=0; i < TOutputImage::ImageDimension; i++ )
    {
    os << m_Origin[i] << ", ";
    }
  os << m_Origin[i] << "]" << std::endl;

  os << indent << "Spacing: [";
  for ( i=0; i < TOutputImage::ImageDimension; i++ )
    {
    os << m_Spacing[i] << ", ";
    }
  os << m_Spacing[i] << "] (nanometers)" << std::endl;

  os << indent << "Size: [";
  for ( i=0; i < TOutputImage::ImageDimension; i++ )
    {
    os << m_Size[i] << ", ";
    }
  os << m_Size[i] << "]" << std::endl;

  os << indent << "BeadCenter: [";
  for ( i=0; i < TOutputImage::ImageDimension; i++ )
    {
    os << m_BeadCenter[i] << ", ";
    }
  os << m_BeadCenter[i] << "]" << std::endl;

  os << m_PSFSource << std::endl;
}

//----------------------------------------------------------------------------
template <typename TOutputImage>
void
GibsonLanniBSFImageSource<TOutputImage>
::GenerateOutputInformation()
{
  TOutputImage *output;
  typename TOutputImage::IndexType index = {{0}};
  typename TOutputImage::SizeType size = {{0}};
  size.SetSize( m_Size );

  output = this->GetOutput(0);

  typename TOutputImage::RegionType largestPossibleRegion;
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
  float psfSpacing[3], psfOrigin[3];

  // Determine Nyquist sampling (taken from Heintzmann, R. and Sheppard, C.
  // (2007). The sampling limit in fluorescence microscopy. Micron,
  // 38(2):145–149. Actually, sample at twice this rate.
  float NA     = GetNumericalAperture();
  float lambda = GetEmissionWavelength();
  float n      = GetDesignImmersionOilRefractiveIndex();
  double alpha = asin(NA/n);
  psfSpacing[0] = psfSpacing[1] = 0.5*lambda / (4.0 * NA);
  psfSpacing[2] = 0.5*(2.0 * psfSpacing[0] * sin(alpha)) / (1.0 - cos(alpha));

  // Determine necessary spatial extent of PSF table.
  unsigned long psfSize[3];
  for (int i = 0; i < 3; i++) {
    // First calculate extent of BSF in this dimension.
    float minExtent = GetOrigin()[i];
    float maxExtent = static_cast<float>(GetSize()[i]-1)*GetSpacing()[i] +
      GetOrigin()[i];

    // Now modify calculated PSF dimension to account for bead shift and radius
    minExtent += -GetBeadCenter()[i] - m_BeadRadius;
    maxExtent += -GetBeadCenter()[i] + m_BeadRadius;

    // Determine logical extent of the PSF table for the min and max extents.
    long iDimMin = static_cast<long>(floor(minExtent / psfSpacing[i]));
    psfOrigin[i] = static_cast<float>(iDimMin)*psfSpacing[i];
    long iDimMax = static_cast<long>(ceil(maxExtent / psfSpacing[i]));

    // Determine the logical extent of the PSF table in this dimension.
    psfSize[i] = iDimMax - iDimMin + 1;
  }

  m_PSFSource->SetSize(psfSize);
  m_PSFSource->SetSpacing(psfSpacing);
  m_PSFSource->SetOrigin(psfOrigin);
  m_PSFSource->UpdateLargestPossibleRegion();

  m_Convolver->UpdateLargestPossibleRegion();

  m_RescaleFilter->GraftOutput(this->GetOutput());
  m_RescaleFilter->SetOutputMinimum(m_BackgroundIntensity);
  m_RescaleFilter->SetOutputMaximum(m_MaximumIntensity);
  m_RescaleFilter->UpdateLargestPossibleRegion();
  this->GraftOutput(m_RescaleFilter->GetOutput());
}


} // end namespace itk

#endif
