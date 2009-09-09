/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniBSFImageSource.cxx,v $
  Language:  C++
  Date:      $Date: 2009/09/09 20:57:09 $
  Version:   $Revision: 1.3 $

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
  for (unsigned int i=0; i<TOutputImage::GetImageDimension(); i++)
    {
    m_Size[i] = 63;
    m_Spacing[i] = 65.0f;
    m_Origin[i] = 0.0f;
    m_BeadCenter[i] = 0.0f;
    }

  m_BeadRadius = 200.0;
  m_PSFSource = PSFSourceType::New();
  m_Convolver = ConvolverType::New();
  m_Convolver->SetInput(m_PSFSource->GetOutput());
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
::SetParameters(const ParametersType& parameters) {
  Array<float> floatParams(GetNumberOfParameters());
  for (unsigned int i = 0; i < GetNumberOfParameters(); i++) {
    floatParams[i] = static_cast<float>(parameters[i]);
  }

  int index = 0;
  float spacing[3];
  spacing[0] = floatParams[index++];
  spacing[1] = floatParams[index++];
  spacing[2] = floatParams[index++];
  SetSpacing(spacing);

  float center[3];
  center[0] = floatParams[index++];
  center[1] = floatParams[index++];
  center[2] = floatParams[index++];
  SetBeadCenter(center);

  SetBeadRadius(floatParams[index++]);

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

  float* beadCenter = GetBeadCenter();
  floatParams[index++] = beadCenter[0];
  floatParams[index++] = beadCenter[1];
  floatParams[index++] = beadCenter[2];

  floatParams[index++] = GetBeadRadius();

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

  ParametersType parameters(GetNumberOfParameters());
  for (unsigned int i = 0; i < GetNumberOfParameters(); i++) {
    parameters[i] = static_cast<double>(floatParams[i]);
  }

  return parameters;
}


template <class TOutputImage>
unsigned int
GibsonLanniBSFImageSource<TOutputImage>
::GetNumberOfParameters() const {
  return 22;
}


template <class TOutputImage>
void 
GibsonLanniBSFImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  unsigned int i;
  os << indent << "Origin: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Origin[i] << ", ";
    }
  os << m_Origin[i] << "]" << std::endl;

  os << indent << "Spacing: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Spacing[i] << ", ";
    }
  os << m_Spacing[i] << "] (nanometers)" << std::endl;

  os << indent << "Size: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Size[i] << ", ";
    }
  os << m_Size[i] << "]" << std::endl;

  os << indent << "BeadCenter: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
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
::GenerateData() {

  // Update the PSF source. Let's super sample it by 2 times in x and  y
  // and 4 times in z and dilated by at least the sphere radius.
  unsigned long superX = 2;
  unsigned long superY = 2;
  unsigned long superZ = 4;

  // Set the PSF size parameters and update.
  unsigned long psfSize[3];
  psfSize[0] = GetSize()[0] * superX;
  psfSize[1] = GetSize()[1] * superY;
  psfSize[2] = GetSize()[2] * superZ;

  //m_PSFSource->SetSize(this->GetSize());
  m_PSFSource->SetSize(psfSize);

  float psfSpacing[3];
  psfSpacing[0] = GetSpacing()[0] / static_cast<float>(superX);
  psfSpacing[1] = GetSpacing()[1] / static_cast<float>(superY);
  psfSpacing[2] = GetSpacing()[2] / static_cast<float>(superZ);
  m_PSFSource->SetSpacing(psfSpacing);

  m_PSFSource->SetOrigin(this->GetOrigin());
  m_PSFSource->Update();

  m_Convolver->GraftOutput(this->GetOutput());
  m_Convolver->Update();
  this->GraftOutput(m_Convolver->GetOutput());

}


} // end namespace itk

#endif
