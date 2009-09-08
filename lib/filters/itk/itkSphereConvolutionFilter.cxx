/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSphereConvolutionFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2009/09/08 21:33:37 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSphereConvolutionFilter_cxx
#define __itkSphereConvolutionFilter_cxx

#include "itkSphereConvolutionFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk {

template <class TInputImage, class TOutputImage>
SphereConvolutionFilter<TInputImage,TOutputImage>
::SphereConvolutionFilter() {
  this->SetNumberOfRequiredInputs(1);

  m_Size    = new unsigned long [TOutputImage::GetImageDimension()];
  m_Spacing = new float [TOutputImage::GetImageDimension()];
  m_Origin  = new float [TOutputImage::GetImageDimension()];
  m_SphereCenter = new float[TOutputImage::GetImageDimension()];

  for (unsigned int i = 0; i < TOutputImage::GetImageDimension(); i++) {
    m_Size[i] = 1;
    m_Spacing[i] = 1.0f;
    m_Origin[i] = 0.0f;
    m_SphereCenter[i] = 0.0f;
  }
  m_SphereRadius = 1.0f;

  m_ScanImageFilter = ScanImageFilterType::New();
  m_ScanImageFilter->SetScanDimension(2);
  m_ScanImageFilter->SetScanOrderToIncreasing();
}


template <class TInputImage, class TOutputImage>
SphereConvolutionFilter<TInputImage,TOutputImage>
::~SphereConvolutionFilter() {
  delete [] m_Size;
  delete [] m_Spacing;
  delete [] m_Origin;
  delete [] m_SphereCenter;
}


template <class TInputImage, class TOutputImage>
unsigned int
SphereConvolutionFilter<TInputImage,TOutputImage>
::IntersectWithVerticalLine(float x, float y, float& z1, float& z2) {
  float cx = m_SphereCenter[0];
  float cy = m_SphereCenter[1];
  float cz = m_SphereCenter[2];
  float r  = m_SphereRadius;
  float sqrtTerm = -(cx*cx)-(cy*cy)+(r*r)+(2*cx*x)-(x*x)+(2*cy*y)-(y*y);
  
  if (sqrtTerm < 0) {
    z1 = z2 = 0.0f;
    return 0; // no solutions
  } else if (sqrtTerm == 0) {
    z1 = z1 = cz;
    return 1; // one solution
  } else {
    z1 = cz - sqrt(sqrtTerm);
    z2 = cz + sqrt(sqrtTerm);
    if (z1 > z2) {
      float tmp = z1;
      z1 = z2;
      z2 = tmp;
    }
    return 2; // two solutions;
  }

}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::GenerateInputRequestedRegion() {
  Superclass::GenerateInputRequestedRegion();

  if (this->GetInput(0)) {
    InputImagePointer input =
      const_cast< TInputImage * > ( this->GetInput(0) );
    InputImageRegionType inputRegion;
    inputRegion = input->GetLargestPossibleRegion();
    input->SetRequestedRegion( inputRegion );
  }
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::GenerateOutputInformation() {
  OutputImageType *output;
  OutputImageIndexType index = {{0}};
  OutputImageSizeType size = {{0}};
  size.SetSize( m_Size );
  
  output = this->GetOutput(0);

  typename TOutputImage::RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( size );
  largestPossibleRegion.SetIndex( index );
  output->SetLargestPossibleRegion( largestPossibleRegion );

  output->SetSpacing(m_Spacing);
  output->SetOrigin(m_Origin);
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::BeforeThreadedGenerateData() {
  // Compute the scan of the convolution kernel.
  m_ScanImageFilter->SetInput(this->GetInput());
  m_ScanImageFilter->Update();
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage> 
::ThreadedGenerateData
(const OutputImageRegionType& outputRegionForThread, int threadId) {
  itkDebugMacro(<<"Generating a random image of scalars");

  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  OutputImagePointer image = this->GetOutput(0);

  ImageRegionIteratorWithIndex<TOutputImage> it(image, outputRegionForThread);

  for (; !it.IsAtEnd(); ++it)
    {
    OutputImageIndexType index = it.GetIndex();
    OutputImagePointType point;
    image->TransformIndexToPhysicalPoint(index, point);
    std::cout << "SphereConvFilter: " << point[0] << ", " << point[1]
	      << ", " << point[2] << std::endl;

    it.Set( ComputeIntegratedVoxelValue(point) );
    progress.CompletedPixel();
    }  
}


template <class TInputImage, class TOutputImage>
float
SphereConvolutionFilter<TInputImage,TOutputImage>
::ComputeSampleValue(OutputImagePointType& point) {
  float value = 0.0f;

  // Find the intersection z-coordinate values, if they exist.
  unsigned int intersections;
  float x = point[0];
  float y = point[1];
  float z1 = 0.0f, z2 = 0.0f;
  intersections = IntersectWithVerticalLine(x, y, z1, z2);

  if (intersections == 2) {
    OutputImagePointType p1, p2;
    p1[0] = x;   p1[1] = y;   p1[2] = z1;
    p2[0] = x;   p2[1] = y;   p2[2] = z2;

    // Get values from the pre-integrated table
    InputImageIndexType index;
    InputImagePointer lookupTable = m_ScanImageFilter->GetOutput();

    lookupTable->TransformPhysicalPointToIndex(p1, index);
    InputImagePixelType v1 = lookupTable->GetPixel(index);
    lookupTable->TransformPhysicalPointToIndex(p2, index);
    InputImagePixelType v2 = lookupTable->GetPixel(index);

    // z1 is always less than z2, and integration goes along positive z,
    // so we return v2 - v1.
    value = v2 - v1;

  } else if (intersections == 1) {
    OutputImagePointType p;
    p[0] = x;   p[1] = y;  p[2] = z1;

    // Get value from the single intersection in the original kernel.
    // An interpolator could go here.
    InputImageIndexType index;
    this->GetInput()->TransformPhysicalPointToIndex(p, index);
    value = this->GetInput()->GetPixel(index);
  }

  return value; // No contribution
}


template <class TInputImage, class TOutputImage>
float
SphereConvolutionFilter<TInputImage,TOutputImage>
::ComputeIntegratedVoxelValue(OutputImagePointType& point) {
  // TODO - integrate over voxel xy-plane area
  return ComputeSampleValue(point);
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const {
  Superclass::PrintSelf(os, indent);
  unsigned int i;
  os << indent << "Origin: [";
  for (i = 0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Origin[i] << ", ";
    }
  os << m_Origin[i] << "]" << std::endl;

  os << indent << "Spacing: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Spacing[i] << ", ";
    }
  os << m_Spacing[i] << "]" << std::endl;

  os << indent << "Size: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Size[i] << ", ";
    }
  os << m_Size[i] << "]" << std::endl;

  os << indent << "SphereCenter: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_SphereCenter[i] << ", ";
    }
  os << m_SphereCenter[i] << "]" << std::endl;

  os << indent << "SphereRadius: " << m_SphereRadius << std::endl;

}


} // end namespace itk

#endif // __itkSphereConvolutionFilter_cxx
