/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSphereConvolutionFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2009/09/14 18:37:18 $
  Version:   $Revision: 1.4 $

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

  m_KernelInterpolator = InterpolatorType::New();
  m_TableInterpolator = InterpolatorType::New();
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
    z1 = z2 = cz;
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

  // Set the inputs for the interpolators
  m_KernelInterpolator->SetInputImage(this->GetInput());
  m_TableInterpolator->SetInputImage(m_ScanImageFilter->GetOutput());
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

    it.Set( ComputeIntegratedVoxelValue(point) );
    progress.CompletedPixel();
    }  
}


template <class TInputImage, class TOutputImage>
float
SphereConvolutionFilter<TInputImage,TOutputImage>
::ComputeSampleValue(OutputImagePointType& point) {
  float value = 0.0f;

  // Compute bounds of geometry sampling region
  float xMin = m_SphereCenter[0] - m_SphereRadius;
  float xMax = m_SphereCenter[0] + m_SphereRadius;
  float yMin = m_SphereCenter[1] - m_SphereRadius;
  float yMax = m_SphereCenter[1] + m_SphereRadius;

  // Compute maximum z value in the pre-integrated table.
  InputImagePointer scannedImage = m_ScanImageFilter->GetOutput();
  float tableZMin = scannedImage->GetOrigin()[2];

  float samplesPerDim = 10.0;
  float diameter = 2.0f*m_SphereRadius;
  float sampleIncrX = diameter / samplesPerDim;
  float sampleIncrY = diameter / samplesPerDim;
  for (float ys = yMin; ys <= yMax; ys += sampleIncrY) {
    for (float xs = xMin; xs <= xMax; xs += sampleIncrX) {

      // Find the intersection z-coordinate values, if they exist.
      float x = point[0];
      float y = point[1];
      float z = point[2];
      float z1 = 0.0f, z2 = 0.0f;
      unsigned int intersections;
      intersections = IntersectWithVerticalLine(xs, ys, z1, z2);
      
      if (intersections == 2) {
	OutputImagePointType p1, p2;
	p1[0] = x - xs;   p1[1] = y - ys;   p1[2] = z - z2;
	p2[0] = x - xs;   p2[1] = y - ys;   p2[2] = z - z1;
	
	// Get values from the pre-integrated table
	bool v1Inside = m_TableInterpolator->IsInsideBuffer(p1);
	bool v2Inside = m_TableInterpolator->IsInsideBuffer(p2);
	InputImagePixelType v1 = 0.0;
	InputImagePixelType v2 = 0.0;
	if (v1Inside) v1 = m_TableInterpolator->Evaluate(p1);
	if (v2Inside) v2 = m_TableInterpolator->Evaluate(p2);

	// If p2 is outside the pre-integrated table, then leaving v2 at 0
	// is fine (the parts of the vertical line that actually contribute
	// will be accounted for.
	//
	// However, if p1 is outside the pre-integrated table, then we need
	// to move p1 to the low z-boundary of the pre-integration table. Leaving
	// v1 = 0.0 is the wrong thing to do.
	if (!v1Inside && v2Inside && p1[2] < tableZMin) {
	  p1[2] = tableZMin + 1e-5;
	  v1 = m_TableInterpolator->Evaluate(p1);
	}
	
	// z1 is always less than z2, and integration goes along positive z,
	// so we return v2 - v1.
	value += v2 - v1;

      } else if (intersections == 1) {
	OutputImagePointType p;
	p[0] = x - xs;   p[1] = y - ys;   p[2] = z - z1;
	
	if (m_KernelInterpolator->IsInsideBuffer(p))
	  value += m_KernelInterpolator->Evaluate(p);
      }

      
    }
  }

  return value;
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
