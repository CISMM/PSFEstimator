/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSphereConvolutionFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2010/03/30 02:53:22 $
  Version:   $Revision: 1.12 $

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
::SphereConvolutionFilter()
{
  this->SetNumberOfRequiredInputs(1);

  m_Size.Fill(1);
  m_Spacing.Fill(1.0);
  m_Origin.Fill(0.0);
  m_SphereCenter.Fill(0.0);
  m_SphereRadius = 1.0f;
  m_ShearX = 0.0f;
  m_ShearY = 0.0f;
  m_UseCustomZCoordinates = false;

  m_ScanImageFilter = ScanImageFilterType::New();
  m_ScanImageFilter->SetScanDimension(2);
  m_ScanImageFilter->SetScanOrderToIncreasing();

  m_KernelInterpolator = InterpolatorType::New();
  m_TableInterpolator = InterpolatorType::New();
}


template <class TInputImage, class TOutputImage>
SphereConvolutionFilter<TInputImage,TOutputImage>
::~SphereConvolutionFilter()
{
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::SetZCoordinate(unsigned int index, double coordinate)
{
  if (index >= m_ZCoordinate.size())
    {
    m_ZCoordinate.resize(index+1);
    }

  m_ZCoordinate[index] = coordinate;
  this->Modified();
}


template <class TInputImage, class TOutputImage>
double
SphereConvolutionFilter<TInputImage,TOutputImage>
::GetZCoordinate(unsigned int index)
{
  if (index < m_ZCoordinate.size())
    {
    return m_ZCoordinate[index];
    }

  return 0.0;
}


template <class TInputImage, class TOutputImage>
unsigned int
SphereConvolutionFilter<TInputImage,TOutputImage>
::IntersectWithVerticalLine(double x, double y, double& z1, double& z2)
{
  double cx = m_SphereCenter[0];
  double cy = m_SphereCenter[1];
  double cz = m_SphereCenter[2];
  double r  = m_SphereRadius;
  double sqrtTerm = -(cx*cx)-(cy*cy)+(r*r)+(2*cx*x)-(x*x)+(2*cy*y)-(y*y);

  if (sqrtTerm < 0)
    {
    z1 = z2 = 0.0f;
    return 0; // no solutions
    }
  else if (sqrtTerm == 0)
    {
    z1 = z2 = cz;
    return 1; // one solution
    }
  else
    {
    z1 = cz - sqrt(sqrtTerm);
    z2 = cz + sqrt(sqrtTerm);
    if (z1 > z2)
      {
      double tmp = z1;
      z1 = z2;
      z2 = tmp;
      }
    return 2; // two solutions;
    }
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  if (this->GetInput(0))
    {
    InputImagePointer input = const_cast< TInputImage * > ( this->GetInput(0) );
    InputImageRegionType inputRegion;
    inputRegion = input->GetLargestPossibleRegion();
    input->SetRequestedRegion( inputRegion );
    }
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::GenerateOutputInformation()
{
  OutputImageType *output;
  OutputImageIndexType index = {{0}};
  OutputImageSizeType size(m_Size);

  output = this->GetOutput(0);

  typename TOutputImage::RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize(size);
  largestPossibleRegion.SetIndex(index);
  output->SetLargestPossibleRegion(largestPossibleRegion);

  output->SetSpacing(m_Spacing);
  output->SetOrigin(m_Origin);
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::BeforeThreadedGenerateData()
{
  // Compute the scan of the convolution kernel.
  m_ScanImageFilter->SetInput(this->GetInput());
  m_ScanImageFilter->UpdateLargestPossibleRegion();

  // Set the inputs for the interpolators
  m_KernelInterpolator->SetInputImage(this->GetInput());
  m_TableInterpolator->SetInputImage(m_ScanImageFilter->GetOutput());
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::ThreadedGenerateData
(const OutputImageRegionType& outputRegionForThread, int threadId)
{
  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  OutputImagePointer image = this->GetOutput(0);

  ImageRegionIteratorWithIndex<TOutputImage> it(image, outputRegionForThread);

  for (; !it.IsAtEnd(); ++it)
    {
    OutputImageIndexType index = it.GetIndex();
    OutputImagePointType point;
    image->TransformIndexToPhysicalPoint(index, point);

    // Change the z coordinate here if using custom z coordinates
    if (m_UseCustomZCoordinates)
      {
      point[2] = GetZCoordinate(index[2]);
      }

    // Apply shear here
    point[0] -= m_ShearX * (point[2] - m_SphereCenter[2]);
    point[1] -= m_ShearY * (point[2] - m_SphereCenter[2]);

    it.Set( ComputeIntegratedVoxelValue(point) );
    progress.CompletedPixel();
    }
}


template <class TInputImage, class TOutputImage>
double
SphereConvolutionFilter<TInputImage,TOutputImage>
::ComputeSampleValue(OutputImagePointType& point)
{
  double value = 0.0f;

  if (m_SphereRadius <= 0.0)
    {
    return value;
    }

  // Compute bounds of geometry sampling region
  double xMin = m_SphereCenter[0] - m_SphereRadius;
  double xMax = m_SphereCenter[0] + m_SphereRadius;
  double yMin = m_SphereCenter[1] - m_SphereRadius;
  double yMax = m_SphereCenter[1] + m_SphereRadius;

  // Compute maximum z value in the pre-integrated table.
  InputImagePointer scannedImage = m_ScanImageFilter->GetOutput();
  double tableZMax = (scannedImage->GetSpacing()[2] *
              static_cast<double>(scannedImage->GetLargestPossibleRegion().GetSize()[2]-1)) + scannedImage->GetOrigin()[2];

  double samplesPerDim = 10.0;
  double diameter = 2.0f*m_SphereRadius;
  double sampleIncrX = diameter / samplesPerDim;
  double sampleIncrY = diameter / samplesPerDim;
  for ( double  ys = yMin; ys <= yMax; ys += sampleIncrY )
    {
    for ( double xs = xMin; xs <= xMax; xs += sampleIncrX )
      {

      // Find the intersection z-coordinate values, if they exist.
      double x = point[0];
      double y = point[1];
      double z = point[2];
      double z1 = 0.0f, z2 = 0.0f;
      unsigned int intersections;
      intersections = IntersectWithVerticalLine(xs, ys, z1, z2);

      if (intersections == 2)
        {
        OutputImagePointType p1, p2;
        p1[0] = x - xs;   p1[1] = y - ys;   p1[2] = z - z1;
        p2[0] = x - xs;   p2[1] = y - ys;   p2[2] = z - z2;

        // Important: z1 is always less than z2, so p1 is always above p2

        // Subtract z-voxel spacing from p2[2] to get the proper behavior
        // in the pre-integrated PSF table.
        p2[2] -= scannedImage->GetSpacing()[2];

        // Get values from the pre-integrated table
        bool v1Inside = m_TableInterpolator->IsInsideBuffer(p1);
        bool v2Inside = m_TableInterpolator->IsInsideBuffer(p2);
        InputImagePixelType v1 = 0.0;
        InputImagePixelType v2 = 0.0;
        if (v1Inside) v1 = m_TableInterpolator->Evaluate(p1);
        if (v2Inside) v2 = m_TableInterpolator->Evaluate(p2);

        if (!v1Inside && v2Inside && p1[2] > tableZMax)
          {
          p1[2] = tableZMax - 1e-5;
          v1 = m_TableInterpolator->Evaluate(p1);
          }
        // z - z1 is always larger than z - z2, and integration goes along
        // positive z, so we return v1 - v2.
        value += v1 - v2;
        }
      else if (intersections == 1)
        {
        OutputImagePointType p;
        p[0] = x - xs;   p[1] = y - ys;   p[2] = z - z1;

        if (m_KernelInterpolator->IsInsideBuffer(p))
          {
          value += m_KernelInterpolator->Evaluate(p);
          }
        }
      }
    }

  return value;
}


template <class TInputImage, class TOutputImage>
double
SphereConvolutionFilter<TInputImage,TOutputImage>
::ComputeIntegratedVoxelValue(OutputImagePointType& point)
{
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
  for ( i = 0; i < m_Origin.Size() - 1; i++ )
    {
    os << m_Origin[i] << ", ";
    }
  os << m_Origin[i] << "]" << std::endl;

  os << indent << "Spacing: [";
  for ( i=0; i < m_Spacing.Size() - 1; i++ )
    {
    os << m_Spacing[i] << ", ";
    }
  os << m_Spacing[i] << "]" << std::endl;

  os << indent << "Size: [";
  for ( i=0; i < m_Size.GetSizeDimension() - 1; i++ )
    {
    os << m_Size[i] << ", ";
    }
  os << m_Size[i] << "]" << std::endl;

  os << indent << "SphereCenter: [";
  for ( i=0; i < m_SphereCenter.Size() - 1; i++ )
    {
    os << m_SphereCenter[i] << ", ";
    }
  os << m_SphereCenter[i] << "]" << std::endl;

  os << indent << "SphereRadius: " << m_SphereRadius << std::endl;

}


} // end namespace itk

#endif // __itkSphereConvolutionFilter_cxx
