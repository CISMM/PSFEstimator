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
  m_VoxelSamplesPerDimension = 2;

  m_ScanImageFilter = ScanImageFilterType::New();
  m_ScanImageFilter->SetScanDimension(2);
  m_ScanImageFilter->SetScanOrderToIncreasing();

  m_TableInterpolator = InterpolatorType::New();

  m_LineSampleSpacing = 10; // 10 nm line spacing
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
  double r  = m_SphereRadius;
  double sqrtTerm = (r*r)-(x*x)-(y*y);

  if (sqrtTerm < 0)
    {
    z1 = z2 = 0.0;
    return 0; // no solutions
    }
  else if (sqrtTerm == 0.0)
    {
    z1 = z2 = 0.0;
    return 1; // one solution
    }
  else
    {
    z1 = -sqrt(sqrtTerm);
    z2 =  sqrt(sqrtTerm);
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
::ComputeIntersections()
{
  // Clear the intersection list
  m_IntersectionArray.clear();

  // Add intersection at origin point
  AddIntersection(0.0, 0.0);

  // Iterate over the top left quadrant and use reflections to
  // determine the other points.
  double eps = 1e-6;
  for ( double ys = m_LineSampleSpacing; ys < m_SphereRadius - eps; ys += m_LineSampleSpacing)
    {
    for (double xs = m_LineSampleSpacing; xs < m_SphereRadius - eps; xs += m_LineSampleSpacing)
      {
      AddIntersection( xs,  ys);
      AddIntersection(-xs,  ys);
      AddIntersection( xs, -ys);
      AddIntersection(-xs, -ys);
      }
    }
}


template <class TInputImage, class TOutputImage>
void
SphereConvolutionFilter<TInputImage,TOutputImage>
::AddIntersection(double xs, double ys) {
  // Find the intersection z-coordinate values, if they exist.
  double z1 = 0.0f, z2 = 0.0f;
  unsigned int numIntersections;
  numIntersections = IntersectWithVerticalLine(xs, ys, z1, z2);

  if ( numIntersections > 0 )
    {
    SphereIntersection intersection;
    intersection.x = xs;
    intersection.y = ys;
    intersection.z1 = z1;
    intersection.z2 = z2;
    intersection.numIntersections = numIntersections;
    m_IntersectionArray.push_back(intersection);
    }
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
  m_TableInterpolator->SetInputImage(m_ScanImageFilter->GetOutput());

  // Generate the list of intersections of vertical lines and the
  // sphere.
  ComputeIntersections();
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

  if (m_SphereRadius < 0.0)
    {
    return value;
    }

  // Compute maximum z value in the pre-integrated table.
  InputImagePointer scannedImage = m_ScanImageFilter->GetOutput();
  double tableZMax = (scannedImage->GetSpacing()[2] *
              static_cast<double>(scannedImage->GetLargestPossibleRegion().GetSize()[2]-1)) + scannedImage->GetOrigin()[2];

  IntersectionArrayConstIterator iter;
  for ( IntersectionArrayConstIterator iter = m_IntersectionArray.begin();
        iter != m_IntersectionArray.end();
        iter++)
    {
    SphereIntersection intersection = *iter;
    double xs = intersection.x  + m_SphereCenter[0];
    double ys = intersection.y  + m_SphereCenter[1];
    double z1 = intersection.z1 + m_SphereCenter[2];
    double z2 = intersection.z2 + m_SphereCenter[2];
    int numIntersections = intersection.numIntersections;

    // Find the intersection z-coordinate values, if they exist.
    double x = point[0];
    double y = point[1];
    double z = point[2];

    if (numIntersections == 2)
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
    }

  return value;
}


template <class TInputImage, class TOutputImage>
double
SphereConvolutionFilter<TInputImage,TOutputImage>
::ComputeIntegratedVoxelValue(OutputImagePointType& point)
{
  // Take 9 samples over the CCD element centered at the point
  // parameter.
  double sum = 0.0;

  int numSamples = this->m_VoxelSamplesPerDimension;
  double dx = this->GetSpacing()[0] /
    static_cast<double>(numSamples);
  double dy = this->GetSpacing()[1] /
    static_cast<double>(numSamples);

  OutputImagePointType samplePoint;
  samplePoint[2] = point[2];
  for (int j = 0; j < numSamples; j++)
    {
    samplePoint[1] = point[1]-(0.5*this->GetSpacing()[1]) + (j+0.5)*dy;
    for (int i = 0; i < numSamples; i++)
      {
      samplePoint[0] = point[0]-(0.5*this->GetSpacing()[0]) + (i+0.5)*dx;

      sum += ComputeSampleValue(samplePoint);
      }
    }

  return sum;
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
