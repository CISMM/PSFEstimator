/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSphereConvolutionFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/09/09 20:35:46 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSphereConvolutionFilter_h
#define __itkSphereConvolutionFilter_h

#include "itkImageToImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkScanImageFilter.h"
#include "itkSumProjectionImageFilter.h"

namespace itk
{

/** \class SphereConvolutionFilter
 *
 * \brief Generate an image of a sphere convolved with the input image.
 * Assumes the input image represents  a convolution kernel. The sphere is
 * assumed to define an intensity volume of unit intensity.
 *
 * This class uses the ScanImageFilter to precompute a lookup table for the
 * convolution operation. Each z-plane of the lookup table represents the
 * convolution of the kernel with a line from the z-plane extending to
 * negative infinity. By evaluating where this line intersects the geometry,
 * the contribution from that portion of the line to the current image plane
 * can be quickly computed with two lookups and a subtraction. For greater
 * precision, the input kernel image should be more finely sampled than
 * the output image.
 *
 * \author Cory Quammen. Department of Computer Science, UNC Chapel Hill.
 *
 * \ingroup Multithreaded
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT SphereConvolutionFilter :
  public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef SphereConvolutionFilter   Self;
  typedef ImageSource<TOutputImage> Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Typedef for the output image PixelType. */
  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename InputImageType::IndexType   InputImageIndexType;
  typedef typename InputImageType::SizeType    InputImageSizeType;
  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename InputImageType::PixelType   InputImagePixelType;
  typedef typename InputImageType::PointType   InputImagePointType;

  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::IndexType  OutputImageIndexType;
  typedef typename OutputImageType::SizeType   OutputImageSizeType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;
  typedef typename OutputImageType::PointType  OutputImagePointType;

  typedef Function::SumAccumulator<InputImagePixelType,OutputImagePixelType>
    AccumulatorType;
  typedef ScanImageFilter<InputImageType, InputImageType, AccumulatorType>
    ScanImageFilterType;
  typedef typename ScanImageFilterType::Pointer
    ScanImageFilterPointer;

  typedef LinearInterpolateImageFunction<InputImageType, float>
    InterpolatorType;
  typedef typename InterpolatorType::Pointer
    InterpolatorPointer;

  itkStaticConstMacro(ImageDimension,
		      unsigned int,
		      TOutputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SphereConvolutionFilter,ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Specify the size of the output image. */
  itkSetVectorMacro(Size,unsigned long,TOutputImage::ImageDimension);

  /** Get the size of the output image. */
  itkGetVectorMacro(Size,unsigned long,TOutputImage::ImageDimension);
  
  /** Specify the spacing of the output image. */
  itkSetVectorMacro(Spacing,float,TOutputImage::ImageDimension);

  /** Get the spacing of the output image. */
  itkGetVectorMacro(Spacing,float,TOutputImage::ImageDimension);

  /** Specify the origin of the output image. */
  itkSetVectorMacro(Origin,float,TOutputImage::ImageDimension);

  /** Get the origin of the output image. */
  itkGetVectorMacro(Origin,float,TOutputImage::ImageDimension);

  /** Specify the sphere center. */
  itkSetVectorMacro(SphereCenter,float,TOutputImage::ImageDimension);

  /** Get the sphere center. */
  itkGetVectorMacro(SphereCenter,float,TOutputImage::ImageDimension);

  /** Specify the sphere radius. */
  itkSetMacro(SphereRadius,float);

  /** Get the sphere radius. */
  itkGetMacro(SphereRadius,float);

protected:
  SphereConvolutionFilter();
  ~SphereConvolutionFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  unsigned long *m_Size;         // the number of voxels in each dimension
  float         *m_Spacing;      // the spacing of the voxels
  float         *m_Origin;       // the origin of the image
  float         *m_SphereCenter; // the center of the sphere
  float         m_SphereRadius;  // the radius of the sphere

  ScanImageFilterPointer  m_ScanImageFilter;
  
  InterpolatorPointer     m_KernelInterpolator;
  InterpolatorPointer     m_TableInterpolator;

  /** Gets the z-coordinate(s) of the intersection of a sphere with a line
   * parallel to the z-axis specified by the x- and y-coordinates. z1 and z2
   * are set to the z-coordinates if there are two intersections, only z1 is 
   * set to the z-coordinate if there is one intersection, and neither is set
   * if there is no intersection. The method returns the number of 
   * intersections.
   */
  unsigned int IntersectWithVerticalLine(float x, float y, float& z1, float& z2);

  virtual void GenerateInputRequestedRegion();

  virtual void GenerateOutputInformation();

  virtual void BeforeThreadedGenerateData();

  virtual void ThreadedGenerateData
    (const OutputImageRegionType& outputRegionForThread, int threadId);

  /** Computes the light intensity at a specified point. */
  float ComputeSampleValue(OutputImagePointType& point);

  /** Computes the integrated light intensity over multipe samples per voxel.*/
  float ComputeIntegratedVoxelValue(OutputImagePointType& point);

private:
  SphereConvolutionFilter(const SphereConvolutionFilter&); // purposely not implemented
  void operator=(const SphereConvolutionFilter&); //purposely not implemented

}; // end class SphereConvolutionFilter

} // end namespace itk

#include "itkSphereConvolutionFilter.cxx"

#endif // __itkSphereConvolutionFilter_h
