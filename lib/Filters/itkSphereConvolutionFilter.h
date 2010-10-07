/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSphereConvolutionFilter.h,v $
  Language:  C++
  Date:      $Date: 2010/03/30 03:58:43 $
  Version:   $Revision: 1.6 $

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
 * This filter assumes the input image represents a convolution
 * kernel. The sphere is assumed to define an intensity volume of unit
 * intensity.
 *
 * This class uses the ScanImageFilter to precompute a lookup table for the
 * convolution operation. Each z-plane of the lookup table represents the
 * convolution of the kernel with a line from the z-plane extending to
 * negative infinity. By evaluating where this line intersects the geometry,
 * the contribution from that portion of the line to the current image plane
 * can be quickly computed with two lookups and a subtraction. For greater
 * accuracy, the input kernel image should be more finely sampled than
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

  typedef TOutputImage                          OutputImageType;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename OutputImageType::IndexType   OutputImageIndexType;
  typedef typename OutputImageType::SizeType    OutputImageSizeType;
  typedef typename OutputImageType::RegionType  OutputImageRegionType;
  typedef typename OutputImageType::PixelType   OutputImagePixelType;
  typedef typename OutputImageType::SpacingType OutputImageSpacingType;
  typedef typename OutputImageType::PointType   OutputImagePointType;

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

  itkStaticConstMacro(ImageDimension, unsigned int,
		      TOutputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SphereConvolutionFilter, ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Specify the size of the output image. */
  virtual void SetSize(const OutputImageSizeType & size)
  {
    if (size != m_Size)
      {
      m_Size = size;
      this->Modified();
      }
    m_ZCoordinate.resize(size[TOutputImage::ImageDimension-1]);
  }

  /** Get the size of the output image. */
  itkGetConstReferenceMacro(Size, OutputImageSizeType);

  /** Specify the spacing of the output image. */
  virtual void SetSpacing(const OutputImageSpacingType & spacing)
  {
    if (spacing != m_Spacing)
      {
      m_Spacing = spacing;
      this->Modified();
      }
  }

  /** Get the spacing of the output image. */
  itkGetConstReferenceMacro(Spacing, OutputImageSpacingType);

  /** Specify the origin of the output image. */
  virtual void SetOrigin(const OutputImagePointType & origin)
  {
    if (origin != m_Origin)
      {
      m_Origin = origin;
      this->Modified();
      }
  }

  /** Get the origin of the output image. */
  itkGetConstReferenceMacro(Origin, OutputImagePointType);

  /** Specify the sphere center. */
  virtual void SetSphereCenter(const OutputImagePointType & center)
  {
    if (center != m_SphereCenter)
      {
      m_SphereCenter = center;
      this->Modified();
      }
  }

  /** Get the sphere center. */
  itkGetConstReferenceMacro(SphereCenter, OutputImagePointType);

  /** Specify the sphere radius. */
  itkSetMacro(SphereRadius, double);

  /** Get the sphere radius. */
  itkGetMacro(SphereRadius, double);

  /** Specify the shear in the X direction. */
  itkSetMacro(ShearX, double);

  /** Get the shear in the X direction. */
  itkGetMacro(ShearX, double);

  /** Specify the shear in the Y direction. */
  itkSetMacro(ShearY, double);

  /** Get the shear in the Y direction. */
  itkGetMacro(ShearY, double);

  /** Get/set the z-coordinate of the image z-plane at the given index. */
  void SetZCoordinate(unsigned int index, double coordinate);
  double GetZCoordinate(unsigned int index);

    /** Get/set use of custom z coordinates. */
  itkSetMacro(UseCustomZCoordinates, bool);
  itkGetMacro(UseCustomZCoordinates, bool);

protected:
  SphereConvolutionFilter();
  ~SphereConvolutionFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;

  OutputImageSizeType    m_Size;         // the number of voxels in each dimension
  OutputImageSpacingType m_Spacing;      // the spacing of the voxels
  OutputImagePointType   m_Origin;       // the origin of the image
  OutputImagePointType   m_SphereCenter; // the center of the sphere
  double                 m_SphereRadius; // the radius of the sphere

  double              m_ShearX;      // shear in the x direction w.r.t. z
  double              m_ShearY;      // shear in the y direction w.r.t. z
  std::vector<double> m_ZCoordinate; // z-slice coordinates
  bool                m_UseCustomZCoordinates;

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
  unsigned int IntersectWithVerticalLine(double x, double y, double& z1, double& z2);

  virtual void GenerateInputRequestedRegion();

  virtual void GenerateOutputInformation();

  virtual void BeforeThreadedGenerateData();

  virtual void ThreadedGenerateData
    (const OutputImageRegionType& outputRegionForThread, int threadId);

  /** Computes the light intensity at a specified point. */
  double ComputeSampleValue(OutputImagePointType& point);

  /** Computes the integrated light intensity over multipe samples per voxel.*/
  double ComputeIntegratedVoxelValue(OutputImagePointType& point);

private:
  SphereConvolutionFilter(const SphereConvolutionFilter&); // purposely not implemented
  void operator=(const SphereConvolutionFilter&); //purposely not implemented

}; // end class SphereConvolutionFilter
} // end namespace itk

#include "itkSphereConvolutionFilter.txx"

#endif // __itkSphereConvolutionFilter_h
