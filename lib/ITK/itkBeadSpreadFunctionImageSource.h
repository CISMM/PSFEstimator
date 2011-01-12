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
#ifndef _itkBeadSpreadFunctionImageSource_h
#define _itkBeadSpreadFunctionImageSource_h

#include "itkParametricImageSource.h"
#include "itkShiftScaleImageFilter.h"
#include "itkSphereConvolutionFilter.h"

namespace itk
{

/** \class BeadSpreadFunctionImageSource
 *
 * \brief Generates a synthetic bead-spread function that is the
 * convolution of a sphere with a ParametricImageSource that generates
 * a convolution kernel.
 *
 * \ingroup DataSources Multithreaded
*/
template < class TOutputImage >
class ITK_EXPORT BeadSpreadFunctionImageSource :
    public ParametricImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef BeadSpreadFunctionImageSource         Self;
  typedef ParametricImageSource< TOutputImage > Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /** Typedef for output types. */
  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::PixelType      PixelType;
  typedef typename OutputImageType::RegionType     RegionType;
  typedef typename OutputImageType::PointType      PointType;
  typedef typename OutputImageType::PointValueType PointValueType;
  typedef typename PointType::VectorType           VectorType;
  typedef typename OutputImageType::SpacingType    SpacingType;
  typedef typename OutputImageType::IndexType      IndexType;
  typedef typename OutputImageType::SizeType       SizeType;
  typedef typename OutputImageType::SizeValueType  SizeValueType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef ParametricImageSource< TOutputImage >
    KernelImageSourceType;
  typedef typename KernelImageSourceType::Pointer
    KernelImageSourcePointer;
  typedef SphereConvolutionFilter< TOutputImage, TOutputImage >
    ConvolverType;
  typedef typename ConvolverType::Pointer
    ConvolverPointer;
  typedef ShiftScaleImageFilter< TOutputImage, TOutputImage >
    RescaleImageFilterType;
  typedef typename RescaleImageFilterType::Pointer
    RescaleImageFilterPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BeadSpreadFunctionImageSource,ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::ParametersType      ParametersType;

  /** Specify the size of the output image. */
  void SetSize(const SizeType & size)
  {
    if (size != m_Convolver->GetSize())
      {
      this->Modified();
      }
    m_Convolver->SetSize(size);
  }

  /** Get the size of the output image. */
  const SizeType & GetSize() const
  {
    return m_Convolver->GetSize();
  }

  /** Specify the spacing of the output image (in nanometers). */
  void SetSpacing(const SpacingType & spacing)
  {
    if (spacing != m_Convolver->GetSpacing())
      {
      this->Modified();
      }
    m_Convolver->SetSpacing(spacing);
  }

  /** Get the spacing of the output image (in nanometers). */
  const SpacingType & GetSpacing() const
  {
    return m_Convolver->GetSpacing();
  }

  /** Specify the origin of the output image (in nanometers). */
  virtual void SetOrigin(const PointType & origin)
  {
    if (origin != m_Convolver->GetOrigin())
      {
      this->Modified();
      }
    m_Convolver->SetOrigin(origin);
  }

  /** Get the origin of the output image (in nanometers). */
  const PointType & GetOrigin() const
  {
    return m_Convolver->GetOrigin();
  }

  /** Specify the point source center (in nanometers). */
  virtual void SetBeadCenter(const PointType & center)
  {
    if (center != m_Convolver->GetSphereCenter())
      {
      m_Convolver->SetSphereCenter(center);
      this->Modified();
      }
  }

  /** Get the point source center (in nanometers). */
  const PointType & GetBeadCenter() const
  {
    return m_Convolver->GetSphereCenter();
  }

  /** Specify the bead radius (in nanometers). */
  void SetBeadRadius(double radius)
  {
    if (radius != m_Convolver->GetSphereRadius())
      {
      m_Convolver->SetSphereRadius(radius);
      this->Modified();
      }
  }

  /** Get the bead radius. */
  double GetBeadRadius() const
  {
    return m_Convolver->GetSphereRadius();
  }

  /** Specify the shear in the X direction. */
  void SetShearX(double shear)
  {
    if (shear != m_Convolver->GetShearX())
      {
      m_Convolver->SetShearX(shear);
      this->Modified();
      }
  }

  /** Get the shear in the X direction. */
  double GetShearX() const
  {
    return m_Convolver->GetShearX();
  }

  /** Specify the shear in the Y direction. */
  void SetShearY(double shear)
  {
    if (shear != m_Convolver->GetShearY())
      {
      m_Convolver->SetShearY(shear);
      this->Modified();
      }
  }

  /** Get the shear in the Y direction. */
  double GetShearY() const
  {
    return m_Convolver->GetShearY();
  }

  /** Set/get the background value. */
  itkSetMacro(IntensityShift, double);
  itkGetConstMacro(IntensityShift, double);

  /** Set/get the maximum intensity. */
  itkSetMacro(IntensityScale, double);
  itkGetConstMacro(IntensityScale, double);

  /** Set/get the convolution kernel source. */
  itkSetObjectMacro(KernelSource, KernelImageSourceType);
  itkGetObjectMacro(KernelSource, KernelImageSourceType);

  /** Set/get a single parameter value. */
  virtual void SetParameter(unsigned int index, ParametersValueType value);
  virtual ParametersValueType GetParameter(unsigned int index) const;

  /** Expects the parameters argument to contain values for ALL parameters. */
  virtual void SetParameters(const ParametersType& parameters);

  /** Gets the full parameters list. */
  virtual ParametersType GetParameters() const;

  /** Gets the total number of parameters. */
  virtual unsigned int GetNumberOfParameters() const;

  /** Get/set the z-coordinate of the image z-plane at the given index. */
  void SetZCoordinate(unsigned int index, double coordinate);
  double GetZCoordinate(unsigned int);

  /** Get/set use of custom z coordinates. */
  void SetUseCustomZCoordinates(bool use);
  bool GetUseCustomZCoordinates();

protected:
  BeadSpreadFunctionImageSource();
  virtual ~BeadSpreadFunctionImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** This class is implicitly multi-threaded because its member filters
   * are mulithreaded, so we go with a "single-threaded"
   * implementation here. */
  virtual void GenerateData();

  virtual void GenerateOutputInformation();

private:
  BeadSpreadFunctionImageSource(const BeadSpreadFunctionImageSource&); // purposely not implemented
  void operator=(const BeadSpreadFunctionImageSource&); // purposely not implemented

  double m_IntensityShift; // Additive background constant
  double m_IntensityScale; // The maximum intensity value

  KernelImageSourcePointer  m_KernelSource;
  ConvolverPointer          m_Convolver;
  RescaleImageFilterPointer m_RescaleFilter;
};

} // end namespace itk


#endif // _itkBeadSpreadFunctionImageSource_h
