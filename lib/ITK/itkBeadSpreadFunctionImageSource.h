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
#include "itkCommand.h"

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

  /** Set/get the size of the output image. */
  void SetSize(const SizeType & size);
  const SizeType & GetSize() const;

  /** Set/get the spacing of the output image (in nanometers). */
  void SetSpacing(const SpacingType & spacing);
  const SpacingType & GetSpacing() const;

  /** Set/get the origin of the output image (in nanometers). */
  virtual void SetOrigin(const PointType & origin);
  const PointType & GetOrigin() const;

  /** Set/get the point source center (in nanometers). */
  virtual void SetBeadCenter(const PointType & center);
  const PointType & GetBeadCenter() const;

  /** Set/get the bead radius (in nanometers). */
  void SetBeadRadius(double radius);
  double GetBeadRadius() const;

  /** Set/get the shear in the X direction. */
  void SetShearX(double shear);
  double GetShearX() const;

  /** Set/get the shear in the Y direction. */
  void SetShearY(double shear);
  double GetShearY() const;

  /** Set/get the background value. */
  itkSetMacro(IntensityShift, double);
  itkGetConstMacro(IntensityShift, double);

  /** Set/get the maximum intensity. */
  itkSetMacro(IntensityScale, double);
  itkGetConstMacro(IntensityScale, double);

  /** Set/get the convolution kernel source. */
  virtual void SetKernelSource( KernelImageSourceType* source );
  itkGetObjectMacro(KernelSource, KernelImageSourceType);

  /** Set/get kernel radial symmetry flag. If this flag is set to
   * true, then only a single slice of the kernel corresponding to a
   * radial profile of the kernel. */
  itkSetMacro(KernelIsRadiallySymmetric, bool);
  itkGetMacro(KernelIsRadiallySymmetric, bool);
  itkBooleanMacro(KernelIsRadiallySymmetric);

  /** Set/get a single parameter value. */
  virtual void SetParameter(unsigned int index, ParametersValueType value);
  virtual ParametersValueType GetParameter(unsigned int index) const;

  /** Expects the parameters argument to contain values for ALL parameters. */
  virtual void SetParameters(const ParametersType& parameters);

  /** Gets the full parameters list. */
  virtual ParametersType GetParameters() const;

  /** Gets the total number of parameters. */
  virtual unsigned int GetNumberOfParameters() const;

  /** Gets the number of bead-spread function parameters. */
  virtual unsigned int GetNumberOfBeadSpreadFunctionParameters() const;

  /** Get/set the z-coordinate of the image z-plane at the given index. */
  void SetZCoordinate(unsigned int index, double coordinate);
  double GetZCoordinate(unsigned int);

  /** Get/set use of custom z coordinates. */
  void SetUseCustomZCoordinates(bool use);
  bool GetUseCustomZCoordinates();

  /** Callback evoked whenever the KernelSource is modified. */
  virtual void KernelModified();

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
  bool                      m_KernelIsRadiallySymmetric;
  ConvolverPointer          m_Convolver;
  RescaleImageFilterPointer m_RescaleFilter;

  typedef SimpleMemberCommand< Self > MemberCommandType;
  typedef typename MemberCommandType::Pointer MemberCommandPointer;
  MemberCommandPointer m_ModifiedEventCommand;
  unsigned long        m_ObserverTag;
};

} // end namespace itk


#endif // _itkBeadSpreadFunctionImageSource_h
