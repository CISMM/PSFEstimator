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
#ifndef __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_h
#define __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_h

#include <complex>

#include "itkWidefieldMicroscopePointSpreadFunctionImageSource.h"

namespace itk
{

/** \class OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource
 * \brief Base class for point-spread function image sources based on
 * mathematical models that incorporate an optical path difference
 * (OPD) term computed from the mismatch between design parameters
 * assumed by an objective lens and actual parameters.
 *
 * \ingroup DataSources Multithreaded
 */
template< class TOutputImage >
class ITK_EXPORT OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource :
  public WidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource         Self;
  typedef WidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage > Superclass;
  typedef SmartPointer< Self >                                              Pointer;
  typedef SmartPointer< const Self >                                        ConstPointer;

  /** Typedef for the output image. */
  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::PixelType      PixelType;
  typedef typename OutputImageType::IndexType      IndexType;
  typedef typename OutputImageType::RegionType     RegionType;
  typedef typename OutputImageType::PointType      PointType;
  typedef typename OutputImageType::PointValueType PointValueType;
  typedef typename OutputImageType::SpacingType    SpacingType;
  typedef typename OutputImageType::SizeType       SizeType;
  typedef typename OutputImageType::SizeValueType  SizeValueType;
  typedef std::complex<double>                     ComplexType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource,
               WidefieldMicroscopePointSpreadFunctionImageSource);

  /** Specify the design cover slip refractive index (unitless). */
  itkSetMacro(DesignCoverSlipRefractiveIndex, double);

  /** Get the design cover slip refractive index (unitless). */
  itkGetConstMacro(DesignCoverSlipRefractiveIndex, double);

  /** Specify the actual cover slip refractive index (unitless). */
  itkSetMacro(ActualCoverSlipRefractiveIndex, double);

  /** Get the actual cover slip refractive index (unitless). */
  itkGetConstMacro(ActualCoverSlipRefractiveIndex, double);

  /** Specify the design cover slip thickness (in micrometers). */
  itkSetMacro(DesignCoverSlipThickness, double);

  /** Get the design cover slip thickness (in micrometers). */
  itkGetConstMacro(DesignCoverSlipThickness, double);

  /** Specify the actual cover slip thickness (in micrometers). */
  itkSetMacro(ActualCoverSlipThickness, double);

  /** Get the actual cover slip thickness (in micrometers). */
  itkGetConstMacro(ActualCoverSlipThickness, double);

  /** Specify the design immersion oil refractive index (unitless). */
  itkSetMacro(DesignImmersionOilRefractiveIndex, double);

  /** Get the design immersion oil refractive index (unitless). */
  itkGetConstMacro(DesignImmersionOilRefractiveIndex, double);

  /** Specify the actual immersion oil refractive index (unitless). */
  itkSetMacro(ActualImmersionOilRefractiveIndex, double);

  /** Get the actual immersion oil refractive index (unitless). */
  itkGetConstMacro(ActualImmersionOilRefractiveIndex, double);

  /** Specify the design immersion oil thickness (in micrometers). */
  itkSetMacro(DesignImmersionOilThickness, double);

  /** Get the design immersion oil refractive index (in micrometers). */
  itkGetConstMacro(DesignImmersionOilThickness, double);

  /** Specify the design specimen layer refractive index (unitless). */
  itkSetMacro(DesignSpecimenLayerRefractiveIndex, double);

  /** Get the design specimen layer refractive index (unitless). */
  itkGetConstMacro(DesignSpecimenLayerRefractiveIndex, double);

  /** Specify the actual specimen layer refractive index (unitless). */
  itkSetMacro(ActualSpecimenLayerRefractiveIndex, double);

  /** Get the actual specimen layer refractive index (unitless). */
  itkGetConstMacro(ActualSpecimenLayerRefractiveIndex, double);

  /** Specify the actual point source depth in the specimen layer (in nanometers). */
  itkSetMacro(ActualPointSourceDepthInSpecimenLayer, double);

  /** Get the actual point source depth in the specimen layer (in nanometers). */
  itkGetConstMacro(ActualPointSourceDepthInSpecimenLayer, double);

protected:
  OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource();
  ~OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Integrates a one-dimensional complex-valued function defined by
   *  a functor from a to b using 2*m+1 subdivisions. */
  template< class TFunctor >
  ComplexType Integrate(const TFunctor& functor, double a, double b,
                        int m, double x, double y, double z)
  {
    int n = 2*m + 1;
    double h = 1.0 / static_cast<double>(n-1);

    double r = sqrt(x*x + y*y);

    // Initialize accumulator for integration.
    ComplexType sum(0.0, 0.0);

    // Compute initial terms in Simpson quadrature method.
    sum += functor(r, z, a);
    sum += functor(r, z, b);

    for (int k = 1; k <= m-1; k++)
      {
      sum += 2.0 * functor(r, z, (2*k)*h);
      }

    for (int k = 1; k <= m; k++)
      {
      sum += 4.0 * functor(r, z, (2*k-1)*h);
      }

    sum *= h / 3.0;

    return sum;
  }

  double m_DesignCoverSlipRefractiveIndex;
  double m_ActualCoverSlipRefractiveIndex;
  double m_DesignCoverSlipThickness;
  double m_ActualCoverSlipThickness;
  double m_DesignImmersionOilRefractiveIndex;
  double m_ActualImmersionOilRefractiveIndex;
  double m_DesignImmersionOilThickness;
  double m_DesignSpecimenLayerRefractiveIndex;
  double m_ActualSpecimenLayerRefractiveIndex;
  double m_ActualPointSourceDepthInSpecimenLayer;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource.txx"
#endif

#endif // __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_h
