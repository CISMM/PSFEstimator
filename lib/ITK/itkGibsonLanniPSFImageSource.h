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
#ifndef __itkGibsonLanniPSFImageSource_h
#define __itkGibsonLanniPSFImageSource_h

#include <complex>

#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class GibsonLanniPSFImageSource
 * \brief Generate a synthetic point-spread function according to the
 * Gibson-Lanni model.
 *
 * The Gibson-Lanni point-spread function model takes into account optical
 * path differences from the design conditions of an objective in a
 * widefield fluorescence microscope. This image source generates images
 * according to this model. IMPORTANT: Please pay attention to the units
 * each method expects. Some take nanometers, some take micrometers, and some
 * take millimeters.
 *
 * \ingroup DataSources Multithreaded
 */
template< class TOutputImage >
class ITK_EXPORT GibsonLanniPSFImageSource :
  public OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef GibsonLanniPSFImageSource                              Self;
  typedef OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage > Superclass;
  typedef SmartPointer< Self >                                   Pointer;
  typedef SmartPointer< const Self >                             ConstPointer;

  /** Typedef for the output image PixelType. */
  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::PixelType      PixelType;
  typedef typename OutputImageType::IndexType      IndexType;
  typedef typename OutputImageType::RegionType     RegionType;
  typedef typename OutputImageType::PointType      PointType;
  typedef typename OutputImageType::PointValueType PointValueType;
  typedef typename OutputImageType::SpacingType    SpacingType;
  typedef typename OutputImageType::SizeType       SizeType;
  typedef typename OutputImageType::SizeValueType  SizeValueType;


  itkStaticConstMacro(ImageDimension, unsigned int,
		      TOutputImage::ImageDimension);

  /** Typedef for complex type. */
  typedef std::complex<double> ComplexType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GibsonLanniPSFImageSource, ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::ParametersType      ParametersType;

  /** Set a single parameter value. */
  virtual void SetParameter(unsigned int index, ParametersValueType value);

  /** Get a single parameter value. */
  virtual ParametersValueType GetParameter(unsigned int index) const;

  /** Expects the parameters argument to contain values for ALL parameters. */
  virtual void SetParameters(const ParametersType& parameters);

  /** Gets the full parameters list. */
  virtual ParametersType GetParameters() const;

  /** Gets the total number of parameters. */
  virtual unsigned int GetNumberOfParameters() const;

protected:
  GibsonLanniPSFImageSource();
  ~GibsonLanniPSFImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void
  ThreadedGenerateData(const RegionType& outputRegionForThread, int threadId );
  virtual void GenerateOutputInformation();

  ComplexType OPD_term(double NA, double n_oil, double rho, double n, double t);

  ComplexType OPD(double rho, double delta_z, double a);

  void PrecomputeOPDTerms(ComplexType* opdCache, double z_o);

  inline ComplexType IntegralTerm(ComplexType* opdCache, double K, double a, double z_d,
                                  int rhoIndex, double h, double r_o, double z_o);

  /** Computes the light intensity at a specified point. */
  double ComputeSampleValue(ComplexType* opdCache, PointType& point);

  /** Computes the integrated light intensity over a CCD pixel centered at
      point. */
  double ComputeIntegratedPixelValue(ComplexType* opdCache,
				    PointType& point);

private:
  GibsonLanniPSFImageSource(const GibsonLanniPSFImageSource&); //purposely not implemented
  void operator=(const GibsonLanniPSFImageSource&); //purposely not implemented

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGibsonLanniPSFImageSource.txx"
#endif

#endif
