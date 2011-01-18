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
#ifndef __itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand_h
#define __itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand_h

#include <complex>

namespace itk
{
namespace Functor
{

/** \class OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand
 * \brief Functor class that defines a method for computing the
 * optical path difference (OPD) in several widefield point-spread
 * function models.
 */

class OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand {
public:
  typedef std::complex<double> ComplexType;

  OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand();
  virtual ~OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand();

  template< class TSource >
  void CopySettings(const TSource* source);

  /** Point-spread function model parameters. */
  double    m_EmissionWavelength;
  double    m_NumericalAperture;
  double    m_Magnification;
  double    m_DesignCoverSlipRefractiveIndex;
  double    m_ActualCoverSlipRefractiveIndex;
  double    m_DesignCoverSlipThickness;
  double    m_ActualCoverSlipThickness;
  double    m_DesignImmersionOilRefractiveIndex;
  double    m_ActualImmersionOilRefractiveIndex;
  double    m_DesignImmersionOilThickness;
  double    m_DesignSpecimenLayerRefractiveIndex;
  double    m_ActualSpecimenLayerRefractiveIndex;
  double    m_ActualPointSourceDepthInSpecimenLayer;

  /** Precomputed values. */
  double m_K; // Wavenumber
  double m_A; // Radius of projection of the limiting aperture onto
              // the back focal plane of the objective lens

protected:
  /** Computes the optical path difference for a ray terminating at
  *   a normalized distance rho from the center of the back focal
  *   plane aperture. */
  ComplexType OPD(double rho, double dz) const;

  /** Common terms for computing the optical path difference term in
   *  point-spread function models descended from this class. */
  ComplexType OPDTerm(double rho, double n, double t) const;

};

} // end namespace Functor

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand.txx"
#endif

#endif // __itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand_h
