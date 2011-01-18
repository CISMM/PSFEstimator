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
#ifndef __itkOPDBasedWidefieldMicroscopePointSpreadFunctionFunctor_txx
#define __itkOPDBasedWidefieldMicroscopePointSpreadFunctionFunctor_txx

#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand.h"

namespace itk
{
namespace Functor
{


OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand
::OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand()
{
}


OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand
::~OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand()
{
}


template< class TSource >
void
OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand
::CopySettings(const TSource* source)
{
  this->m_EmissionWavelength                    = source->GetEmissionWavelength();
  this->m_NumericalAperture                     = source->GetNumericalAperture();
  this->m_Magnification                         = source->GetMagnification();
  this->m_DesignCoverSlipRefractiveIndex        = source->GetDesignCoverSlipRefractiveIndex();
  this->m_ActualCoverSlipRefractiveIndex        = source->GetActualCoverSlipRefractiveIndex();
  this->m_DesignCoverSlipThickness              = source->GetDesignCoverSlipThickness();
  this->m_ActualCoverSlipThickness              = source->GetActualCoverSlipThickness();
  this->m_DesignImmersionOilRefractiveIndex     = source->GetDesignImmersionOilRefractiveIndex();
  this->m_ActualImmersionOilRefractiveIndex     = source->GetActualImmersionOilRefractiveIndex();
  this->m_DesignImmersionOilThickness           = source->GetDesignImmersionOilThickness();
  this->m_DesignSpecimenLayerRefractiveIndex    = source->GetDesignSpecimenLayerRefractiveIndex();
  this->m_ActualSpecimenLayerRefractiveIndex    = source->GetActualSpecimenLayerRefractiveIndex();
  this->m_ActualPointSourceDepthInSpecimenLayer = source->GetActualPointSourceDepthInSpecimenLayer();

  this->m_K = 2.0 * itk::Math::pi / (this->m_EmissionWavelength * 1e-9);

  double NA = this->m_NumericalAperture;
  double M  = this->m_Magnification;

  // Assumes classical 160mm tube length.
  this->m_A = 0.160 * NA / sqrt(M*M - NA*NA);
}


inline
OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand::ComplexType
OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand
::OPD(double rho, double dz) const
{
  double NA      = this->m_NumericalAperture;
  double n_oil_d = this->m_DesignImmersionOilRefractiveIndex;
  double n_oil   = this->m_ActualImmersionOilRefractiveIndex;
  double t_oil_d = this->m_DesignImmersionOilThickness * 1e-6;
  double n_s     = this->m_ActualSpecimenLayerRefractiveIndex;
  double t_s     = this->m_ActualPointSourceDepthInSpecimenLayer * 1e-6;
  double n_g_d   = this->m_DesignCoverSlipRefractiveIndex;
  double n_g     = this->m_ActualCoverSlipRefractiveIndex;
  double t_g_d   = this->m_DesignCoverSlipThickness * 1e-6;
  double t_g     = this->m_ActualCoverSlipThickness * 1e-6;

  ComplexType t1 = this->OPDTerm(rho, n_s,     t_s);
  ComplexType t2 = this->OPDTerm(rho, n_g,     t_g);
  ComplexType t3 = this->OPDTerm(rho, n_g_d,   t_g_d);
  ComplexType t4 = this->OPDTerm(rho, n_oil_d, t_oil_d);

  ComplexType c1 = n_oil * dz * sqrt(1.0 - ((NA*NA*rho*rho)/(n_oil*n_oil)));
  ComplexType result = c1 + t1 + t2 - t3 - t4;

  return result;
}


inline
OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand::ComplexType
OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand
::OPDTerm(double rho, double n, double t) const
{
  double NA    = this->m_NumericalAperture;
  double n_oil = this->m_ActualImmersionOilRefractiveIndex;
  double NA_rho_sq = NA*NA*rho*rho;
  double NA_rho_over_n_sq = NA_rho_sq/(n*n);
  double n_oil_over_n_sq = (n_oil*n_oil) / (n*n);
  double NA_rho_over_n_oil_sq = NA_rho_sq/(n_oil*n_oil);

  ComplexType sq1(1.0 - NA_rho_over_n_sq);
  sq1 = sqrt(sq1);

  ComplexType sq2(1.0 - NA_rho_over_n_oil_sq);
  sq2 = n_oil_over_n_sq * sqrt(sq2);

  ComplexType result = n*t*(sq1 - sq2);
  return result;
}


} // end namespace Functor

} // end namespace itk

#endif // __itkOPDBasedWidefieldMicroscopePointSpreadFunctionFunctor_txx
