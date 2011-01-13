/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniPSFImageSource.cxx,v $
  Language:  C++
  Date:      $Date: 2010/05/17 15:41:35 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGibsonLanniPSFImageSource_txx
#define __itkGibsonLanniPSFImageSource_txx

#include "itkGibsonLanniPSFImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMath.h"
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"

#define INTEGRATE_M 20
#define INTEGRATE_N (2*INTEGRATE_M+1)

namespace itk
{

//----------------------------------------------------------------------------
template <class TOutputImage>
GibsonLanniPSFImageSource<TOutputImage>
::GibsonLanniPSFImageSource()
{
  this->m_Size.Fill(32);
  this->m_Spacing.Fill(65.0);
  this->m_Origin.Fill(0.0);

  // Set default PSF model parameters.
  this->m_EmissionWavelength   = 550.0; // in nanometers
  this->m_NumericalAperture    = 1.4;   // unitless
  this->m_Magnification        = 60.0;  // unitless

  this->m_DesignCoverSlipRefractiveIndex    = 1.522; // unitless
  this->m_ActualCoverSlipRefractiveIndex    = 1.522; // unitless
  this->m_DesignCoverSlipThickness          = 170.0; // in micrometers
  this->m_ActualCoverSlipThickness          = 170.0; // in micrometers
  this->m_DesignImmersionOilRefractiveIndex = 1.515; // unitless
  this->m_ActualImmersionOilRefractiveIndex = 1.515; // unitless
  this->m_DesignImmersionOilThickness       = 100.0; // in micrometers

  this->m_DesignSpecimenLayerRefractiveIndex         =  1.33; // unitless
  this->m_ActualSpecimenLayerRefractiveIndex         =  1.33; // unitless
  this->m_ActualPointSourceDepthInSpecimenLayer      =   0.0; // in micrometers
  this->m_DesignDistanceFromBackFocalPlaneToDetector = 160.0; // in millimeters
  this->m_ActualDistanceFromBackFocalPlaneToDetector = 160.0; // in millimeters
}


//----------------------------------------------------------------------------
template <class TOutputImage>
GibsonLanniPSFImageSource<TOutputImage>
::~GibsonLanniPSFImageSource()
{
}


//----------------------------------------------------------------------------
template <class TOutputImage>
void
GibsonLanniPSFImageSource<TOutputImage>
::SetParameter(unsigned int index, ParametersValueType value)
{
  switch (index)
    {
    case 0:
      this->SetEmissionWavelength(value);
      break;

    case 1:
      this->SetNumericalAperture(value);
      break;

    case 2:
      this->SetMagnification(value);
      break;

    case 3:
      this->SetDesignCoverSlipRefractiveIndex(value);
      break;

    case 4:
      this->SetActualCoverSlipRefractiveIndex(value);
      break;

    case 5:
      this->SetDesignCoverSlipThickness(value);
      break;

    case 6:
      this->SetActualCoverSlipThickness(value);
      break;

    case 7:
      this->SetDesignImmersionOilRefractiveIndex(value);
      break;

    case 8:
      this->SetActualImmersionOilRefractiveIndex(value);
      break;

    case 9:
      this->SetDesignImmersionOilThickness(value);
      break;

    case 10:
      this->SetDesignSpecimenLayerRefractiveIndex(value);
      break;

    case 11:
      this->SetActualSpecimenLayerRefractiveIndex(value);
      break;

    case 12:
      this->SetActualPointSourceDepthInSpecimenLayer(value);
      break;

    case 13:
      this->SetDesignDistanceFromBackFocalPlaneToDetector(value);
      break;

    case 14:
      this->SetActualDistanceFromBackFocalPlaneToDetector(value);
      break;
    }
}


//----------------------------------------------------------------------------
template <class TOutputImage>
typename GibsonLanniPSFImageSource<TOutputImage>::ParametersValueType
GibsonLanniPSFImageSource<TOutputImage>
::GetParameter(unsigned int index) const
{
  switch (index)
    {
    case 0:
      return this->GetEmissionWavelength();
      break;

    case 1:
      return this->GetNumericalAperture();
      break;

    case 2:
      return this->GetMagnification();
      break;

    case 3:
      return this->GetDesignCoverSlipRefractiveIndex();
      break;

    case 4:
      return this->GetActualCoverSlipRefractiveIndex();
      break;

    case 5:
      return this->GetDesignCoverSlipThickness();
      break;

    case 6:
      return this->GetActualCoverSlipThickness();
      break;

    case 7:
      return this->GetDesignImmersionOilRefractiveIndex();
      break;

    case 8:
      return this->GetActualImmersionOilRefractiveIndex();
      break;

    case 9:
      return this->GetDesignImmersionOilThickness();
      break;

    case 10:
      return this->GetDesignSpecimenLayerRefractiveIndex();
      break;

    case 11:
      return this->GetActualSpecimenLayerRefractiveIndex();
      break;

    case 12:
      return this->GetActualPointSourceDepthInSpecimenLayer();
      break;

    case 13:
      return this->GetDesignDistanceFromBackFocalPlaneToDetector();
      break;

    case 14:
      return this->GetActualDistanceFromBackFocalPlaneToDetector();
      break;

    default:
      return 0.0;
      break;
    }
}


//----------------------------------------------------------------------------
template <class TOutputImage>
void
GibsonLanniPSFImageSource<TOutputImage>
::SetParameters(const ParametersType& parameters)
{
  int index = 0;

  SetEmissionWavelength(parameters[index++]);
  SetNumericalAperture(parameters[index++]);
  SetMagnification(parameters[index++]);

  SetDesignCoverSlipRefractiveIndex(parameters[index++]);
  SetActualCoverSlipRefractiveIndex(parameters[index++]);
  SetDesignCoverSlipThickness(parameters[index++]);
  SetActualCoverSlipThickness(parameters[index++]);
  SetDesignImmersionOilRefractiveIndex(parameters[index++]);
  SetActualImmersionOilRefractiveIndex(parameters[index++]);
  SetDesignImmersionOilThickness(parameters[index++]);

  SetDesignSpecimenLayerRefractiveIndex(parameters[index++]);
  SetActualSpecimenLayerRefractiveIndex(parameters[index++]);
  SetActualPointSourceDepthInSpecimenLayer(parameters[index++]);
  SetDesignDistanceFromBackFocalPlaneToDetector(parameters[index++]);
  SetActualDistanceFromBackFocalPlaneToDetector(parameters[index++]);
}


template <class TOutputImage>
typename GibsonLanniPSFImageSource<TOutputImage>::ParametersType
GibsonLanniPSFImageSource<TOutputImage>
::GetParameters() const
{
  ParametersType parameters(GetNumberOfParameters());

  int index = 0;
  parameters[index++] = this->GetEmissionWavelength();
  parameters[index++] = this->GetNumericalAperture();
  parameters[index++] = this->GetMagnification();

  parameters[index++] = this->GetDesignCoverSlipRefractiveIndex();
  parameters[index++] = this->GetActualCoverSlipRefractiveIndex();
  parameters[index++] = this->GetDesignCoverSlipThickness();
  parameters[index++] = this->GetActualCoverSlipThickness();
  parameters[index++] = this->GetDesignImmersionOilRefractiveIndex();
  parameters[index++] = this->GetActualImmersionOilRefractiveIndex();
  parameters[index++] = this->GetDesignImmersionOilThickness();

  parameters[index++] = this->GetDesignSpecimenLayerRefractiveIndex();
  parameters[index++] = this->GetActualSpecimenLayerRefractiveIndex();
  parameters[index++] = this->GetActualPointSourceDepthInSpecimenLayer();
  parameters[index++] = this->GetDesignDistanceFromBackFocalPlaneToDetector();
  parameters[index++] = this->GetActualDistanceFromBackFocalPlaneToDetector();

  return parameters;
}


//----------------------------------------------------------------------------
template <class TOutputImage>
unsigned int
GibsonLanniPSFImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return 15;
}


//----------------------------------------------------------------------------
template <class TOutputImage>
void
GibsonLanniPSFImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
void
GibsonLanniPSFImageSource<TOutputImage>
::ThreadedGenerateData(const RegionType& outputRegionForThread, int threadId )
{
  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  typename TOutputImage::Pointer image = this->GetOutput(0);

  ImageRegionIteratorWithIndex<OutputImageType> it(image, outputRegionForThread);

  int zSlice = -1;
  ComplexType opdCache[INTEGRATE_N];

  for (; !it.IsAtEnd(); ++it)
    {
    IndexType index = it.GetIndex();
    PointType point;
    image->TransformIndexToPhysicalPoint(index, point);

    // Apply x and y shear here
    //point[0] = point[0] - this->m_ShearX*(point[2] - this->m_PointCenter[2]);
    //point[1] = point[1] - this->m_ShearY*(point[2] - this->m_PointCenter[2]);

    // Shift the center of the point
    //for (int i = 0; i < 3; i++) point[i] -= this->m_PointCenter[i];

    // See if we have switched slices. If so, we need to precompute some
    // integral terms for the new slice.
    if (zSlice != index[2])
      {
      PrecomputeOPDTerms(opdCache, point[2] * 1e-9);
      zSlice = index[2];
      }

    it.Set( ComputeSampleValue(opdCache, point) );
    progress.CompletedPixel();
    }
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
typename GibsonLanniPSFImageSource<TOutputImage>::ComplexType
GibsonLanniPSFImageSource<TOutputImage>
::OPD_term(double NA, double n_oil, double rho, double n, double t)
{
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

//----------------------------------------------------------------------------
template <typename TOutputImage>
typename GibsonLanniPSFImageSource<TOutputImage>::ComplexType
GibsonLanniPSFImageSource<TOutputImage>
::OPD(double rho, double delta_z, double a)
{
  double NA      = this->m_NumericalAperture;
  double n_oil_d = this->m_DesignImmersionOilRefractiveIndex;
  double n_oil   = this->m_ActualImmersionOilRefractiveIndex;
  double t_oil_d = this->m_DesignImmersionOilThickness * 1e-6;
  double z_d_d   = this->m_DesignDistanceFromBackFocalPlaneToDetector * 1e-3;
  double z_d     = this->m_ActualDistanceFromBackFocalPlaneToDetector * 1e-3;
  double n_s     = this->m_ActualSpecimenLayerRefractiveIndex;
  double t_s     = this->m_ActualPointSourceDepthInSpecimenLayer * 1e-6;
  double n_g_d   = this->m_DesignCoverSlipRefractiveIndex;
  double n_g     = this->m_ActualCoverSlipRefractiveIndex;
  double t_g_d   = this->m_DesignCoverSlipThickness * 1e-6;
  double t_g     = this->m_ActualCoverSlipThickness * 1e-6;

  double r = n_oil *
    (delta_z + (((z_d_d - z_d)*a*a*n_oil) / (z_d_d*z_d*NA*NA)));

  ComplexType c1(1.0 - ((NA*NA*rho*rho)/(n_oil*n_oil)));
  c1 = sqrt(c1);

  ComplexType c2((a*a*rho*rho*(z_d_d - z_d)) /
               (2.0f*n_oil*z_d_d*z_d));

  ComplexType t1 = OPD_term(NA, n_oil, rho, n_s, t_s);
  ComplexType t2 = OPD_term(NA, n_oil, rho, n_g, t_g);
  ComplexType t3 = OPD_term(NA, n_oil, rho, n_g_d, t_g_d);
  ComplexType t4 = OPD_term(NA, n_oil, rho, n_oil_d, t_oil_d);

  ComplexType result = r*c1 + c2 + t1 + t2 - t3 - t4;
  return result;
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
void
GibsonLanniPSFImageSource<TOutputImage>
::PrecomputeOPDTerms(ComplexType* opdCache, double z_o)
{
  double K = 2.0f*itk::Math::pi / (this->m_EmissionWavelength * 1e-9);
  double h = 1.0f / static_cast<double>(INTEGRATE_N-1);
  double NA = this->m_NumericalAperture;
  double mag = this->m_Magnification;
  double z_d_d = this->m_DesignDistanceFromBackFocalPlaneToDetector * 1e-3;
  double a = (z_d_d*NA) / sqrt(mag*mag - NA*NA);

  for (int i = 0; i < INTEGRATE_N; i++) {
    double rho = static_cast<double>(i)*h;
    ComplexType W = OPD(rho, z_o, a) * K;
    ComplexType I(0.0f, 1.0f);
    opdCache[i] = exp(I*W);
  }

}


//----------------------------------------------------------------------------
template <typename TOutputImage>
typename GibsonLanniPSFImageSource<TOutputImage>::ComplexType
GibsonLanniPSFImageSource<TOutputImage>
::IntegralTerm(ComplexType* opdCache, double K, double a, double z_d,
	       int rhoIndex, double h, double r_o, double z_o)
{
  double rho = static_cast<double>(rhoIndex)*h;
  double bessel = j0(K*a*rho*r_o/z_d);

  return bessel*opdCache[rhoIndex]*rho;
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
double
GibsonLanniPSFImageSource<TOutputImage>
::ComputeSampleValue(ComplexType* opdCache, typename TOutputImage::PointType& point)
{
  PixelType px = point[0] * 1e-9;
  PixelType py = point[1] * 1e-9;
  PixelType pz = point[2] * 1e-9;

  /* Compute terms that are independent of terms within the integral. */
  double K = 2.0f*itk::Math::pi / (this->m_EmissionWavelength * 1e-9);
  double NA = this->m_NumericalAperture;
  double mag = this->m_Magnification;
  double z_d = this->m_ActualDistanceFromBackFocalPlaneToDetector * 1e-3;
  double a = (z_d*NA) / sqrt(mag*mag - NA*NA);

  // We have to convert to coordinates of the detector points
  double x_o = px * mag;
  double y_o = py * mag;
  double z_o = pz; // No conversion needed

  // Compute common terms in all steps of the integration
  double r_o = sqrt((x_o*x_o) + (y_o*y_o));

  // Compute integration of the formula
  double h = 1.0 / static_cast<double>(INTEGRATE_N-1);

  // Accumulator for integration.
  ComplexType sum(0.0, 0.0);

  // Compute initial terms in Simpson quadrature method.
  sum += IntegralTerm(opdCache, K, a, z_d, 0, h, r_o, z_o);

  sum += IntegralTerm(opdCache, K, a, z_d, INTEGRATE_N-1, h, r_o, z_o);

  for (int k = 1; k <= INTEGRATE_M-1; k++)
    {
    sum += 2.0*IntegralTerm(opdCache, K, a, z_d, 2*k, h, r_o, z_o);
    }

  for (int k = 1; k <= INTEGRATE_M; k++)
    {
    sum += 4.0*IntegralTerm(opdCache, K, a, z_d, 2*k-1, h, r_o, z_o);
    }

  sum *= (h/3.0f);

  // Return squared magnitude of the integrated value
  return (PixelType) norm(sum);
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
double
GibsonLanniPSFImageSource<TOutputImage>
::ComputeIntegratedPixelValue(ComplexType* opdCache, typename TOutputImage::PointType& point)
{
  double integrated = 0.0f;

  // Evaluate over a grid
  int divs = 1;
  double dx = this->m_Spacing[0] / static_cast<double>(divs);
  double dy = this->m_Spacing[1] / static_cast<double>(divs);
  for (int iy = 0; iy < divs; iy++)
    {
    for (int ix = 0; ix < divs; ix++)
      {
      PointType samplePoint;
      double fx = (static_cast<double>(ix) + 0.5f) * dx;
      double fy = (static_cast<double>(iy) + 0.5f) * dy;
      samplePoint[0] = point[0] - 0.5*this->m_Spacing[0] + fx;
      samplePoint[1] = point[1] - 0.5*this->m_Spacing[1] + fy;
      samplePoint[2] = point[2];

      integrated += ComputeSampleValue(opdCache, samplePoint);
      }
    }

  return integrated;
}


} // end namespace itk

#endif
