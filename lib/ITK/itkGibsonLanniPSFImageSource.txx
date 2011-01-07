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
  m_Size.Fill(32);
  m_Spacing.Fill(65.0);
  m_Origin.Fill(0.0);
  m_PointCenter.Fill(0.0);

  m_ShearX = 0.0;
  m_ShearY = 0.0;

  // Set default PSF model parameters.
  m_EmissionWavelength   = 550.0; // in nanometers
  m_NumericalAperture    = 1.4;   // unitless
  m_Magnification        = 60.0;  // unitless

  m_DesignCoverSlipRefractiveIndex    = 1.522; // unitless
  m_ActualCoverSlipRefractiveIndex    = 1.522; // unitless
  m_DesignCoverSlipThickness          = 170.0; // in micrometers
  m_ActualCoverSlipThickness          = 170.0; // in micrometers
  m_DesignImmersionOilRefractiveIndex = 1.515; // unitless
  m_ActualImmersionOilRefractiveIndex = 1.515; // unitless
  m_DesignImmersionOilThickness       = 100.0; // in micrometers

  m_DesignSpecimenLayerRefractiveIndex         =  1.33; // unitless
  m_ActualSpecimenLayerRefractiveIndex         =  1.33; // unitless
  m_ActualPointSourceDepthInSpecimenLayer      =   0.0; // in micrometers
  m_DesignDistanceFromBackFocalPlaneToDetector = 160.0; // in millimeters
  m_ActualDistanceFromBackFocalPlaneToDetector = 160.0; // in millimeters
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
::SetParameters(const ParametersType& parameters)
{
  int index = 0;
  SpacingType spacing;
  spacing[0] = parameters[index++];
  spacing[1] = parameters[index++];
  spacing[2] = parameters[index++];
  SetSpacing(spacing);

  PointType center;
  center[0] = parameters[index++];
  center[1] = parameters[index++];
  center[2] = parameters[index++];
  SetPointCenter(center);

  SetShearX(parameters[index++]);
  SetShearY(parameters[index++]);

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
  for ( unsigned int i = 0; i < m_Spacing.Size(); i++ )
    {
    parameters[index++] = GetSpacing()[i];
    }

  for ( unsigned int i = 0; i < m_PointCenter.Size(); i++ )
    {
    parameters[index++] = GetPointCenter()[i];
    }

  parameters[index++] = GetShearX();
  parameters[index++] = GetShearY();

  parameters[index++] = GetEmissionWavelength();
  parameters[index++] = GetNumericalAperture();
  parameters[index++] = GetMagnification();

  parameters[index++] = GetDesignCoverSlipRefractiveIndex();
  parameters[index++] = GetActualCoverSlipRefractiveIndex();
  parameters[index++] = GetDesignCoverSlipThickness();
  parameters[index++] = GetActualCoverSlipThickness();
  parameters[index++] = GetDesignImmersionOilRefractiveIndex();
  parameters[index++] = GetActualImmersionOilRefractiveIndex();
  parameters[index++] = GetDesignImmersionOilThickness();

  parameters[index++] = GetDesignSpecimenLayerRefractiveIndex();
  parameters[index++] = GetActualSpecimenLayerRefractiveIndex();
  parameters[index++] = GetActualPointSourceDepthInSpecimenLayer();
  parameters[index++] = GetDesignDistanceFromBackFocalPlaneToDetector();
  parameters[index++] = GetActualDistanceFromBackFocalPlaneToDetector();

  return parameters;
}


//----------------------------------------------------------------------------
template <class TOutputImage>
unsigned int
GibsonLanniPSFImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return 23;
}


//----------------------------------------------------------------------------
template <class TOutputImage>
void
GibsonLanniPSFImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  unsigned int i;
  os << indent << "Origin: [";
  for ( i=0; i < m_Origin.Size() - 1; i++ )
    {
    os << m_Origin[i] << ", ";
    }
  os << m_Origin[i] << "]" << std::endl;

  os << indent << "Spacing: [";
  for ( i=0; i < m_Spacing.Size() - 1; i++ )
    {
    os << m_Spacing[i] << ", ";
    }
  os << m_Spacing[i] << "] (nanometers)" << std::endl;

  os << indent << "Size: [";
  for ( i=0; i < m_Size.GetSizeDimension() - 1; i++ )
    {
    os << m_Size[i] << ", ";
    }
  os << m_Size[i] << "]" << std::endl;

  os << indent << "PointCenter: [";
  for ( i=0; i < m_PointCenter.Size() - 1; i++ )
    {
    os << m_PointCenter[i] << ", ";
    }
  os << m_PointCenter[i] << "]" << std::endl;

  os << "ShearX: " << m_ShearX << std::endl;
  os << "ShearY: " << m_ShearY << std::endl;

  os << indent << "EmissionWavelength (nanometers): "
     << m_EmissionWavelength << std::endl;

  os << indent << "NumericalAperture: "
     << m_NumericalAperture << std::endl;

  os << indent << "Magnification: "
     << m_Magnification << std::endl;

  os << indent << "DesignCoverSlipRefractiveIndex: "
     << m_DesignCoverSlipRefractiveIndex << std::endl;

  os << indent << "ActualCoverSlipRefractiveIndex: "
     << m_ActualCoverSlipRefractiveIndex << std::endl;

  os << indent << "DesignCoverSlipThickness (micrometers): "
     << m_DesignCoverSlipThickness << std::endl;

  os << indent << "ActualCoverSlipThickness (micrometers): "
     << m_ActualCoverSlipThickness << std::endl;

  os << indent << "DesignImmersionOilRefractiveIndex: "
     << m_DesignImmersionOilRefractiveIndex << std::endl;

  os << indent << "ActualImmersionOilRefractiveIndex: "
     << m_ActualImmersionOilRefractiveIndex << std::endl;

  os << indent << "DesignImmersionOilThickness (micrometers): "
     << m_DesignImmersionOilThickness << std::endl;

  os << indent << "DesignSpecimenLayerRefractiveIndex: "
     << m_DesignSpecimenLayerRefractiveIndex << std::endl;

  os << indent << "ActualSpecimenLayerRefractiveIndex: "
     << m_ActualSpecimenLayerRefractiveIndex << std::endl;

  os << indent << "ActualPointSourceDepthInSpecimenLayer (nanometers): "
     << m_ActualPointSourceDepthInSpecimenLayer << std::endl;

  os << indent << "DesignDistanceFromBackFocalPlaneToDetector (millimeters): "
     << m_DesignDistanceFromBackFocalPlaneToDetector << std::endl;

  os << indent << "ActualDistanceFromBackFocalPlaneToDetector (millimeters): "
     << m_ActualDistanceFromBackFocalPlaneToDetector << std::endl;

}

//----------------------------------------------------------------------------
template <typename TOutputImage>
void
GibsonLanniPSFImageSource<TOutputImage>
::GenerateOutputInformation()
{
  OutputImageType *output;
  IndexType index = {{0}};
  SizeType size( m_Size );

  output = this->GetOutput(0);

  RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( size );
  largestPossibleRegion.SetIndex( index );
  output->SetLargestPossibleRegion( largestPossibleRegion );

  output->SetSpacing(m_Spacing);
  output->SetOrigin(m_Origin);
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
    point[0] = point[0] - m_ShearX*(point[2] - m_PointCenter[2]);
    point[1] = point[1] - m_ShearY*(point[2] - m_PointCenter[2]);

    // Shift the center of the point
    for (int i = 0; i < 3; i++) point[i] -= m_PointCenter[i];

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
  double NA      = m_NumericalAperture;
  double n_oil_d = m_DesignImmersionOilRefractiveIndex;
  double n_oil   = m_ActualImmersionOilRefractiveIndex;
  double t_oil_d = m_DesignImmersionOilThickness * 1e-6;
  double z_d_d   = m_DesignDistanceFromBackFocalPlaneToDetector * 1e-3;
  double z_d     = m_ActualDistanceFromBackFocalPlaneToDetector * 1e-3;
  double n_s     = m_ActualSpecimenLayerRefractiveIndex;
  double t_s     = m_ActualPointSourceDepthInSpecimenLayer * 1e-6;
  double n_g_d   = m_DesignCoverSlipRefractiveIndex;
  double n_g     = m_ActualCoverSlipRefractiveIndex;
  double t_g_d   = m_DesignCoverSlipThickness * 1e-6;
  double t_g     = m_ActualCoverSlipThickness * 1e-6;

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
  double K = 2.0f*itk::Math::pi / (m_EmissionWavelength * 1e-9);
  double h = 1.0f / static_cast<double>(INTEGRATE_N-1);
  double NA = m_NumericalAperture;
  double mag = m_Magnification;
  double z_d_d = m_DesignDistanceFromBackFocalPlaneToDetector * 1e-3;
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
  double K = 2.0f*itk::Math::pi / (m_EmissionWavelength * 1e-9);
  double NA = m_NumericalAperture;
  double mag = m_Magnification;
  double z_d = m_ActualDistanceFromBackFocalPlaneToDetector * 1e-3;
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
  double dx = m_Spacing[0] / static_cast<double>(divs);
  double dy = m_Spacing[1] / static_cast<double>(divs);
  for (int iy = 0; iy < divs; iy++)
    {
    for (int ix = 0; ix < divs; ix++)
      {
      PointType samplePoint;
      double fx = (static_cast<double>(ix) + 0.5f) * dx;
      double fy = (static_cast<double>(iy) + 0.5f) * dy;
      samplePoint[0] = point[0] - 0.5*m_Spacing[0] + fx;
      samplePoint[1] = point[1] - 0.5*m_Spacing[1] + fy;
      samplePoint[2] = point[2];

      integrated += ComputeSampleValue(opdCache, samplePoint);
      }
    }

  return integrated;
}


} // end namespace itk

#endif
