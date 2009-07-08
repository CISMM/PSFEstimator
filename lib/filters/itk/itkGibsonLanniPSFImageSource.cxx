/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniPSFImageSource.cxx,v $
  Language:  C++
  Date:      $Date: 2009/07/08 18:55:48 $
  Version:   $Revision: 1.2 $

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
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"

#ifndef M_PI
#define M_PI 3.14159265f
#endif // M_PI

#define INTEGRATE_M 50
#define INTEGRATE_N (2*INTEGRATE_M+1)
 
namespace itk
{

/**
 *
 */
template <class TOutputImage>
GibsonLanniPSFImageSource<TOutputImage>
::GibsonLanniPSFImageSource()
{
  m_Size    = new unsigned long [TOutputImage::GetImageDimension()];
  m_Spacing = new float [TOutputImage::GetImageDimension()];
  m_Origin  = new float [TOutputImage::GetImageDimension()];  

  //Initial image is 63 wide in each direction.
  for (unsigned int i=0; i<TOutputImage::GetImageDimension(); i++)
    {
    m_Size[i] = 63;
    m_Spacing[i] = 65.0;
    m_Origin[i] = 0.0;
    }
  m_Size[1] = 1;

  // Set default PSF model parameters.
  m_EmissionWavelength   = 550.0f; // in nanometers
  m_NumericalAperture    = 1.4f;   // unitless
  m_Magnification        = 60.0f;  // unitless
  m_MechanicalTubeLength = 160.0f; // in millimeters

  m_DesignCoverSlipRefractiveIndex    = 1.522f; // unitless
  m_ActualCoverSlipRefractiveIndex    = 1.522f; // unitless
  m_DesignCoverSlipThickness          = 170.0f; // in micrometers
  m_ActualCoverSlipThickness          = 170.0f; // in micrometers
  m_DesignImmersionOilRefractiveIndex = 1.515f; // unitless 
  m_ActualImmersionOilRefractiveIndex = 1.515f; // unitless
  m_DesignImmersionOilThickness       = 100.0f; // in micrometers

  m_DesignSpecimenLayerRefractiveIndex         =  1.33f; // unitless
  m_ActualSpecimenLayerRefractiveIndex         =  1.33f; // unitless
  m_ActualPointSourceDepthInSpecimenLayer      =   0.0f; // in micrometers
  m_DesignDistanceFromBackFocalPlaneToDetector = 170.0f; // in millimeters
  m_ActualDistanceFromBackFocalPlaneToDetector = 170.0f; // in millimeters

}


template <class TOutputImage>
GibsonLanniPSFImageSource<TOutputImage>
::~GibsonLanniPSFImageSource()
{
  delete [] m_Size;
  delete [] m_Spacing;
  delete [] m_Origin;
}


template <class TOutputImage>
float
GibsonLanniPSFImageSource<TOutputImage>
::BesselFunctionZeroOrderFirstKind(float x) {
  float ax, z;
  float xx,y,ans,ans1,ans2;

  if ((ax=fabs(x)) < 8.0) {
    float p1 = 57568490574.0, p2 = -13362590354.0,
      p3 = 651619640.7, p4 = -11214424.18,
      p5 = 77392.33017, p6 = -184.9052456,
      p7 = 57568490411.0, p8 = 1029532985.0,
      p9 = 9494680.718, p10 = 59272.64853,
      p11 = 267.8532712;

    y = x*x;
    ans1=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*p6))));
    ans2=p7+y*(p8+y*(p9+y*(p10+y*(p11+y))));
    ans=ans1/ans2;

  } else {
    float p1 = 0.785398164, p2 = -0.1098628627e-2,
      p3 = 0.2734510407e-4, p4 = -0.2073370639e-5,
      p5 = 0.2093887211e-6, p6 = -0.1562499995e-1,
      p7 = 0.1430488765e-3, p8 = -0.6911147651e-5,
      p9 = 0.7621095161e-6, p10 = 0.934945152e-7,
      p11 = 0.636619772;

    z=8.0/ax;
    y=z*z;
    xx=ax-p1;
    ans1=1.0+y*(p2+y*(p3+
		      y*(p4+y*p5)));
    ans2=p6+y*(p7+
	       y*(p8+y*(p9 -
			y*p10)));
    ans=sqrt(p11/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}


/**
 *
 */
template <class TOutputImage>
void 
GibsonLanniPSFImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  unsigned int i;
  os << indent << "Origin: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Origin[i] << ", ";
    }
  os << m_Origin[i] << "]" << std::endl;

  os << indent << "Spacing: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Spacing[i] << ", ";
    }
  os << m_Spacing[i] << "] (nanometers)" << std::endl;

  os << indent << "Size: [";
  for (i=0; i < TOutputImage::ImageDimension - 1; i++)
    {
    os << m_Size[i] << ", ";
    }
  os << m_Size[i] << "]" << std::endl;
  
  os << indent << "EmissionWavelength (nanometers): "
     << m_EmissionWavelength << std::endl;
  
  os << indent << "NumericalAperture: "
     << m_NumericalAperture << std::endl;

  os << indent << "Magnification: "
     << m_Magnification << std::endl;

  os << indent << "MechanicalTubeLength (millimeters): "
     << m_MechanicalTubeLength << std::endl;

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
  TOutputImage *output;
  typename TOutputImage::IndexType index = {{0}};
  typename TOutputImage::SizeType size = {{0}};
  size.SetSize( m_Size );
  
  output = this->GetOutput(0);

  typename TOutputImage::RegionType largestPossibleRegion;
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
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId )
{
  itkDebugMacro(<<"Generating a random image of scalars");

  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
       
  typedef typename TOutputImage::PixelType ScalarType;
  typename TOutputImage::Pointer image=this->GetOutput(0);

  ImageRegionIteratorWithIndex<TOutputImage> it(image, outputRegionForThread);

  int zSlice = -1;
  complex_t opdCache[INTEGRATE_N];

  for (; !it.IsAtEnd(); ++it)
    {
    typename TOutputImage::IndexType index = it.GetIndex();
    typename TOutputImage::PointType point;
    image->TransformIndexToPhysicalPoint(index, point);

    // See if we have switched slices. If so, we need to precompute some
    // integral terms for the new slice.
    if (zSlice != index[2]) {
      PrecomputeOPDTerms(opdCache, point[2] * 1e-9);
      zSlice = index[2];
    }

    it.Set( ComputePixelValue(opdCache, point) );
    progress.CompletedPixel();
    }
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
typename GibsonLanniPSFImageSource<TOutputImage>::complex_t
GibsonLanniPSFImageSource<TOutputImage>
::OPD_term(float NA, float n_oil, float rho, float n, float t) {
  float NA_rho_sq = NA*NA*rho*rho;
  float NA_rho_over_n_sq = NA_rho_sq/(n*n);
  float n_oil_over_n_sq = (n_oil*n_oil) / (n*n);
  float NA_rho_over_n_oil_sq = NA_rho_sq/(n_oil*n_oil);

  complex_t sq1(1.0f - NA_rho_over_n_sq);
  sq1 = sqrt(sq1);

  complex_t sq2(1.0f - NA_rho_over_n_oil_sq);
  sq2 = n_oil_over_n_sq * sqrt(sq2);

  complex_t result = n*t*(sq1 - sq2);
  return result;
}

//----------------------------------------------------------------------------
template <typename TOutputImage>
typename GibsonLanniPSFImageSource<TOutputImage>::complex_t
GibsonLanniPSFImageSource<TOutputImage>
::OPD(float rho, float delta_z, float a) {
  float NA      = m_NumericalAperture;
  float n_oil_d = m_DesignImmersionOilRefractiveIndex;
  float n_oil   = m_ActualImmersionOilRefractiveIndex;
  float t_oil_d = m_DesignImmersionOilThickness;
  float z_d_d   = m_DesignDistanceFromBackFocalPlaneToDetector * 1e-3;
  float z_d     = m_ActualDistanceFromBackFocalPlaneToDetector * 1e-3;
  float n_s     = m_ActualSpecimenLayerRefractiveIndex;
  float t_s     = m_ActualPointSourceDepthInSpecimenLayer * 1e-9;
  float n_g_d   = m_DesignCoverSlipRefractiveIndex;
  float n_g     = m_ActualCoverSlipRefractiveIndex;
  float t_g_d   = m_DesignCoverSlipThickness * 1e-6;
  float t_g     = m_ActualCoverSlipThickness * 1e-6;

  float r = n_oil *
    (delta_z + (((z_d_d - z_d)*a*a*n_oil) / (z_d_d*z_d*NA*NA)));

  complex_t c1(1.0f - ((NA*NA*rho*rho)/(n_oil*n_oil)));
  c1 = sqrt(c1);

  complex_t c2((a*a*rho*rho*(z_d_d - z_d)) /
               (2.0f*n_oil*z_d_d*z_d));

  complex_t t1 = OPD_term(NA, n_oil, rho, n_s, t_s);
  complex_t t2 = OPD_term(NA, n_oil, rho, n_g, t_g);
  complex_t t3 = OPD_term(NA, n_oil, rho, n_g_d, t_g_d);
  complex_t t4 = OPD_term(NA, n_oil, rho, n_oil_d, t_oil_d);

  complex_t result = r*c1 + c2 + t1 + t2 - t3 - t4;
  return result;
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
void
GibsonLanniPSFImageSource<TOutputImage>
::PrecomputeOPDTerms(complex_t* opdCache, float z_o) {

  float K = 2.0f*M_PI / (m_EmissionWavelength * 1e-9);
  float h = 1.0f / static_cast<float>(INTEGRATE_N-1);
  float NA = m_NumericalAperture;
  float mag = m_Magnification;
  float z_d_d = m_DesignDistanceFromBackFocalPlaneToDetector * 1e-3;
  float a = (z_d_d*NA) / sqrt(mag*mag - NA*NA);

  for (int i = 0; i < INTEGRATE_N; i++) {
    float rho = static_cast<float>(i)*h;
    complex_t W = OPD(rho, z_o, a) * K;
    complex_t I(0.0f, 1.0f);
    opdCache[i] = exp(I*W);
  }

}


//----------------------------------------------------------------------------
template <typename TOutputImage>
typename GibsonLanniPSFImageSource<TOutputImage>::complex_t
GibsonLanniPSFImageSource<TOutputImage>
::IntegralTerm(complex_t* opdCache, float K, float a, float z_d,
	       int rhoIndex, float h, float r_o, float z_o) {
  float rho = static_cast<float>(rhoIndex)*h;
  float bessel = BesselFunctionZeroOrderFirstKind(K*a*rho*r_o/z_d);

  return bessel*opdCache[rhoIndex]*rho;
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
float
GibsonLanniPSFImageSource<TOutputImage>
::ComputePixelValue(complex_t* opdCache, typename TOutputImage::PointType& point) {
  typedef typename TOutputImage::PixelType ScalarType;

  ScalarType px = point[0] * 1e-9;
  ScalarType py = point[1] * 1e-9;
  ScalarType pz = point[2] * 1e-9;

  /* Compute terms that are independent of terms within the integral. */
  float K = 2.0f*M_PI / (m_EmissionWavelength * 1e-9);
  float NA = m_NumericalAperture;
  float mag = m_Magnification;
  float z_d = m_ActualDistanceFromBackFocalPlaneToDetector * 1e-3;
  float a = (z_d*NA) / sqrt(mag*mag - NA*NA);

  // We have to convert to coordinates of the detector points
  float x_o = px * mag;
  float y_o = py * mag;
  float z_o = pz; // No conversion needed

  // Compute common terms in all steps of the integration
  float r_o = sqrt((x_o*x_o) + (y_o*y_o));

  // Compute integration of the formula
  float h = 1.0f / static_cast<float>(INTEGRATE_N-1);

  // Accumulator for integration.
  complex<ScalarType> sum(0.0f, 0.0f);

  // Compute initial terms in Simpson quadrature method.
  sum += IntegralTerm(opdCache, K, a, z_d, 0, h, r_o, z_o);

  sum += IntegralTerm(opdCache, K, a, z_d, INTEGRATE_N-1, h, r_o, z_o);

  for (int k = 1; k <= INTEGRATE_M-1; k++) {
    sum += 2.0f*IntegralTerm(opdCache, K, a, z_d, 2*k, h, r_o, z_o);
  }

  for (int k = 1; k <= INTEGRATE_M; k++) {
    sum += 4.0f*IntegralTerm(opdCache, K, a, z_d, 2*k-1, h, r_o, z_o);
  }

  sum *= (h/3.0f);

  // Return squared magnitude of the integrated value
  return (ScalarType) norm(sum);
}


} // end namespace itk

#endif
