/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniPSFImageSource.h,v $
  Language:  C++
  Date:      $Date: 2009/09/14 13:57:30 $
  Version:   $Revision: 1.10 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGibsonLanniPSFImageSource_h
#define __itkGibsonLanniPSFImageSource_h

#include <complex>

#include "itkParameterizedImageSource.h"
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
template <class TOutputImage>
class ITK_EXPORT GibsonLanniPSFImageSource : 
  public ParameterizedImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GibsonLanniPSFImageSource         Self;
  typedef ParameterizedImageSource<TOutputImage> Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Typedef for the output image PixelType. */
  typedef TOutputImage                        OutputImageType;
  typedef typename OutputImageType::PixelType PixelType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  itkStaticConstMacro(ImageDimension,
		      unsigned int,
		      TOutputImage::ImageDimension);

  /** Typedef for complex type. */
  typedef std::complex<PixelType> complex_t;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GibsonLanniPSFImageSource,ParameterizedImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::ParametersType      ParametersType;
  
  /** Specify the size of the output image. */
  itkSetVectorMacro(Size,unsigned long,TOutputImage::ImageDimension);

  /** Get the size of the output image. */
  itkGetVectorMacro(Size,unsigned long,TOutputImage::ImageDimension);
  
  /** Specify the spacing of the output image (in nanometers). */
  itkSetVectorMacro(Spacing,float,TOutputImage::ImageDimension);

  /** Get the spacing of the output image (in nanometers). */
  itkGetVectorMacro(Spacing,float,TOutputImage::ImageDimension);

  /** Specify the origin of the output image (in nanometers). */
  itkSetVectorMacro(Origin,float,TOutputImage::ImageDimension);

  /** Get the origin of the output image (in nanometers). */
  itkGetVectorMacro(Origin,float,TOutputImage::ImageDimension);

  /** Specify the point source center (in nanometers). */
  itkSetVectorMacro(PointCenter,float,TOutputImage::ImageDimension);

  /** Get the point source center (in nanometers). */
  itkGetVectorMacro(PointCenter,float,TOutputImage::ImageDimension);

  /** Specify the CCD border width (between the outer edge of the 
      sensing portion of the CCD and the outer edge of the 
      CCD element itself). */
  itkSetVectorMacro(CCDBorderWidth,float,2);

  /** Get the CCD border width. */
  itkGetVectorMacro(CCDBorderWidth,float,2);
  
  /** Specify the emission wavelength (in nanometers). */
  itkSetMacro(EmissionWavelength,float);

  /** Get the emission wavelength (in nanometers). */
  itkGetConstMacro(EmissionWavelength,float);

  /** Specify the numerical aperture (unitless). */
  itkSetMacro(NumericalAperture,float);

  /** Get the numerical aperture (unitless). */
  itkGetConstMacro(NumericalAperture,float);

  /** Specify the magnification (unitless). */
  itkSetMacro(Magnification,float);

  /** Get the magnification (unitless). */
  itkGetConstMacro(Magnification,float);

  /** Specify the design cover slip refractive index (unitless). */
  itkSetMacro(DesignCoverSlipRefractiveIndex,float);

  /** Get the design cover slip refractive index (unitless). */
  itkGetConstMacro(DesignCoverSlipRefractiveIndex,float);

  /** Specify the actual cover slip refractive index (unitless). */
  itkSetMacro(ActualCoverSlipRefractiveIndex,float);

  /** Get the actual cover slip refractive index (unitless). */
  itkGetConstMacro(ActualCoverSlipRefractiveIndex,float);

  /** Specify the design cover slip thickness (in micrometers). */
  itkSetMacro(DesignCoverSlipThickness,float);

  /** Get the design cover slip thickness (in micrometers). */
  itkGetConstMacro(DesignCoverSlipThickness,float);

  /** Specify the actual cover slip thickness (in micrometers). */
  itkSetMacro(ActualCoverSlipThickness,float);

  /** Get the actual cover slip thickness (in micrometers). */
  itkGetConstMacro(ActualCoverSlipThickness,float);

  /** Specify the design immersion oil refractive index (unitless). */
  itkSetMacro(DesignImmersionOilRefractiveIndex,float);

  /** Get the design immersion oil refractive index (unitless). */
  itkGetConstMacro(DesignImmersionOilRefractiveIndex,float);

  /** Specify the actual immersion oil refractive index (unitless). */
  itkSetMacro(ActualImmersionOilRefractiveIndex,float);

  /** Get the actual immersion oil refractive index (unitless). */
  itkGetConstMacro(ActualImmersionOilRefractiveIndex,float);

  /** Specify the design immersion oil thickness (in micrometers). */
  itkSetMacro(DesignImmersionOilThickness,float);

  /** Get the actual immersion oil refractive index (in micrometers). */
  itkGetConstMacro(DesignImmersionOilThickness,float);

  /** Specify the design specimen layer refractive index (unitless). */
  itkSetMacro(DesignSpecimenLayerRefractiveIndex,float);

  /** Get the design specimen layer refractive index (unitless). */
  itkGetConstMacro(DesignSpecimenLayerRefractiveIndex,float);

  /** Specify the actual specimen layer refractive index (unitless). */
  itkSetMacro(ActualSpecimenLayerRefractiveIndex,float);

  /** Get the actual specimen layer refractive index (unitless). */
  itkGetConstMacro(ActualSpecimenLayerRefractiveIndex,float);

  /** Specify the actual point source depth in the specimen layer (in nanometers). */
  itkSetMacro(ActualPointSourceDepthInSpecimenLayer,float);

  /** Get the actual point source depth in the specimen layer (in nanometers). */
  itkGetConstMacro(ActualPointSourceDepthInSpecimenLayer,float);

  /** Specify the design distance from the back focal plane to the detector (in millimeters). */
  itkSetMacro(DesignDistanceFromBackFocalPlaneToDetector,float);

  /** Get the design distance from the back focal plane to the detector (in millimeters). */
  itkGetConstMacro(DesignDistanceFromBackFocalPlaneToDetector,float);

  /** Specify the actual distance from the back focal plane to the detector (in millimeters). */
  itkSetMacro(ActualDistanceFromBackFocalPlaneToDetector,float);

  /** Get the actual distance from the back focal plane to the detector (in millimeters). */
  itkGetConstMacro(ActualDistanceFromBackFocalPlaneToDetector,float);

  /** Expects the parameters argument to contain values for ALL parameters. */
  virtual void SetParameters(const ParametersType& parameters);

  /** Gets the full parameters list. */
  virtual ParametersType GetParameters() const;

  /** Gets the total number of parameters. */
  virtual unsigned int GetNumberOfParameters() const;

  static float BesselFunctionZeroOrderFirstKind(float x);

protected:
  GibsonLanniPSFImageSource();
  ~GibsonLanniPSFImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  virtual void 
  ThreadedGenerateData(const OutputImageRegionType& 
                       outputRegionForThread, int threadId );
  virtual void GenerateOutputInformation();

  complex_t OPD_term(float NA, float n_oil, float rho, float n, float t);

  complex_t OPD(float rho, float delta_z, float a);
  
  void PrecomputeOPDTerms(complex_t* opdCache, float z_o);

  complex_t IntegralTerm(complex_t* opdCache, float K, float a, float z_d,
			 int rhoIndex, float h, float r_o, float z_o);

  /** Computes the light intensity at a specified point. */
  float ComputeSampleValue(complex_t* opdCache,
			   typename TOutputImage::PointType& point);

  /** Computes the integrated light intensity over a CCD pixel centered at
      point. */
  float ComputeIntegratedPixelValue(complex_t* opdCache,
				    typename TOutputImage::PointType& point);

private:
  GibsonLanniPSFImageSource(const GibsonLanniPSFImageSource&); //purposely not implemented
  void operator=(const GibsonLanniPSFImageSource&); //purposely not implemented

  unsigned long *m_Size;         //size of the output image
  float         *m_Spacing;      //spacing
  float         *m_Origin;       //origin
  float         *m_PointCenter;  // the center of the point source

  mutable float m_CCDBorderWidth[2]; // size of border around CCD (x,y)

  /** Point-spread function model parameters. */
  float m_EmissionWavelength;
  float m_NumericalAperture;
  float m_Magnification;
  float m_DesignCoverSlipRefractiveIndex;
  float m_ActualCoverSlipRefractiveIndex;
  float m_DesignCoverSlipThickness;
  float m_ActualCoverSlipThickness;
  float m_DesignImmersionOilRefractiveIndex;
  float m_ActualImmersionOilRefractiveIndex;
  float m_DesignImmersionOilThickness;
  float m_DesignSpecimenLayerRefractiveIndex;
  float m_ActualSpecimenLayerRefractiveIndex;
  float m_ActualPointSourceDepthInSpecimenLayer;
  float m_DesignDistanceFromBackFocalPlaneToDetector;
  float m_ActualDistanceFromBackFocalPlaneToDetector;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGibsonLanniPSFImageSource.cxx"
#endif

#endif
