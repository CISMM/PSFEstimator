/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniPSFImageSource.h,v $
  Language:  C++
  Date:      $Date: 2010/04/19 18:50:02 $
  Version:   $Revision: 1.12 $

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

#include "itkParametricImageSource.h"
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
  public ParametricImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GibsonLanniPSFImageSource           Self;
  typedef ParametricImageSource<TOutputImage> Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

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

  /** Specify the size of the output image. */
  virtual void SetSize(const SizeType & size)
  {
    if (size != m_Size)
      {
      m_Size = size;
      this->Modified();
      }
  }

  /** Get the size of the output image. */
  itkGetConstReferenceMacro(Size, SizeType);

  /** Specify the spacing of the output image (in nanometers). */
  virtual void SetSpacing(const SpacingType & spacing)
  {
    if (spacing != m_Spacing)
      {
      m_Spacing = spacing;
      this->Modified();
      }
  }

  /** Get the spacing of the output image (in nanometers). */
  itkGetConstReferenceMacro(Spacing, SpacingType);

  /** Specify the origin of the output image (in nanometers). */
  virtual void SetOrigin(const PointType & origin)
  {
    if (origin != m_Origin)
      {
      m_Origin = origin;
      this->Modified();
      }
  }

  /** Get the origin of the output image (in nanometers). */
  itkGetConstReferenceMacro(Origin, PointType);

  /** Specify the point source center (in nanometers). */
  virtual void SetPointCenter(const PointType & center)
  {
    if (center != m_PointCenter)
      {
      m_PointCenter = center;
      this->Modified();
      }
  }

  /** Get the point source center (in nanometers). */
  itkGetConstReferenceMacro(PointCenter, PointType);

  /** Specify the X shear. */
  itkSetMacro(ShearX, double);

  /** Get the X shear. */
  itkGetConstMacro(ShearX, double);

  /** Specify the Y shear. */
  itkSetMacro(ShearY, double);

  /** Get the Y shear. */
  itkGetConstMacro(ShearY, double);

  /** Specify the emission wavelength (in nanometers). */
  itkSetMacro(EmissionWavelength, double);

  /** Get the emission wavelength (in nanometers). */
  itkGetConstMacro(EmissionWavelength, double);

  /** Specify the numerical aperture (unitless). */
  itkSetMacro(NumericalAperture, double);

  /** Get the numerical aperture (unitless). */
  itkGetConstMacro(NumericalAperture, double);

  /** Specify the magnification (unitless). */
  itkSetMacro(Magnification, double);

  /** Get the magnification (unitless). */
  itkGetConstMacro(Magnification, double);

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

  /** Specify the design distance from the back focal plane to the detector (in millimeters). */
  itkSetMacro(DesignDistanceFromBackFocalPlaneToDetector, double);

  /** Get the design distance from the back focal plane to the detector (in millimeters). */
  itkGetConstMacro(DesignDistanceFromBackFocalPlaneToDetector, double);

  /** Specify the actual distance from the back focal plane to the detector (in millimeters). */
  itkSetMacro(ActualDistanceFromBackFocalPlaneToDetector, double);

  /** Get the actual distance from the back focal plane to the detector (in millimeters). */
  itkGetConstMacro(ActualDistanceFromBackFocalPlaneToDetector, double);

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

  SizeType    m_Size;        // size of the output image
  SpacingType m_Spacing;     // spacing
  PointType   m_Origin;      // origin
  PointType   m_PointCenter; // the center of the point source
  double      m_ShearX;      // Shear in the x-direction with respect to z
  double      m_ShearY;      // Shear in the y-direction with respect to z

  /** Point-spread function model parameters. */
  double m_EmissionWavelength;
  double m_NumericalAperture;
  double m_Magnification;
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
  double m_DesignDistanceFromBackFocalPlaneToDetector;
  double m_ActualDistanceFromBackFocalPlaneToDetector;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGibsonLanniPSFImageSource.txx"
#endif

#endif
