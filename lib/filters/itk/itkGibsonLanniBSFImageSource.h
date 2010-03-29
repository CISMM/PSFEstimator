/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniBSFImageSource.h,v $
  Language:  C++
  Date:      $Date: 2010/03/29 05:36:32 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGibsonLanniBSFImageSource_h
#define __itkGibsonLanniBSFImageSource_h

#include "itkGibsonLanniPSFImageSource.h"
#include "itkSphereConvolutionFilter.h"
#include "itkParameterizedImageSource.h"
#include "itkNumericTraits.h"


#define DelegateSetMacro(name, type) \
virtual void Set##name(type data) \
{ \
  m_PSFSource->Set##name(data); \
  this->Modified(); \
}

#define DelegateGetMacro(name, type) \
virtual type Get##name() const \
{ \
  return m_PSFSource->Get##name(); \
}


namespace itk
{

/** \class GibsonLanniBSFImageSource
 * \brief Generate a synthetic point-spread function according to the
 * Gibson-Lanni model and convolve it with a spherical bead shape.
 *
 * The Gibson-Lanni point-spread function model takes into account optical
 * path differences caused by imaing conditions differing from the design 
 * conditions of an objective in a widefield fluorescence microscope. This 
 * image source generates images according to this model. IMPORTANT: Please 
 * pay attention to the units each method expects. Some take nanometers, some
 * take micrometers, and some take millimeters.
 *
 * Typically, measuring the point-spread function of an image is done by
 * taking an image stack of a sub-diffraction limit bead (100-150nm). This
 * results in a point-spread function that is convolved with the spherical
 * bead shape. This class calculates the theoretical bead-spread function.
 * It is essentially a wrapper around the GibsonLanniPSFImageSource and
 * SphereConvolutionFilter.
 *
 * \ingroup DataSources Multithreaded
 * \sa GibsonLanniPSFImageSource.
 */
template <class TOutputImage>
class ITK_EXPORT GibsonLanniBSFImageSource : 
  public ParameterizedImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef GibsonLanniBSFImageSource              Self;
  typedef ParameterizedImageSource<TOutputImage> Superclass;
  typedef SmartPointer<Self>                     Pointer;
  typedef SmartPointer<const Self>               ConstPointer;

  /** Typedef for output types. */
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::PixelType  PixelType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  itkStaticConstMacro(ImageDimension,
		      unsigned int,
		      TOutputImage::ImageDimension);

  typedef GibsonLanniPSFImageSource<TOutputImage> PSFSourceType;
  typedef typename PSFSourceType::Pointer         PSFSourcePointer;
  typedef SphereConvolutionFilter<TOutputImage,TOutputImage> ConvolverType;
  typedef typename ConvolverType::Pointer                    ConvolverPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GibsonLanniBSFImageSource,ParameterizedImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef typename Superclass::ParametersValueType ParametersValueType;
  typedef typename Superclass::ParametersType      ParametersType;
  
  /** Specify the size of the output image. */
  void SetSize(unsigned long size[TOutputImage::ImageDimension]) {
    for (int i = 0; i < TOutputImage::ImageDimension; i++) {
      m_Size[i] = size[i];
    }

    m_Convolver->SetSize(size);
    this->Modified();
  }

  /** Get the size of the output image. */
  itkGetVectorMacro(Size,unsigned long,TOutputImage::ImageDimension);
  
  /** Specify the spacing of the output image (in nanometers). */
  void SetSpacing(float spacing[TOutputImage::ImageDimension]) {
    for (int i = 0; i < TOutputImage::ImageDimension; i++) {
      m_Spacing[i] = spacing[i];
    }
    m_Convolver->SetSpacing(spacing);
    this->Modified();
  }

  /** Get the spacing of the output image (in nanometers). */
  itkGetVectorMacro(Spacing,float,TOutputImage::ImageDimension);

  /** Specify the origin of the output image (in nanometers). */
  void SetOrigin(float origin[TOutputImage::ImageDimension]) {
    for (int i = 0; i < TOutputImage::ImageDimension; i++) {
      m_Origin[i] = origin[i];
    }

    m_Convolver->SetOrigin(origin);
    this->Modified();
  }

  /** Get the origin of the output image (in nanometers). */
  itkGetVectorMacro(Origin,float,TOutputImage::ImageDimension);

  /** Specify the point source center (in nanometers). */
  void SetBeadCenter(float center[TOutputImage::ImageDimension]) {
    m_Convolver->SetSphereCenter(center);
    this->Modified();
  }

  /** Get the point source center (in nanometers). */
  float * GetBeadCenter() const {
    return m_Convolver->GetSphereCenter();
  }

  /** Specify the bead radius (in nanometers). */
  void SetBeadRadius(float radius) {
    m_Convolver->SetSphereRadius(radius);
    this->Modified();
  }

  /** Get the bead radius. */
  float GetBeadRadius() const {
    return m_Convolver->GetSphereRadius();
  }

  /** Specify the shear in the X direction. */
  void SetShearX(float shear) {
    m_Convolver->SetShearX(shear);
  }

  /** Get the shear in the X direction. */
  float GetShearX() const {
    return m_Convolver->GetShearX();
  }

  /** Specify the shear in the Y direction. */
  void SetShearY(float shear) {
    m_Convolver->SetShearY(shear);
  }

  /** Get the shear in the Y direction. */
  float GetShearY() const {
    return m_Convolver->GetShearY();
  }

  /** Specify the emission wavelength (in nanometers). */
  DelegateSetMacro(EmissionWavelength,float);

  /** Get the emission wavelength (in nanometers). */
  DelegateGetMacro(EmissionWavelength,float);

  /** Specify the numerical aperture (unitless). */
  DelegateSetMacro(NumericalAperture,float);

  /** Get the numerical aperture (unitless). */
  DelegateGetMacro(NumericalAperture,float);

  /** Specify the magnification (unitless). */
  DelegateSetMacro(Magnification,float);

  /** Get the magnification (unitless). */
  DelegateGetMacro(Magnification,float);

  /** Specify the design cover slip refractive index (unitless). */
  DelegateSetMacro(DesignCoverSlipRefractiveIndex,float);

  /** Get the design cover slip refractive index (unitless). */
  DelegateGetMacro(DesignCoverSlipRefractiveIndex,float);

  /** Specify the actual cover slip refractive index (unitless). */
  DelegateSetMacro(ActualCoverSlipRefractiveIndex,float);

  /** Get the actual cover slip refractive index (unitless). */
  DelegateGetMacro(ActualCoverSlipRefractiveIndex,float);

  /** Specify the design cover slip thickness (in micrometers). */
  DelegateSetMacro(DesignCoverSlipThickness,float);

  /** Get the design cover slip thickness (in micrometers). */
  DelegateGetMacro(DesignCoverSlipThickness,float);

  /** Specify the actual cover slip thickness (in micrometers). */
  DelegateSetMacro(ActualCoverSlipThickness,float);

  /** Get the actual cover slip thickness (in micrometers). */
  DelegateGetMacro(ActualCoverSlipThickness,float);

  /** Specify the design immersion oil refractive index (unitless). */
  DelegateSetMacro(DesignImmersionOilRefractiveIndex,float);

  /** Get the design immersion oil refractive index (unitless). */
  DelegateGetMacro(DesignImmersionOilRefractiveIndex,float);

  /** Specify the actual immersion oil refractive index (unitless). */
  DelegateSetMacro(ActualImmersionOilRefractiveIndex,float);

  /** Get the actual immersion oil refractive index (unitless). */
  DelegateGetMacro(ActualImmersionOilRefractiveIndex,float);

  /** Specify the design immersion oil thickness (in micrometers). */
  DelegateSetMacro(DesignImmersionOilThickness,float);

  /** Get the actual immersion oil refractive index (in micrometers). */
  DelegateGetMacro(DesignImmersionOilThickness,float);

  /** Specify the design specimen layer refractive index (unitless). */
  DelegateSetMacro(DesignSpecimenLayerRefractiveIndex,float);

  /** Get the design specimen layer refractive index (unitless). */
  DelegateGetMacro(DesignSpecimenLayerRefractiveIndex,float);

  /** Specify the actual specimen layer refractive index (unitless). */
  DelegateSetMacro(ActualSpecimenLayerRefractiveIndex,float);

  /** Get the actual specimen layer refractive index (unitless). */
  DelegateGetMacro(ActualSpecimenLayerRefractiveIndex,float);

  /** Specify the actual point source depth in the specimen layer (in nanometers). */
  DelegateSetMacro(ActualPointSourceDepthInSpecimenLayer,float);

  /** Get the actual point source depth in the specimen layer (in nanometers). */
  DelegateGetMacro(ActualPointSourceDepthInSpecimenLayer,float);

  /** Specify the design distance from the back focal plane to the detector (in millimeters). */
  DelegateSetMacro(DesignDistanceFromBackFocalPlaneToDetector,float);

  /** Get the design distance from the back focal plane to the detector (in millimeters). */
  DelegateGetMacro(DesignDistanceFromBackFocalPlaneToDetector,float);

  /** Specify the actual distance from the back focal plane to the detector (in millimeters). */
  DelegateSetMacro(ActualDistanceFromBackFocalPlaneToDetector,float);

  /** Get the actual distance from the back focal plane to the detector (in millimeters). */
  DelegateGetMacro(ActualDistanceFromBackFocalPlaneToDetector,float);

  /** Expects the parameters argument to contain values for ALL parameters. */
  virtual void SetParameters(const ParametersType& parameters);

  /** Gets the full parameters list. */
  virtual ParametersType GetParameters() const;

  /** Gets the total number of parameters. */
  virtual unsigned int GetNumberOfParameters() const;

  /** Get/set the z-coordinate of the image z-plane at the given index. */
  void SetZCoordinate(unsigned int index, double coordinate);
  double GetZCoordinate(unsigned int);

protected:
  GibsonLanniBSFImageSource();
  ~GibsonLanniBSFImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** This class is implicitly multi-threaded because its member filters
   * are mulithreaded, so we go with a "single-threaded" implementation here. */
  virtual void GenerateData();

  virtual void GenerateOutputInformation();


private:
  GibsonLanniBSFImageSource(const GibsonLanniBSFImageSource&); //purposely not implemented
  void operator=(const GibsonLanniBSFImageSource&); //purposely not implemented

  unsigned long *m_Size;        // Size of the output image
  float         *m_Spacing;     // Spacing
  float         *m_Origin;      // Origin
  float         *m_BeadCenter;  // The center of the bead
  float          m_BeadRadius;  // The radius of the bead

  PSFSourcePointer m_PSFSource;
  ConvolverPointer  m_Convolver;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGibsonLanniBSFImageSource.cxx"
#endif

#endif
