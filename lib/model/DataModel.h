#ifndef _DATA_MODEL_H_
#define _DATA_MODEL_H_

#include <string>

#include "Configuration.h"

#include <itkGibsonLanniBSFImageSource.h>
#include <itkGibsonLanniPSFImageSource.h>
#include <itkScanImageFilter.h>
#include <itkSumProjectionImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkShiftScaleImageFilter.h>

#include <itkAmoebaOptimizer.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkImageToParameterizedImageSourceMetric.h>

#include <vtkAlgorithmOutput.h>

#include "ITKImageToVTKImage.h"

// This is the data model for the PSF Optimizer library.
class DataModel {

public:
  typedef float FloatPixelType;
  static const unsigned int Dimension3 = 3;
  typedef itk::Image<FloatPixelType, Dimension3> 
    Float3DImageType;
  typedef itk::ImageFileReader<Float3DImageType>
    Float3DImageReaderType;
  typedef itk::VTKImageExport<Float3DImageType>
    Float3DExporterType;
  typedef Float3DImageType::PointType
    Float3DPointType;

  typedef itk::GibsonLanniPSFImageSource<Float3DImageType>
    GibsonLanniPSFImageSourceType;
  typedef GibsonLanniPSFImageSourceType::Pointer
    GibsonLanniPSFImageSourcePointer;
  typedef itk::GibsonLanniBSFImageSource<Float3DImageType>
    GibsonLanniBSFImageSourceType;
  typedef GibsonLanniBSFImageSourceType::Pointer
    GibsonLanniBSFImageSourcePointer;
  
  typedef Float3DImageType TImage;
  typedef TImage InputImageType;

  typedef itk::MinimumMaximumImageCalculator<TImage>
    MinMaxType;
  typedef itk::ShiftScaleImageFilter<TImage, TImage>
    ScaleFilterType;

  typedef itk::ImageFileReader<TImage>
    ScalarFileReaderType;
  typedef itk::ImageFileWriter<TImage>
    ScalarFileWriterType;

  // For writing out to TIFF images.
  typedef itk::Image<unsigned short, 3>
    TIFFOutputImageType;
  typedef itk::ShiftScaleImageFilter<TImage,TIFFOutputImageType>
    TIFFScaleType;
  typedef itk::ImageFileWriter<TIFFOutputImageType>
    TIFFWriterType;

  // Types for optimization.
  typedef itk::ImageToParameterizedImageSourceMetric<TImage, GibsonLanniBSFImageSourceType>
    ParameterizedCostFunctionType;
  typedef itk::NormalizedCorrelationImageToImageMetric<TImage, TImage>
    ImageToImageCostFunctionType;
  typedef itk::AmoebaOptimizer OptimizerType;

  DataModel();
  virtual ~DataModel();

  void LoadImageFile(std::string fileName);
  void SavePSFImageFile(std::string fileName);
  void SaveBSFImageFile(std::string fileName);

  void SetConfiguration(Configuration & configuration);
  void GetConfiguration(Configuration & configuration);

  std::string GetMeasuredImageFileName();

  // Number of threads to run for all the multithreaded ITK algorithms
  void SetNumberOfThreads(int threads);
  int  GetNumberOfThreads();

  void SetMeasuredImageData(TImage::Pointer image);
  TImage::Pointer GetMeasuredImageData();

  TImage::Pointer GetPSFImageData();
  TImage::Pointer GetBSFImageData();

  // Returns the VTK output port for the original scalar image data.
  vtkAlgorithmOutput* GetMeasuredImageOutputPort();

  // Returns the VTK output port for the generated point-spread function
  // image data.
  vtkAlgorithmOutput* GetPSFImageOutputPort();

  // Returns the VTK output port for the generated bead-spread function
  // image data.
  vtkAlgorithmOutput* GetBSFImageOutputPort();

  double GetMeasuredImageDataMinimum();
  double GetMeasuredImageDataMaximum();
  Float3DPointType GetMeasuredImageDataMaximumCoordinates();

  void GetMeasuredImageDimensions(int dimensions[3]);

  void SetMeasuredImageVoxelSpacing(double spacing[3]);
  void SetMeasuredImageVoxelSpacing(int dimension, double spacing);
  void GetMeasuredImageVoxelSpacing(double spacing[3]);

  void SetMeasuredImageOrigin(double origin[3]);
  void GetMeasuredImageOrigin(double origin[3]);

  double GetPSFImageDataMinimum();
  double GetPSFImageDataMaximum();

  double GetBSFImageDataMinimum();
  double GetBSFImageDataMaximum();

  void SetPSFImageDimensions(int dimensions[3]);
  void SetPSFImageDimension(int index, int dimension);
  void GetPSFImageDimensions(int dimensions[3]);

  void SetBSFImageDimensions(int dimensions[3]);
  void SetBSFImageDimension(int index, int dimension);
  void GetBSFImageDimensions(int dimensions[3]);

  void SetPSFImageVoxelSpacing(double spacing[3]);
  void SetPSFImageVoxelSpacing(int dimension, double spacing);
  void GetPSFImageVoxelSpacing(double spacing[3]);

  void SetBSFImageVoxelSpacing(double spacing[3]);
  void SetBSFImageVoxelSpacing(int dimension, double spacing);
  void GetBSFImageVoxelSpacing(double spacing[3]);

  void SetCCDBorderWidth(double borderWidth[2]);
  void GetCCDBorderWidth(double borderWidth[2]);

  void SetPSFImageOrigin(double origin[3]);
  void GetPSFImageOrigin(double origin[3]);

  void SetBSFImageOrigin(double origin[3]);
  void GetBSFImageOrigin(double origin[3]);

  // Recenters the PSF image origin to the center of the image bounds.
  void RecenterPSFImageOrigin();

  // Sets the PSF center
  void SetPSFPointCenter(double center[3]);
  void GetPSFPointCenter(double center[3]);
  void SetBSFPointCenter(double center[3]);
  void GetBSFPointCenter(double center[3]);

  void UpdateGibsonLanniPSFImage();

  void SetBeadRadius(double radius);
  double GetBeadRadius();

  void  SetGLEmissionWavelength(float wavelength) {
    m_GibsonLanniPSFSource->SetEmissionWavelength(wavelength);
    m_GibsonLanniBSFSource->SetEmissionWavelength(wavelength);
  }
  float GetGLEmissionWavelength() { 
    return m_GibsonLanniPSFSource->GetEmissionWavelength(); }

  void  SetGLNumericalAperture(float na) {
    m_GibsonLanniPSFSource->SetNumericalAperture(na);
    m_GibsonLanniBSFSource->SetNumericalAperture(na);
  }
  float GetGLNumericalAperture() {
    return m_GibsonLanniPSFSource->GetNumericalAperture(); }

  void  SetGLMagnification(float magnification) {
    m_GibsonLanniPSFSource->SetMagnification(magnification);
    m_GibsonLanniBSFSource->SetMagnification(magnification);
  }
  float GetGLMagnification() {
    return m_GibsonLanniPSFSource->GetMagnification(); }

  void  SetGLDesignCoverSlipRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetDesignCoverSlipRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetDesignCoverSlipRefractiveIndex(ri);
  }
  float GetGLDesignCoverSlipRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetDesignCoverSlipRefractiveIndex(); }
  void  SetGLActualCoverSlipRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetActualCoverSlipRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetActualCoverSlipRefractiveIndex(ri);
  }
  float GetGLActualCoverSlipRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetActualCoverSlipRefractiveIndex(); }

  void  SetGLDesignCoverSlipThickness(float thickness) {
    m_GibsonLanniPSFSource->SetDesignCoverSlipThickness(thickness);
    m_GibsonLanniBSFSource->SetDesignCoverSlipThickness(thickness);
  }
  float GetGLDesignCoverSlipThickness() {
    return m_GibsonLanniPSFSource->GetDesignCoverSlipThickness(); }
  void  SetGLActualCoverSlipThickness(float thickness) {
    m_GibsonLanniPSFSource->SetActualCoverSlipThickness(thickness);
    m_GibsonLanniBSFSource->SetActualCoverSlipThickness(thickness);
  }
  float GetGLActualCoverSlipThickness() {
    return m_GibsonLanniPSFSource->GetActualCoverSlipThickness(); }

  void  SetGLDesignImmersionOilRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetDesignImmersionOilRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetDesignImmersionOilRefractiveIndex(ri);
  }
  float GetGLDesignImmersionOilRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetDesignImmersionOilRefractiveIndex(); }
  void  SetGLActualImmersionOilRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetActualImmersionOilRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetActualImmersionOilRefractiveIndex(ri);
  }
  float GetGLActualImmersionOilRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetActualImmersionOilRefractiveIndex(); }

  void  SetGLDesignImmersionOilThickness(float thickness) {
    m_GibsonLanniPSFSource->SetDesignImmersionOilThickness(thickness);
    m_GibsonLanniBSFSource->SetDesignImmersionOilThickness(thickness);
  }
  float GetGLDesignImmersionOilThickness() {
    return m_GibsonLanniPSFSource->GetDesignImmersionOilThickness(); }

  void  SetGLDesignSpecimenLayerRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetDesignSpecimenLayerRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetDesignSpecimenLayerRefractiveIndex(ri);
  }
  float GetGLDesignSpecimenLayerRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetDesignSpecimenLayerRefractiveIndex(); }
  void  SetGLActualSpecimenLayerRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetActualSpecimenLayerRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetActualSpecimenLayerRefractiveIndex(ri);
  }
  float GetGLActualSpecimenLayerRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetActualSpecimenLayerRefractiveIndex(); }

  void  SetGLActualPointSourceDepthInSpecimenLayer(float depth) {
    m_GibsonLanniPSFSource->SetActualPointSourceDepthInSpecimenLayer(depth);
    m_GibsonLanniBSFSource->SetActualPointSourceDepthInSpecimenLayer(depth);
  }
  float GetGLActualPointSourceDepthInSpecimenLayer() {
    return m_GibsonLanniPSFSource->GetActualPointSourceDepthInSpecimenLayer(); }

  void  SetGLDesignDistanceFromBackFocalPlaneToDetector(float distance) {
    m_GibsonLanniPSFSource->SetDesignDistanceFromBackFocalPlaneToDetector(distance);
    m_GibsonLanniBSFSource->SetDesignDistanceFromBackFocalPlaneToDetector(distance);
  }
  float GetGLDesignDistanceFromBackFocalPlaneToDetector() {
    return m_GibsonLanniPSFSource->GetDesignDistanceFromBackFocalPlaneToDetector(); }
  void  SetGLActualDistanceFromBackFocalPlaneToDetector(float distance) {
    m_GibsonLanniPSFSource->SetActualDistanceFromBackFocalPlaneToDetector(distance);
    m_GibsonLanniBSFSource->SetActualDistanceFromBackFocalPlaneToDetector(distance);
  }
  float GetGLActualDistanceFromBackFocalPlaneToDetector() {
    return m_GibsonLanniPSFSource->GetActualDistanceFromBackFocalPlaneToDetector(); }

  void SetGLParameterEnabled(unsigned int index, bool enabled);
  bool GetGLParameterEnabled(unsigned int index);

  double GetImageComparisonMetric();

  void Optimize();

protected:
  Configuration m_Configuration;

  std::string m_ImageFileName;

  TImage::Pointer m_MeasuredImageData;

  GibsonLanniPSFImageSourcePointer m_GibsonLanniPSFSource;
  GibsonLanniBSFImageSourcePointer m_GibsonLanniBSFSource;

  MinMaxType::Pointer m_MeasuredImageMinMaxFilter;
  MinMaxType::Pointer m_PSFImageMinMaxFilter;
  MinMaxType::Pointer m_BSFImageMinMaxFilter;

  ITKImageToVTKImage<TImage>* m_MeasuredImageITKToVTKFilter;
  ITKImageToVTKImage<TImage>* m_PSFImageITKToVTKFilter;
  ITKImageToVTKImage<TImage>* m_BSFImageITKToVTKFilter;

  // The cost function used by the optimizer.
  ParameterizedCostFunctionType::Pointer m_CostFunction;

  // The delegate cost function used by m_CostFunction
  ImageToImageCostFunctionType::Pointer  m_ImageToImageCostFunction;
  
  // The optimizer
  OptimizerType::Pointer m_Optimizer;
};

// _DATA_MODEL_H_
#endif
