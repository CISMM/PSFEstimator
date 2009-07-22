#ifndef _DATA_MODEL_H_
#define _DATA_MODEL_H_

#include <string>

#include <itkGibsonLanniPSFImageSource.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkShiftScaleImageFilter.h>

#include <itkAmoebaOptimizer.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkImageToParameterizedImageSourceMetric.h>

#include <vtkAlgorithmOutput.h>

#include "ITKImageToVTKImage.h"

// This is the data model for the Fibrin Analysis library.
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
  typedef itk::ImageToParameterizedImageSourceMetric<TImage, GibsonLanniPSFImageSourceType>
    ParameterizedCostFunctionType;
  typedef itk::NormalizedCorrelationImageToImageMetric<TImage, TImage>
    ImageToImageCostFunctionType;
  typedef itk::AmoebaOptimizer OptimizerType;

  DataModel();
  virtual ~DataModel();

  void LoadImageFile(std::string fileName);
  void SavePSFImageFile(std::string fileName);

  std::string GetMeasuredImageFileName();

  // Number of threads to run for all the multithreaded ITK algorithms
  void SetNumberOfThreads(int threads);
  int  GetNumberOfThreads();

  void SetMeasuredImageData(TImage::Pointer image);
  TImage::Pointer GetMeasuredImageData();

  TImage::Pointer GetPSFImageData();

  // Returns the VTK output port for the original scalar image data.
  vtkAlgorithmOutput* GetMeasuredImageOutputPort();

  // Returns the VTK output port for the generated PSF image data.
  vtkAlgorithmOutput* GetPSFImageOutputPort();

  double GetMeasuredImageDataMinimum();
  double GetMeasuredImageDataMaximum();
  Float3DPointType GetMeasuredImageDataMaximumCoordinates();

  void GetMeasuredImageDimensions(int dimensions[3]);

  void SetMeasuredImageVoxelSpacing(double spacing[3]);
  void SetMeasuredImageVoxelSpacing(int dimension, double spacing);
  void SetMeasuredImageVoxelXSpacing(double spacing);
  void SetMeasuredImageVoxelYSpacing(double spacing);
  void SetMeasuredImageVoxelZSpacing(double spacing);
  void GetMeasuredImageVoxelSpacing(double spacing[3]);

  void SetMeasuredImageOrigin(double origin[3]);
  void GetMeasuredImageOrigin(double origin[3]);

  double GetPSFImageDataMinimum();
  double GetPSFImageDataMaximum();

  void SetPSFImageDimensions(int dimensions[3]);
  void SetPSFImageDimension(int index, int dimension);
  void SetPSFImageXDimension(int dimension);
  void SetPSFImageYDimension(int dimension);
  void SetPSFImageZDimension(int dimension);
  void GetPSFImageDimensions(int dimensions[3]);

  void SetPSFImageVoxelSpacing(double spacing[3]);
  void SetPSFImageVoxelSpacing(int dimension, double spacing);
  void SetPSFImageVoxelXSpacing(double spacing);
  void SetPSFImageVoxelYSpacing(double spacing);
  void SetPSFImageVoxelZSpacing(double spacing);
  void GetPSFImageVoxelSpacing(double spacing[3]);

  void SetPSFImageOrigin(double origin[3]);
  void GetPSFImageOrigin(double origin[3]);

  // Recenters the PSF image origin to the center of the image bounds.
  void RecenterPSFImageOrigin();

  // Sets the PSF center
  void SetPSFPointCenter(double center[3]);
  void GetPSFPointCenter(double center[3]);

  void UpdateGibsonLanniPSFImage();

  void  SetGLEmissionWavelength(float wavelength) {
    m_GibsonLanniPSFSource->SetEmissionWavelength(wavelength); }    
  float GetGLEmissionWavelength() { 
    return m_GibsonLanniPSFSource->GetEmissionWavelength(); }

  void  SetGLNumericalAperture(float na) {
    m_GibsonLanniPSFSource->SetNumericalAperture(na); }
  float GetGLNumericalAperture() {
    return m_GibsonLanniPSFSource->GetNumericalAperture(); }

  void  SetGLMagnification(float magnification) {
    m_GibsonLanniPSFSource->SetMagnification(magnification); }
  float GetGLMagnification() {
    return m_GibsonLanniPSFSource->GetMagnification(); }

  void  SetGLDesignCoverSlipRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetDesignCoverSlipRefractiveIndex(ri); }
  float GetGLDesignCoverSlipRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetDesignCoverSlipRefractiveIndex(); }
  void  SetGLActualCoverSlipRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetActualCoverSlipRefractiveIndex(ri); }
  float GetGLActualCoverSlipRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetActualCoverSlipRefractiveIndex(); }

  void  SetGLDesignCoverSlipThickness(float thickness) {
    m_GibsonLanniPSFSource->SetDesignCoverSlipThickness(thickness); }
  float GetGLDesignCoverSlipThickness() {
    return m_GibsonLanniPSFSource->GetDesignCoverSlipThickness(); }
  void  SetGLActualCoverSlipThickness(float thickness) {
    m_GibsonLanniPSFSource->SetActualCoverSlipThickness(thickness); }
  float GetGLActualCoverSlipThickness() {
    return m_GibsonLanniPSFSource->GetActualCoverSlipThickness(); }

  void  SetGLDesignImmersionOilRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetDesignImmersionOilRefractiveIndex(ri); }
  float GetGLDesignImmersionOilRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetDesignImmersionOilRefractiveIndex(); }
  void  SetGLActualImmersionOilRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetActualImmersionOilRefractiveIndex(ri); }
  float GetGLActualImmersionOilRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetActualImmersionOilRefractiveIndex(); }

  void  SetGLDesignImmersionOilThickness(float thickness) {
    m_GibsonLanniPSFSource->SetDesignImmersionOilThickness(thickness); }
  float GetGLDesignImmersionOilThickness() {
    return m_GibsonLanniPSFSource->GetDesignImmersionOilThickness(); }

  void  SetGLDesignSpecimenLayerRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetDesignSpecimenLayerRefractiveIndex(ri); }
  float GetGLDesignSpecimenLayerRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetDesignSpecimenLayerRefractiveIndex(); }
  void  SetGLActualSpecimenLayerRefractiveIndex(float ri) {
    m_GibsonLanniPSFSource->SetActualSpecimenLayerRefractiveIndex(ri); }
  float GetGLActualSpecimenLayerRefractiveIndex() {
    return m_GibsonLanniPSFSource->GetActualSpecimenLayerRefractiveIndex(); }

  void  SetGLActualPointSourceDepthInSpecimenLayer(float depth) {
    m_GibsonLanniPSFSource->SetActualPointSourceDepthInSpecimenLayer(depth); }
  float GetGLActualPointSourceDepthInSpecimenLayer() {
    return m_GibsonLanniPSFSource->GetActualPointSourceDepthInSpecimenLayer(); }

  void  SetGLDesignDistanceFromBackFocalPlaneToDetector(float distance) {
    m_GibsonLanniPSFSource->SetDesignDistanceFromBackFocalPlaneToDetector(distance); }
  float GetGLDesignDistanceFromBackFocalPlaneToDetector() {
    return m_GibsonLanniPSFSource->GetDesignDistanceFromBackFocalPlaneToDetector(); }
  void  SetGLActualDistanceFromBackFocalPlaneToDetector(float distance) {
    m_GibsonLanniPSFSource->SetActualDistanceFromBackFocalPlaneToDetector(distance); }
  float GetGLActualDistanceFromBackFocalPlaneToDetector() {
    return m_GibsonLanniPSFSource->GetActualDistanceFromBackFocalPlaneToDetector(); }

  void SetGLParameterEnabled(unsigned int index, bool enabled);
  bool GetGLParameterEnabled(unsigned int index);

  double GetImageComparisonMetric();

  void Optimize();

protected:
  std::string m_ImageFileName;

  TImage::Pointer m_MeasuredImageData;

  GibsonLanniPSFImageSourcePointer m_GibsonLanniPSFSource;

  MinMaxType::Pointer m_MeasuredImageMinMaxFilter;
  MinMaxType::Pointer m_PSFImageMinMaxFilter;

  ITKImageToVTKImage<TImage>* m_MeasuredImageITKToVTKFilter;
  ITKImageToVTKImage<TImage>* m_PSFImageITKToVTKFilter;

  // The cost function used by the optimizer.
  ParameterizedCostFunctionType::Pointer m_CostFunction;

  // The delegate cost function used by m_CostFunction
  ImageToImageCostFunctionType::Pointer  m_ImageToImageCostFunction;
  
  // The optimizer
  OptimizerType::Pointer m_Optimizer;
};

// _DATA_MODEL_H_
#endif
