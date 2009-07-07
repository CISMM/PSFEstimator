#ifndef _DATA_MODEL_H_
#define _DATA_MODEL_H_

#include <string>

#include <itkCastImageFilter.h>
#include <itkGibsonLanniPSFImageSource.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkShiftScaleImageFilter.h>

#include <vtkAlgorithmOutput.h>

#include "ITKImageToVTKImage.h"

// This is the data model for the Fibrin Analysis library.
class DataModel {

  typedef float FloatPixelType;
  static const unsigned int Dimension3 = 3;
  typedef itk::Image<FloatPixelType, Dimension3> 
    Float3DImageType;
  typedef itk::ImageFileReader<Float3DImageType>
    Float3DImageReaderType;
  typedef itk::VTKImageExport<Float3DImageType>
    Float3DExporterType;

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

public:

  DataModel();
  virtual ~DataModel();

  void LoadImageFile(std::string fileName);
  void SavePSFImageFile(std::string fileName);

  std::string GetImageFileName();

  void SetNumberOfThreads(int threads);
  int  GetNumberOfThreads();

  void SetImageData(TImage::Pointer image);
  TImage::Pointer GetImageData();

  // Returns the VTK output port for the original scalar image data.
  vtkAlgorithmOutput* GetImageOutputPort();

  double GetImageDataMinimum();
  double GetImageDataMaximum();

  void GetDimensions(int dimensions[3]);

  void SetVoxelSpacing(double spacing[3]);
  void SetVoxelSpacing(int dimension, double spacing);
  void SetVoxelXSpacing(double spacing);
  void SetVoxelYSpacing(double spacing);
  void SetVoxelZSpacing(double spacing);
  void GetVoxelSpacing(double spacing[3]);

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

  void  SetGLMechanicalTubeLength(float mtl) {
    m_GibsonLanniPSFSource->SetMechanicalTubeLength(mtl); }
  float GetGLMechanicalTubeLength() {
    return m_GibsonLanniPSFSource->GetMechanicalTubeLength(); }

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

  void  SetGLDesignPointSourceDistanceFromBackFocalPlane(float distance) {
    m_GibsonLanniPSFSource->SetDesignPointSourceDistanceFromBackFocalPlane(distance); }
  float GetGLDesignPointSourceDistanceFromBackFocalPlane() {
    return m_GibsonLanniPSFSource->GetDesignPointSourceDistanceFromBackFocalPlane(); }
  void  SetGLActualPointSourceDistanceFromBackFocalPlane(float distance) {
    m_GibsonLanniPSFSource->SetActualPointSourceDistanceFromBackFocalPlane(distance); }
  float GetGLActualPointSourceDistanceFromBackFocalPlane() {
    return m_GibsonLanniPSFSource->GetActualPointSourceDistanceFromBackFocalPlane(); }

protected:
  std::string m_ImageFileName;

  TImage::Pointer m_ImageData;

  GibsonLanniPSFImageSourcePointer m_GibsonLanniPSFSource;

  MinMaxType::Pointer m_MinMaxFilter;

  ITKImageToVTKImage<TImage>* m_InputImageITKToVTKFilter;

  ITKImageToVTKImage<TImage>* m_PSFImageITKToVTKFilter;

};

// _DATA_MODEL_H_
#endif
