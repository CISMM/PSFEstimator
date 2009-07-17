#ifndef _DATA_MODEL_CXX_
#define _DATA_MODEL_CXX_

#if defined(_WIN32) // Turn off deprecation warnings in Visual Studio
#pragma warning( disable : 4996 )
#endif

#include <cstdlib>

#include "DataModel.h"

#include <itkMultiThreader.h>
#include <itkPoint.h>

DataModel
::DataModel() {
  m_MeasuredImageData = NULL;
  
  m_GibsonLanniPSFSource = GibsonLanniPSFImageSourceType::New();

  m_MeasuredImageMinMaxFilter = MinMaxType::New();
  m_PSFImageMinMaxFilter      = MinMaxType::New();

  m_MeasuredImageITKToVTKFilter = new ITKImageToVTKImage<TImage>();
  m_PSFImageITKToVTKFilter   = new ITKImageToVTKImage<TImage>();

  m_ImageToImageCostFunction = ImageToImageCostFunctionType::New();
  m_CostFunction = ParameterizedCostFunctionType::New();
  m_CostFunction->SetImageToImageMetric(m_ImageToImageCostFunction);

  // ITK will detect the number of cores on the system and set it by default.
  // Here we can override that setting if the proper environment variable is
  // set.
  char *var = getenv("VisualPSFOptimizer_THREADS");
  if (var) {
    int numberOfThreads = atoi(var);
    if (numberOfThreads > 0)
      SetNumberOfThreads(numberOfThreads);
  }

}


DataModel
::~DataModel() {
  delete m_MeasuredImageITKToVTKFilter;
  delete m_PSFImageITKToVTKFilter;
}


void
DataModel
::LoadImageFile(std::string fileName) {
  m_ImageFileName = fileName;
  ScalarFileReaderType::Pointer reader = ScalarFileReaderType::New();
  reader->SetFileName(fileName.c_str());
  reader->Update();
  SetMeasuredImageData(reader->GetOutput());
  
  // Connect this image data to the various pipelines.
  m_MeasuredImageMinMaxFilter->SetImage(m_MeasuredImageData);
  Float3DImageType::RegionType region = 
    m_MeasuredImageData->GetLargestPossibleRegion();
  m_MeasuredImageMinMaxFilter->SetRegion(region);
  m_MeasuredImageMinMaxFilter->Compute();
  
  m_MeasuredImageITKToVTKFilter->Modified();
  m_MeasuredImageITKToVTKFilter->Update();

  // Set up PSF settings to match loaded image settings
  double spacing[3];
  GetMeasuredImageVoxelSpacing(spacing);
  SetPSFImageVoxelSpacing(spacing);

  int size[3];
  GetMeasuredImageDimensions(size);
  SetPSFImageDimensions(size);
  
  float origin[3];
  for (int i = 0; i < 3; i++)
    origin[i] = -spacing[i]*static_cast<float>(size[i])*0.5;
  m_GibsonLanniPSFSource->SetOrigin(origin);
  m_GibsonLanniPSFSource->Update();
  m_MeasuredImageData->SetOrigin(origin);

  m_PSFImageMinMaxFilter = MinMaxType::New();
  m_PSFImageMinMaxFilter->SetImage(m_GibsonLanniPSFSource->GetOutput());
  region = m_GibsonLanniPSFSource->GetOutput()->GetLargestPossibleRegion();
  m_PSFImageMinMaxFilter->SetRegion(region);
  m_PSFImageMinMaxFilter->Compute();

  m_PSFImageITKToVTKFilter->SetInput(m_GibsonLanniPSFSource->GetOutput());
  m_PSFImageITKToVTKFilter->Modified();
  m_PSFImageITKToVTKFilter->Update();

  // Set up cost function
  m_CostFunction->SetFixedImage(m_MeasuredImageData);
  m_CostFunction->SetMovingImageSource(m_GibsonLanniPSFSource);

}


void
DataModel
::SavePSFImageFile(std::string fileName) {
  TIFFScaleType::Pointer scaler = TIFFScaleType::New();
  scaler->SetInput(m_GibsonLanniPSFSource->GetOutput());
  double min = GetPSFImageDataMinimum();
  double max = GetPSFImageDataMaximum();
  scaler->SetShift(-min);
  scaler->SetScale(65535.0f * (max - min));

  TIFFWriterType::Pointer writer = TIFFWriterType::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInput(scaler->GetOutput());
  writer->Update();
}


std::string
DataModel
::GetMeasuredImageFileName() {
  return m_ImageFileName;
}


void
DataModel
::SetNumberOfThreads(int threads) {
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(threads);
  itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads);
}


int
DataModel
::GetNumberOfThreads() {
  return itk::MultiThreader::GetGlobalMaximumNumberOfThreads();
}


void 
DataModel
::SetMeasuredImageData(TImage::Pointer image) {
  m_MeasuredImageData = image;

  // Set image data.
  m_MeasuredImageITKToVTKFilter->SetInput(m_MeasuredImageData);
}


DataModel::TImage::Pointer
DataModel
::GetMeasuredImageData() {
  return m_MeasuredImageData;
}


vtkAlgorithmOutput*
DataModel
::GetMeasuredImageOutputPort() {
  return m_MeasuredImageITKToVTKFilter->GetOutputPort();
}


vtkAlgorithmOutput*
DataModel
::GetPSFImageOutputPort() {
  return m_PSFImageITKToVTKFilter->GetOutputPort();
}


double
DataModel
::GetMeasuredImageDataMinimum() {
  if (!GetMeasuredImageData())
    return 0.0;

  return m_MeasuredImageMinMaxFilter->GetMinimum();
}


double
DataModel
::GetMeasuredImageDataMaximum() {
  if (!GetMeasuredImageData())
    return 0.0;

  return m_MeasuredImageMinMaxFilter->GetMaximum();
}


void
DataModel
::GetMeasuredImageDimensions(int dimensions[3]) {
  if (!GetMeasuredImageData()) {
    for (int i = 0; i < 3; i++)
      dimensions[i] = 0;
    return;
  }

  Float3DImageType::RegionType region 
      = GetMeasuredImageData()->GetLargestPossibleRegion();
  itk::Size<3> size = region.GetSize();

  for (int i = 0; i < 3; i++)
    dimensions[i] = size[i];
}


void
DataModel
::SetMeasuredImageVoxelSpacing(double spacing[3]) {
  float floatSpacing[3];
  for (int i = 0; i < 3; i++)
    floatSpacing[i] = static_cast<float>(spacing[i]);

  m_MeasuredImageData->SetSpacing(floatSpacing);
  m_MeasuredImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetMeasuredImageVoxelSpacing(int dimension, double spacing) {
  if (!m_MeasuredImageData)
    return;
  
  TImage::SpacingType currentSpacing = m_MeasuredImageData->GetSpacing();
  currentSpacing[dimension] = spacing;

  double doubleSpacing[3];
  for (int i = 0; i < 3; i++)
    doubleSpacing[i] = currentSpacing[i];

  SetMeasuredImageVoxelSpacing(doubleSpacing);

  m_MeasuredImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetMeasuredImageVoxelXSpacing(double spacing) {
  SetMeasuredImageVoxelSpacing(0, spacing); 
}


void
DataModel
::SetMeasuredImageVoxelYSpacing(double spacing) {
  SetMeasuredImageVoxelSpacing(1, spacing); 
}


void
DataModel
::SetMeasuredImageVoxelZSpacing(double spacing) {
  SetMeasuredImageVoxelSpacing(2, spacing); 
}


void
DataModel
::GetMeasuredImageVoxelSpacing(double spacing[3]) {
  if (!GetMeasuredImageData()) {
    for (int i = 0; i < 3; i++)
      spacing[i] = 0;
    return;
  }

  itk::Vector<double> thisSpacing = GetMeasuredImageData()->GetSpacing();
  for (int i = 0; i < 3; i++)
    spacing[i] = thisSpacing[i];
}


void
DataModel
::SetMeasuredImageOrigin(double origin[3]) {
  m_MeasuredImageData->SetOrigin(origin);
}


void
DataModel
::GetMeasuredImageOrigin(double origin[3]) {
  for (int i = 0; i < 3; i++) {
    origin[i] = m_MeasuredImageData->GetOrigin()[i];
  }
}


double
DataModel
::GetPSFImageDataMinimum() {
  if (!GetMeasuredImageData()) {
    return 0.0;
  }

  Float3DImageType::RegionType region = m_GibsonLanniPSFSource->GetOutput()
    ->GetLargestPossibleRegion();
  m_PSFImageMinMaxFilter->SetRegion(region);
  m_PSFImageMinMaxFilter->Compute();
  return m_PSFImageMinMaxFilter->GetMinimum();
}


double
DataModel
::GetPSFImageDataMaximum() {
  if (!GetMeasuredImageData()) {
    return 0.0;
  }

  Float3DImageType::RegionType region = m_GibsonLanniPSFSource->GetOutput()
    ->GetLargestPossibleRegion();
  m_PSFImageMinMaxFilter->SetRegion(region);
  m_PSFImageMinMaxFilter->Compute();
  return m_PSFImageMinMaxFilter->GetMaximum();  
}


void
DataModel
::SetPSFImageDimensions(int dimensions[3]) {
  unsigned long ulDimensions[3];
  for (int i = 0; i < 3; i++)
    ulDimensions[i] = static_cast<unsigned long>(dimensions[i]);

  m_GibsonLanniPSFSource->SetSize(ulDimensions);
}


void
DataModel
::SetPSFImageDimension(int index, int dimension) {
  int dimensions[3];
  GetPSFImageDimensions(dimensions);
  if (index >= 0 && index <= 2) {
    dimensions[index] = dimension;
    SetPSFImageDimensions(dimensions);  
  }
}


void
DataModel
::SetPSFImageXDimension(int dimension) {
  SetPSFImageDimension(0, dimension);
}


void
DataModel
::SetPSFImageYDimension(int dimension) {
  SetPSFImageDimension(1, dimension);
}


void
DataModel
::SetPSFImageZDimension(int dimension) {
  SetPSFImageDimension(2, dimension);
}


void
DataModel
::GetPSFImageDimensions(int dimensions[3]) {
  Float3DImageType::RegionType region 
      = GetMeasuredImageData()->GetLargestPossibleRegion();
  itk::Size<3> size = region.GetSize();

  for (int i = 0; i < 3; i++)
    dimensions[i] = size[i];
}


void
DataModel
::SetPSFImageVoxelSpacing(double spacing[3]) {
  float floatSpacing[3];
  for (int i = 0; i < 3; i++)
    floatSpacing[i] = static_cast<float>(spacing[i]);

  m_GibsonLanniPSFSource->SetSpacing(floatSpacing);
  m_PSFImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetPSFImageVoxelSpacing(int dimension, double spacing) {
  if (!GetMeasuredImageData())
    return;
  
  TImage::SpacingType currentSpacing = m_MeasuredImageData->GetSpacing();
  currentSpacing[dimension] = spacing;

  double doubleSpacing[3];
  for (int i = 0; i < 3; i++)
    doubleSpacing[i] = currentSpacing[i];

  SetPSFImageVoxelSpacing(doubleSpacing);
}


void
DataModel
::SetPSFImageVoxelXSpacing(double spacing) {
  SetPSFImageVoxelSpacing(0, spacing); 
}


void
DataModel
::SetPSFImageVoxelYSpacing(double spacing) {
  SetPSFImageVoxelSpacing(1, spacing); 
}


void
DataModel
::SetPSFImageVoxelZSpacing(double spacing) {
  SetPSFImageVoxelSpacing(2, spacing); 
}


void
DataModel
::GetPSFImageVoxelSpacing(double spacing[3]) {
  if (!GetMeasuredImageData()) {
    for (int i = 0; i < 3; i++)
      spacing[i] = 0.0;
    return;
  }

  itk::Vector<double> thisSpacing = GetMeasuredImageData()->GetSpacing();
  for (int i = 0; i < 3; i++)
    spacing[i] = thisSpacing[i];
}


void
DataModel
::SetPSFImageOrigin(double origin[3]) {
  float fOrigin[3];
  for (int i = 0; i < 3; i++)
    fOrigin[i] = static_cast<float>(origin[i]);
  m_GibsonLanniPSFSource->SetOrigin(fOrigin);
  m_GibsonLanniPSFSource->Modified();
}


void
DataModel
::GetPSFImageOrigin(double origin[3]) {
  float* fOrigin = m_GibsonLanniPSFSource->GetOrigin();
  for (int i = 0; i < 3; i++)
    origin[i] = static_cast<double>(fOrigin[i]);
}


void
DataModel
::RecenterPSFImageOrigin() {
  double spacing[3];
  GetPSFImageVoxelSpacing(spacing);

  int size[3];
  GetPSFImageDimensions(size);
  
  float origin[3];
  for (int i = 0; i < 3; i++)
    origin[i] = -spacing[i]*static_cast<float>(size[i]-1)*0.5;
  m_GibsonLanniPSFSource->SetOrigin(origin);
}


void
DataModel
::UpdateGibsonLanniPSFImage() {
  std::cout << m_GibsonLanniPSFSource << std::endl;
  m_GibsonLanniPSFSource->Update();
  m_PSFImageMinMaxFilter->Compute();

  double min = GetPSFImageDataMinimum();
  double max = GetPSFImageDataMaximum();
  std::cout << "Min: " << min << ", " << max << std::endl;
}


double
DataModel
::GetImageComparisonMetric() {
  return m_CostFunction->GetValue(m_GibsonLanniPSFSource->GetParameters());
}

#endif // _DATA_MODEL_CXX_
