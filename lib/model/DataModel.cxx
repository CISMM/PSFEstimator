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
  m_ImageData = NULL;
  
  m_GibsonLanniPSFSource = GibsonLanniPSFImageSourceType::New();

  m_MinMaxFilter = MinMaxType::New();
  m_InputImageITKToVTKFilter = new ITKImageToVTKImage<TImage>();

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
  delete m_InputImageITKToVTKFilter;
}


void
DataModel
::LoadImageFile(std::string fileName) {
  m_ImageFileName = fileName;
  ScalarFileReaderType::Pointer reader = ScalarFileReaderType::New();
  reader->SetFileName(fileName.c_str());
  reader->Update();
  SetImageData(reader->GetOutput());
  
#if 1
  float spacing[3];
  spacing[0] = 65.0f; spacing[1] = 65.0f; spacing[2] = 200.0f;
  m_GibsonLanniPSFSource->SetSpacing(spacing);

  float origin[3];
  origin[0] = -spacing[0]*63*0.5;
  origin[1] = -spacing[1]*63*0.5;
  origin[2] = -spacing[2]*63*0.5;
  m_GibsonLanniPSFSource->SetOrigin(origin);
  m_GibsonLanniPSFSource->Update();
  SetImageData(m_GibsonLanniPSFSource->GetOutput());
#endif

  // Connect this image data to the various pipelines.
  //m_MinMaxFilter->SetImage(m_ImageData);
  m_MinMaxFilter->SetImage(m_GibsonLanniPSFSource->GetOutput());
  m_MinMaxFilter->Compute();
  
  m_InputImageITKToVTKFilter->Modified();
  m_InputImageITKToVTKFilter->Update();
}


void
DataModel
::SavePSFImageFile(std::string fileName) {
  TIFFScaleType::Pointer scaler = TIFFScaleType::New();
  scaler->SetInput(m_ImageData);
  double min = GetImageDataMinimum();
  double max = GetImageDataMaximum();
  scaler->SetShift(-min);
  scaler->SetScale(65535.0f * (max - min));

  TIFFWriterType::Pointer writer = TIFFWriterType::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInput(scaler->GetOutput());
  writer->Update();
}


std::string
DataModel
::GetImageFileName() {
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
::SetImageData(TImage::Pointer image) {
  m_ImageData = image;

  // Set image data.
  m_InputImageITKToVTKFilter->SetInput(m_ImageData);
}


DataModel::TImage::Pointer
DataModel
::GetImageData() {
  return m_ImageData;
}


vtkAlgorithmOutput*
DataModel
::GetImageOutputPort() {
  return m_InputImageITKToVTKFilter->GetOutputPort();
}


double
DataModel
::GetImageDataMinimum() {
  return m_MinMaxFilter->GetMinimum();
}


double
DataModel
::GetImageDataMaximum() {
  return m_MinMaxFilter->GetMaximum();
}


void
DataModel
::GetDimensions(int dimensions[3]) {
  if (!GetImageData()) {
    dimensions[0] = 0;
    dimensions[1] = 0;
    dimensions[2] = 0;
    return;
  }

  Float3DImageType::RegionType region 
      = GetImageData()->GetLargestPossibleRegion();
  itk::Size<3> size = region.GetSize();

  dimensions[0] = size[0];
  dimensions[1] = size[1];
  dimensions[2] = size[2];
}


void
DataModel
::SetVoxelSpacing(double spacing[3]) {
  if (!m_ImageData)
    return;

  m_ImageData->SetSpacing(spacing);

  m_InputImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetVoxelSpacing(int dimension, double spacing) {
  if (!m_ImageData)
    return;

  TImage::SpacingType currentSpacing = m_ImageData->GetSpacing();
  currentSpacing[dimension] = spacing;
  m_ImageData->SetSpacing(currentSpacing);

  m_InputImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetVoxelXSpacing(double spacing) {
  SetVoxelSpacing(0, spacing); 
}


void
DataModel
::SetVoxelYSpacing(double spacing) {
  SetVoxelSpacing(1, spacing); 
}


void
DataModel
::SetVoxelZSpacing(double spacing) {
  SetVoxelSpacing(2, spacing); 
}


void
DataModel
::GetVoxelSpacing(double spacing[3]) {
  if (!GetImageData()) {
    spacing[0] = 0;
    spacing[1] = 0;
    spacing[2] = 0;
    return;
  }

  itk::Vector<double> thisSpacing = GetImageData()->GetSpacing();
  spacing[0] = thisSpacing[0];
  spacing[1] = thisSpacing[1];
  spacing[2] = thisSpacing[2];
}


void
DataModel
::UpdateGibsonLanniPSFImage() {
  m_GibsonLanniPSFSource->Update();
}


#endif // _DATA_MODEL_CXX_
