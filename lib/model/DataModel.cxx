#ifndef _DATA_MODEL_CXX_
#define _DATA_MODEL_CXX_

#if defined(_WIN32) // Turn off deprecation warnings in Visual Studio
#pragma warning( disable : 4996 )
#endif

#include <cstdlib>

#include <itkMultiThreader.h>
#include <itkPoint.h>
#include <itkImageFileWriter.h>

#include <DataModel.h>


DataModel
::DataModel() {
  m_MeasuredImageData = NULL;
  
  m_GibsonLanniPSFSource = GibsonLanniPSFImageSourceType::New();
  m_GibsonLanniBSFSource = GibsonLanniBSFImageSourceType::New();

  m_MeasuredImageMinMaxFilter = MinMaxType::New();
  m_PSFImageMinMaxFilter      = MinMaxType::New();
  m_BSFImageMinMaxFilter      = MinMaxType::New();

  m_MeasuredImageITKToVTKFilter = new ITKImageToVTKImage<TImage>();
  m_PSFImageITKToVTKFilter      = new ITKImageToVTKImage<TImage>();
  m_BSFImageITKToVTKFilter      = new ITKImageToVTKImage<TImage>();

  m_ImageToImageCostFunction = ImageToImageCostFunctionType::New();
  m_CostFunction = ParameterizedCostFunctionType::New();
  m_CostFunction->SetImageToImageMetric(m_ImageToImageCostFunction);
  m_CostFunction->SetMovingImageSource(m_GibsonLanniBSFSource);

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
  delete m_BSFImageITKToVTKFilter;
}


void
DataModel
::CreateImageFile(int xSize, int ySize, int zSize,
                  float xSpacing, float ySpacing, float zSpacing) {

  DummyImageSourcePointer dummy = DummyImageSourceType::New();

  DummyImageSourceType::SizeType dummySize;
  dummySize[0] = xSize;
  dummySize[1] = ySize;
  dummySize[2] = zSize;
  dummy->SetSize(dummySize);

  DummyImageSourceType::SpacingType dummySpacing;
  dummySpacing[0] = xSpacing;
  dummySpacing[1] = ySpacing;
  dummySpacing[2] = zSpacing;
  dummy->SetSpacing(dummySpacing);
  
  dummy->SetScale(0.0);
  dummy->Update();
  SetMeasuredImageData(dummy->GetOutput());
  
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
  SetBSFImageVoxelSpacing(spacing);

  int size[3];
  GetMeasuredImageDimensions(size);
  SetPSFImageDimensions(size);
  SetBSFImageDimensions(size);
  
  float origin[3];
  for (int i = 0; i < 3; i++)
    origin[i] = -spacing[i]*static_cast<float>(size[i])*0.5;
  m_GibsonLanniPSFSource->SetOrigin(origin);
  m_GibsonLanniBSFSource->SetOrigin(origin);
  m_MeasuredImageData->SetOrigin(origin);

  m_PSFImageMinMaxFilter = MinMaxType::New();
  m_PSFImageMinMaxFilter->SetImage(m_GibsonLanniPSFSource->GetOutput());
  m_PSFImageITKToVTKFilter->SetInput(m_GibsonLanniPSFSource->GetOutput());  

  m_BSFImageMinMaxFilter = MinMaxType::New();
  m_BSFImageMinMaxFilter->SetImage(m_GibsonLanniBSFSource->GetOutput());
  m_BSFImageITKToVTKFilter->SetInput(m_GibsonLanniBSFSource->GetOutput());  

  // Set up cost function
  m_CostFunction->SetFixedImage(m_MeasuredImageData);
  m_CostFunction->SetMovingImageSource(m_GibsonLanniBSFSource);

  // Set up optimizer, but don't connect it to the cost function just yet.
  m_Optimizer = OptimizerType::New();  
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
  SetBSFImageVoxelSpacing(spacing);

  int size[3];
  GetMeasuredImageDimensions(size);
  SetPSFImageDimensions(size);
  SetBSFImageDimensions(size);
  
  float origin[3];
  for (int i = 0; i < 3; i++)
    origin[i] = -spacing[i]*static_cast<float>(size[i])*0.5;
  m_GibsonLanniPSFSource->SetOrigin(origin);
  m_GibsonLanniBSFSource->SetOrigin(origin);
  m_MeasuredImageData->SetOrigin(origin);

  // Set the shifting and scaling for the BSF source to that of the measured image
  m_GibsonLanniBSFSource->SetBackgroundIntensity(m_MeasuredImageMinMaxFilter->GetMinimum());
  m_GibsonLanniBSFSource->SetMaximumIntensity(m_MeasuredImageMinMaxFilter->GetMaximum());

  m_PSFImageMinMaxFilter = MinMaxType::New();
  m_PSFImageMinMaxFilter->SetImage(m_GibsonLanniPSFSource->GetOutput());
  m_PSFImageITKToVTKFilter->SetInput(m_GibsonLanniPSFSource->GetOutput());  

  m_BSFImageMinMaxFilter = MinMaxType::New();
  m_BSFImageMinMaxFilter->SetImage(m_GibsonLanniBSFSource->GetOutput());
  m_BSFImageITKToVTKFilter->SetInput(m_GibsonLanniBSFSource->GetOutput());  

  // Set up cost function
  m_CostFunction->SetFixedImage(m_MeasuredImageData);

  // Set up optimizer, but don't connect it to the cost function just yet.
  m_Optimizer = OptimizerType::New();
}


void
DataModel
::SavePSFImageFile(std::string fileName) {
  TIFFScaleType::Pointer scaler = TIFFScaleType::New();
  scaler->SetInput(m_GibsonLanniPSFSource->GetOutput());
  double min = GetPSFImageDataMinimum();
  double max = GetPSFImageDataMaximum();
  scaler->SetShift(-min);
  scaler->SetScale(65535.0f / (max - min));

  TIFFWriterType::Pointer writer = TIFFWriterType::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInput(scaler->GetOutput());
  writer->Update();
}


void
DataModel
::SaveBSFImageFile(std::string fileName) {
  TIFFScaleType::Pointer scaler = TIFFScaleType::New();
  scaler->SetInput(m_GibsonLanniBSFSource->GetOutput());
  double min = GetBSFImageDataMinimum();
  double max = GetBSFImageDataMaximum();
  scaler->SetShift(-min);
  scaler->SetScale(65535.0f / (max - min));

  TIFFWriterType::Pointer writer = TIFFWriterType::New();
  writer->SetFileName(fileName.c_str());
  writer->SetInput(scaler->GetOutput());
  writer->Update();
}


void
DataModel
::SetConfiguration(Configuration & c) {
  // Read the settings from the configuration structure
  std::string sec("FileInfo");

  std::string fileName = c.GetValue(sec, "FileName");
  LoadImageFile(fileName);

  sec = std::string("GibsonLanniPSFSettings");
  
  double vec3[3];
  c.GetValueAsDoubleArray(sec, "VoxelSpacing", vec3, 3);
  SetMeasuredImageVoxelSpacing(vec3);
  SetPSFImageVoxelSpacing(vec3);

  // Set up origin so that (0, 0, 0) is centered in the image volume.
  int dimensions[3];
  double origin[3];
  GetPSFImageDimensions(dimensions);
  for (int i = 0; i < 3; i++) {
    origin[i] = -0.5*static_cast<double>(dimensions[i]-1)*vec3[i];
  }
  SetMeasuredImageOrigin(origin);
  SetPSFImageOrigin(origin);

  c.GetValueAsDoubleArray(sec, "CCDBorderWidth", vec3, 2);
  SetCCDBorderWidth(vec3);

  c.GetValueAsDoubleArray(sec, "BeadCenter", vec3, 3);
  SetBSFPointCenter(vec3);

  double beadRadius = c.GetValueAsDouble(sec, "BeadRadius");
  SetBeadRadius(beadRadius);

  SetGLNumericalAperture(c.GetValueAsFloat(sec, "NumericalAperture"));
  SetGLMagnification(c.GetValueAsFloat(sec, "Magnification"));
  SetGLDesignCoverSlipRefractiveIndex
    (c.GetValueAsFloat(sec, "DesignCoverSlipRefractiveIndex"));
  SetGLActualCoverSlipRefractiveIndex
    (c.GetValueAsFloat(sec, "ActualCoverSlipRefractiveIndex"));
  SetGLDesignCoverSlipThickness
    (c.GetValueAsFloat(sec, "DesignCoverSlipThickness"));
  SetGLActualCoverSlipThickness
    (c.GetValueAsFloat(sec, "ActualCoverSlipThickness"));
  SetGLDesignImmersionOilRefractiveIndex
    (c.GetValueAsFloat(sec, "DesignImmersionOilRefractiveIndex"));
  SetGLActualImmersionOilRefractiveIndex
    (c.GetValueAsFloat(sec, "ActualImmersionOilRefractiveIndex"));
  SetGLDesignImmersionOilThickness
    (c.GetValueAsFloat(sec, "DesignImmersionOilThickness"));
  SetGLDesignSpecimenLayerRefractiveIndex
    (c.GetValueAsFloat(sec, "DesignSpecimenLayerRefractiveIndex"));
  SetGLActualSpecimenLayerRefractiveIndex
    (c.GetValueAsFloat(sec, "ActualSpecimenLayerRefractiveIndex"));
  SetGLActualPointSourceDepthInSpecimenLayer
    (c.GetValueAsFloat(sec, "ActualPointSourceDepthInSpecimenLayer"));
  SetGLDesignDistanceFromBackFocalPlaneToDetector
    (c.GetValueAsFloat(sec, "DesignDistanceFromBackFocalPlaneToDetector"));
  SetGLActualDistanceFromBackFocalPlaneToDetector
    (c.GetValueAsFloat(sec, "ActualDistanceFromBackFocalPlaneToDetector"));

  sec = std::string("ZSliceCoordinates");
  
  SetUseCustomZCoordinates(c.GetValueAsBool(sec, "UseCustomZCoordinates"));

  for (unsigned int i = 0; i < m_GibsonLanniBSFSource->GetSize()[2]; i++) {
    char name[128];
    sprintf(name, "ZCoordinate%d", i);
    SetZCoordinate(i, c.GetValueAsDouble(sec, name));
  }
}


void
DataModel
::GetConfiguration(Configuration & c) {
  // Dump the settings into the configuration structure
  std::string sec("FileInfo");
  c.SetValue(sec, "FileName", m_ImageFileName);

  sec = std::string("GibsonLanniPSFSettings");

  double vec3[3];
  GetPSFImageVoxelSpacing(vec3);
  c.SetValueFromDoubleArray(sec, "VoxelSpacing", vec3, 3);

  GetCCDBorderWidth(vec3);
  c.SetValueFromDoubleArray(sec, "CCDBorderWidth", vec3, 2);

  GetPSFPointCenter(vec3);
  c.SetValueFromDoubleArray(sec, "BeadCenter", vec3, 3);

  c.SetValueFromDouble(sec, "BeadRadius", GetBeadRadius());

  c.SetValueFromFloat(sec, "NumericalAperture",
		      GetGLNumericalAperture());
  c.SetValueFromFloat(sec, "Magnification",
		      GetGLMagnification());
  c.SetValueFromFloat(sec, "DesignCoverSlipRefractiveIndex",
		      GetGLDesignCoverSlipRefractiveIndex());
  c.SetValueFromFloat(sec, "ActualCoverSlipRefractiveIndex",
		      GetGLActualCoverSlipRefractiveIndex());
  c.SetValueFromFloat(sec, "DesignCoverSlipThickness",
		      GetGLDesignCoverSlipThickness());
  c.SetValueFromFloat(sec, "ActualCoverSlipThickness",
		      GetGLActualCoverSlipThickness());
  c.SetValueFromFloat(sec, "DesignImmersionOilRefractiveIndex",
		      GetGLDesignImmersionOilRefractiveIndex());
  c.SetValueFromFloat(sec, "ActualImmersionOilRefractiveIndex",
		      GetGLActualImmersionOilRefractiveIndex());
  c.SetValueFromFloat(sec, "DesignImmersionOilThickness",
		      GetGLDesignImmersionOilThickness());
  c.SetValueFromFloat(sec, "DesignSpecimenLayerRefractiveIndex",
		      GetGLDesignSpecimenLayerRefractiveIndex());
  c.SetValueFromFloat(sec, "ActualSpecimenLayerRefractiveIndex",
		      GetGLActualSpecimenLayerRefractiveIndex());
  c.SetValueFromFloat(sec, "ActualPointSourceDepthInSpecimenLayer",
		      GetGLActualPointSourceDepthInSpecimenLayer());
  c.SetValueFromFloat(sec, "DesignDistanceFromBackFocalPlaneToDetector",
		      GetGLDesignDistanceFromBackFocalPlaneToDetector());
  c.SetValueFromFloat(sec, "ActualDistanceFromBackFocalPlaneToDetector",
		      GetGLActualDistanceFromBackFocalPlaneToDetector());

  sec = std::string("ZSliceCoordinates");
  
  c.SetValueFromBool(sec, "UseCustomZCoordinates", GetUseCustomZCoordinates());

  for (unsigned int i = 0; i < m_GibsonLanniBSFSource->GetSize()[2]; i++) {
    char name[128];
    sprintf(name, "ZCoordinate%d", i);
    c.SetValueFromDouble(sec, name, GetZCoordinate(i));
  }
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


vtkAlgorithmOutput*
DataModel
::GetBSFImageOutputPort() {
  return m_BSFImageITKToVTKFilter->GetOutputPort();
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


DataModel::Float3DPointType
DataModel
::GetMeasuredImageDataMaximumCoordinates() {
  Float3DPointType point;
  Float3DImageType::Pointer image = GetMeasuredImageData();
  if (image) {
    image->TransformIndexToPhysicalPoint(m_MeasuredImageMinMaxFilter->GetIndexOfMaximum(), point);
  }

  return point;
}


int
DataModel
::GetNumberOfProperties() {
  return m_GibsonLanniBSFSource->GetNumberOfParameters();
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

  m_GibsonLanniPSFSource->UpdateLargestPossibleRegion();

  m_PSFImageMinMaxFilter->Compute();
  return m_PSFImageMinMaxFilter->GetMinimum();
}


double
DataModel
::GetPSFImageDataMaximum() {
  if (!GetMeasuredImageData()) {
    return 0.0;
  }

  m_GibsonLanniPSFSource->UpdateLargestPossibleRegion();

  m_PSFImageMinMaxFilter->Compute();
  return m_PSFImageMinMaxFilter->GetMaximum();  
}


double
DataModel
::GetBSFImageDataMinimum() {
  if (!GetMeasuredImageData()) {
    return 0.0;
  }

  m_GibsonLanniBSFSource->UpdateLargestPossibleRegion();

  m_BSFImageMinMaxFilter->Compute();
  return m_BSFImageMinMaxFilter->GetMinimum();
}


double
DataModel
::GetBSFImageDataMaximum() {
  if (!GetMeasuredImageData()) {
    return 0.0;
  }

  m_GibsonLanniBSFSource->UpdateLargestPossibleRegion();

  m_BSFImageMinMaxFilter->Compute();
  return m_BSFImageMinMaxFilter->GetMaximum();  
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
::GetPSFImageDimensions(int dimensions[3]) {
  if (GetMeasuredImageData()) {
    Float3DImageType::RegionType region 
      = GetMeasuredImageData()->GetLargestPossibleRegion();
    itk::Size<3> size = region.GetSize();

    for (int i = 0; i < 3; i++)
      dimensions[i] = size[i];
  } else {
    for (int i = 0; i < 3; i++)
      dimensions[i] = 0;
  }
}


void
DataModel
::SetBSFImageDimensions(int dimensions[3]) {
  unsigned long ulDimensions[3];
  for (int i = 0; i < 3; i++)
    ulDimensions[i] = static_cast<unsigned long>(dimensions[i]);

  m_GibsonLanniBSFSource->SetSize(ulDimensions);
}


void
DataModel
::SetBSFImageDimension(int index, int dimension) {
  int dimensions[3];
  GetBSFImageDimensions(dimensions);
  if (index >= 0 && index <= 2) {
    dimensions[index] = dimension;
    SetBSFImageDimensions(dimensions);  
  }
}


void
DataModel
::GetBSFImageDimensions(int dimensions[3]) {
  if (GetMeasuredImageData()) {
    Float3DImageType::RegionType region 
      = GetMeasuredImageData()->GetLargestPossibleRegion();
    itk::Size<3> size = region.GetSize();

    for (int i = 0; i < 3; i++)
      dimensions[i] = size[i];

  } else {
    for (int i = 0; i < 3; i++)
      dimensions[i] = 0;
  }
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
::SetBSFImageVoxelSpacing(double spacing[3]) {
  float floatSpacing[3];
  for (int i = 0; i < 3; i++)
    floatSpacing[i] = static_cast<float>(spacing[i]);

  m_GibsonLanniBSFSource->SetSpacing(floatSpacing);
  m_BSFImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetBSFImageVoxelSpacing(int dimension, double spacing) {
  if (!GetMeasuredImageData())
    return;
  
  TImage::SpacingType currentSpacing = m_MeasuredImageData->GetSpacing();
  currentSpacing[dimension] = spacing;

  double doubleSpacing[3];
  for (int i = 0; i < 3; i++)
    doubleSpacing[i] = currentSpacing[i];

  SetBSFImageVoxelSpacing(doubleSpacing);
}


void
DataModel
::GetBSFImageVoxelSpacing(double spacing[3]) {
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
::SetCCDBorderWidth(double borderWidth[2]) {
  float width[2];
  width[0] = static_cast<float>(borderWidth[0]);
  width[1] = static_cast<float>(borderWidth[1]);
  m_GibsonLanniPSFSource->SetCCDBorderWidth(width);
}


void
DataModel
::GetCCDBorderWidth(double borderWidth[2]) {
  float* width =  m_GibsonLanniPSFSource->GetCCDBorderWidth();
  borderWidth[0] = static_cast<double>(width[0]);
  borderWidth[1] = static_cast<double>(width[1]);
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
::SetBSFImageOrigin(double origin[3]) {
  float fOrigin[3];
  for (int i = 0; i < 3; i++)
    fOrigin[i] = static_cast<float>(origin[i]);
  m_GibsonLanniBSFSource->SetOrigin(fOrigin);
  m_GibsonLanniBSFSource->Modified();
}


void
DataModel
::GetBSFImageOrigin(double origin[3]) {
  float* fOrigin = m_GibsonLanniBSFSource->GetOrigin();
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
::SetPSFPointCenter(double center[3]) {
  float fCenter[3];
  for (int i = 0; i < 3; i++)
    fCenter[i] = static_cast<float>(center[i]);
  m_GibsonLanniPSFSource->SetPointCenter(fCenter);
  m_GibsonLanniBSFSource->SetBeadCenter(fCenter);
}


void
DataModel
::GetPSFPointCenter(double center[3]) {
  float* fCenter = m_GibsonLanniPSFSource->GetPointCenter();
  for (int i = 0; i < 3; i++)
    center[i] = static_cast<double>(fCenter[i]);
}


void
DataModel
::SetBSFPointCenter(double center[3]) {
  float fCenter[3];
  for (int i = 0; i < 3; i++)
    fCenter[i] = static_cast<float>(center[i]);
  m_GibsonLanniPSFSource->SetPointCenter(fCenter);
}


void
DataModel
::GetBSFPointCenter(double center[3]) {
  float* fCenter = m_GibsonLanniBSFSource->GetBeadCenter();
  for (int i = 0; i < 3; i++)
    center[i] = static_cast<double>(fCenter[i]);
}


void
DataModel
::UpdateGibsonLanniPSFImage() {
  m_GibsonLanniPSFSource->UpdateLargestPossibleRegion();
  m_PSFImageMinMaxFilter->Compute();
}

void
DataModel
::UpdateGibsonLanniBSFImage() {
  m_GibsonLanniBSFSource->UpdateLargestPossibleRegion();
  m_BSFImageMinMaxFilter->Compute();
}


void
DataModel
::SetBeadRadius(double radius) {
  m_GibsonLanniBSFSource->SetBeadRadius(radius);
}


double
DataModel
::GetBeadRadius() {
  return m_GibsonLanniBSFSource->GetBeadRadius();
}


void
DataModel
::SetGLParameterEnabled(unsigned int index, bool enabled) {
  try {
    typedef ParameterizedCostFunctionType::ParametersMaskType
      ParametersMaskType;
    ParametersMaskType* parametersMask = m_CostFunction->GetParametersMask();
    if (index < parametersMask->Size()) {
      parametersMask->SetElement(index, enabled ? 1 : 0);
    }
  } catch (...) {}
}


bool
DataModel
::GetGLParameterEnabled(unsigned int index) {
  bool enabled = false;
  try {
    typedef ParameterizedCostFunctionType::ParametersMaskType
      ParametersMaskType;
    ParametersMaskType* parametersMask = m_CostFunction->GetParametersMask();
    if (index < parametersMask->Size()) {
      enabled = parametersMask->GetElement(index) != 0;
    }
  } catch (...) {}

  return enabled;
}


double
DataModel
::GetImageComparisonMetric() {
  return m_CostFunction->GetValue(m_GibsonLanniBSFSource->GetParameters());
}


void
DataModel
::Optimize() {
  typedef ParameterizedCostFunctionType::ParametersMaskType
    ParametersMaskType;
  ParametersMaskType* mask = m_CostFunction->GetParametersMask();

  // Pluck out the active parameters
  typedef ParameterizedCostFunctionType::ParametersType ParametersType;
  ParametersType activeParameters
    = ParametersType(m_CostFunction->GetNumberOfParameters());
  int activeIndex = 0;
  for (unsigned int i = 0; i < mask->Size(); i++) {
    if (mask->GetElement(i)) {
      activeParameters[activeIndex++] = m_GibsonLanniBSFSource->GetParameters()[i];
    }
  }

  // Connect to the cost function, set the initial parameters, and optimize.
  m_ImageToImageCostFunction
    ->SetFixedImageRegion(m_GibsonLanniBSFSource->GetOutput()->GetLargestPossibleRegion());
  m_Optimizer->SetCostFunction(m_CostFunction);
  m_Optimizer->SetFunctionConvergenceTolerance(1e-3);
  m_Optimizer->SetInitialPosition(activeParameters);
  m_Optimizer->StartOptimization();
}

#endif // _DATA_MODEL_CXX_
