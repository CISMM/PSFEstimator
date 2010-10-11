#ifndef _DATA_MODEL_CXX_
#define _DATA_MODEL_CXX_

#if defined(_WIN32) // Turn off deprecation warnings in Visual Studio
#pragma warning( disable : 4996 )
#endif

#include <cfloat>
#include <cstdlib>

#include <itkMultiThreader.h>
#include <itkPoint.h>
#include <itkImageFileWriter.h>

#include <itkGibsonLanniBSFImageSource.txx>
#include <itkGibsonLanniPSFImageSource.txx>
#include <itkGridImageSource.txx>
#include <itkImageFileReader.txx>
#include <itkImageFileWriter.txx>
#include <itkImageToParametricImageSourceMetric.txx>
#include <itkMeanSquaresImageToImageMetric.txx>
#include <itkMinimumMaximumImageCalculator.txx>
//#include <itkPoissonNoiseImageToImageMetric.cxx>
#include <itkShiftScaleImageFilter.txx>

#include <ITKImageToVTKImage.cxx>

#include <vtkAlgorithm.h>

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
  m_CostFunction = ParametricCostFunctionType::New();
  m_CostFunction->SetInterpolator(InterpolatorType::New());
  m_CostFunction->SetDelegateMetric(m_ImageToImageCostFunction);
  m_CostFunction->SetMovingImageSource(m_GibsonLanniBSFSource);

  SetInitialSimplexDeltas();

  // ITK will detect the number of cores on the system and set it by default.
  // Here we can override that setting if the proper environment variable is
  // set.
  char *var = getenv("PSFEstimator_THREADS");
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


bool
DataModel
::LoadSessionFile(const std::string& fileName) {
  Configuration config;
  config.Parse(fileName);

  // Read the settings from the configuration structure
  LoadImageFile(config.GetValue("FileInfo", "FileName"));

  SetConfiguration(config);

  return true;
}


bool
DataModel
::SaveSessionFile(const std::string& fileName) {
  Configuration config;
  GetConfiguration(config);
  std::ofstream os(fileName.c_str());
  config.Write(os);

  return true;
}


void
DataModel
::SetInitialSimplexDeltas() {
  m_InitialSimplexDelta.
    SetSize(m_GibsonLanniBSFSource->GetNumberOfParameters());

  int index = 0;
  m_InitialSimplexDelta[index++] = 1.0; // X-spacing
  m_InitialSimplexDelta[index++] = 1.0; // Y-spacing
  m_InitialSimplexDelta[index++] = 1.0; // Z-spacing
  m_InitialSimplexDelta[index++] = 10.0; // Bead radius
  m_InitialSimplexDelta[index++] = 100.0; // Bead center X
  m_InitialSimplexDelta[index++] = 100.0; // Bead center Y
  m_InitialSimplexDelta[index++] = 100.0; // Bead center Z
  m_InitialSimplexDelta[index++] = 0.1; // Shear X
  m_InitialSimplexDelta[index++] = 0.1; // Shear Y
  m_InitialSimplexDelta[index++] = 5.0; // Emission wavelength
  m_InitialSimplexDelta[index++] = 0.05; // NA
  m_InitialSimplexDelta[index++] = 1.0; // Magnification
  m_InitialSimplexDelta[index++] = 0.001; // Design cover slip RI
  m_InitialSimplexDelta[index++] = 0.001; // Actual cover slip RI
  m_InitialSimplexDelta[index++] = 1.0;   // Design cover slip thickness
  m_InitialSimplexDelta[index++] = 1.0;   // Actual cover slip thickness
  m_InitialSimplexDelta[index++] = 0.001; // Design immersion oil RI
  m_InitialSimplexDelta[index++] = 0.001; // Actual immersion oil RI
  m_InitialSimplexDelta[index++] = 1.0;   // Design immersion oil thickness
  m_InitialSimplexDelta[index++] = 0.001; // Design specimen layer RI
  m_InitialSimplexDelta[index++] = 0.001; // Actual specimen layer RI
  m_InitialSimplexDelta[index++] = 1.0;   // Actual point source depth
  m_InitialSimplexDelta[index++] = 1.0;   // Design distance from back focal plane
  m_InitialSimplexDelta[index++] = 1.0;   // Actual distance from back focal plane
  m_InitialSimplexDelta[index++] = 1.0;   // Background intensity
  m_InitialSimplexDelta[index++] = 1.0;   // Maximum intensity
}


void
DataModel
::CreateImageFile(int xSize, int ySize, int zSize,
                  double xSpacing, double ySpacing, double zSpacing) {

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

  PointType origin;
  for (int i = 0; i < 3; i++)
    origin[i] = -spacing[i]*static_cast<double>(size[i])*0.5;
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

  PointType origin;
  for (int i = 0; i < 3; i++)
    origin[i] = -spacing[i]*static_cast<double>(size[i])*0.5;
  m_GibsonLanniPSFSource->SetOrigin(origin);
  m_GibsonLanniBSFSource->SetOrigin(origin);
  m_MeasuredImageData->SetOrigin(origin);

  // Set the shifting and scaling for the BSF source to that of the measured image
  //m_GibsonLanniBSFSource->SetBackgroundIntensity(m_MeasuredImageMinMaxFilter->GetMinimum());
  //m_GibsonLanniBSFSource->SetMaximumIntensity(m_MeasuredImageMinMaxFilter->GetMaximum());

  m_PSFImageMinMaxFilter = MinMaxType::New();
  m_PSFImageMinMaxFilter->SetImage(m_GibsonLanniPSFSource->GetOutput());
  m_PSFImageITKToVTKFilter->SetInput(m_GibsonLanniPSFSource->GetOutput());

  m_BSFImageMinMaxFilter = MinMaxType::New();
  m_BSFImageMinMaxFilter->SetImage(m_GibsonLanniBSFSource->GetOutput());
  m_BSFImageITKToVTKFilter->SetInput(m_GibsonLanniBSFSource->GetOutput());

  // Set up cost function
  m_CostFunction->SetFixedImage(m_MeasuredImageData);
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
  std::string sec = std::string("GibsonLanniPSFSettings");

  double vec3[3];
  c.GetValueAsDoubleArray(sec, "VoxelSpacing", vec3, 3);
  SetMeasuredImageVoxelSpacing(vec3);
  SetPSFImageVoxelSpacing(vec3);
  SetBSFImageVoxelSpacing(vec3);

  // Set up origin so that (0, 0, 0) is centered in the image volume.
  int dimensions[3];
  double origin[3];
  GetPSFImageDimensions(dimensions);
  for (int i = 0; i < 3; i++) {
    origin[i] = -0.5*static_cast<double>(dimensions[i]-1)*vec3[i];
  }
  SetMeasuredImageOrigin(origin);
  SetPSFImageOrigin(origin);
  SetBSFImageOrigin(origin);

  c.GetValueAsDoubleArray(sec, "BeadCenter", vec3, 3);
  SetPSFPointCenter(vec3);
  SetBSFPointCenter(vec3);

  double beadRadius = c.GetValueAsDouble(sec, "BeadRadius", GetBeadRadius());
  SetBeadRadius(beadRadius);

  double shearX = c.GetValueAsDouble(sec, "ShearX");
  m_GibsonLanniBSFSource->SetShearX(shearX);
  double shearY = c.GetValueAsDouble(sec, "ShearY");
  m_GibsonLanniBSFSource->SetShearY(shearY);

  SetGLEmissionWavelength
    (c.GetValueAsDouble(sec, "EmissionWavelength", GetGLEmissionWavelength()));

  SetGLNumericalAperture
    (c.GetValueAsDouble(sec, "NumericalAperture", GetGLNumericalAperture()));
  SetGLMagnification
    (c.GetValueAsDouble(sec, "Magnification", GetGLMagnification()));
  SetGLDesignCoverSlipRefractiveIndex
    (c.GetValueAsDouble(sec, "DesignCoverSlipRefractiveIndex",
                       GetGLDesignCoverSlipRefractiveIndex()));
  SetGLActualCoverSlipRefractiveIndex
    (c.GetValueAsDouble(sec, "ActualCoverSlipRefractiveIndex",
                       GetGLActualCoverSlipRefractiveIndex()));
  SetGLDesignCoverSlipThickness
    (c.GetValueAsDouble(sec, "DesignCoverSlipThickness",
                       GetGLDesignCoverSlipThickness()));
  SetGLActualCoverSlipThickness
    (c.GetValueAsDouble(sec, "ActualCoverSlipThickness",
                       GetGLActualCoverSlipThickness()));
  SetGLDesignImmersionOilRefractiveIndex
    (c.GetValueAsDouble(sec, "DesignImmersionOilRefractiveIndex",
                       GetGLDesignImmersionOilRefractiveIndex()));
  SetGLActualImmersionOilRefractiveIndex
    (c.GetValueAsDouble(sec, "ActualImmersionOilRefractiveIndex",
                       GetGLActualImmersionOilRefractiveIndex()));
  SetGLDesignImmersionOilThickness
    (c.GetValueAsDouble(sec, "DesignImmersionOilThickness",
                       GetGLDesignImmersionOilThickness()));
  SetGLDesignSpecimenLayerRefractiveIndex
    (c.GetValueAsDouble(sec, "DesignSpecimenLayerRefractiveIndex",
                       GetGLDesignSpecimenLayerRefractiveIndex()));
  SetGLActualSpecimenLayerRefractiveIndex
    (c.GetValueAsDouble(sec, "ActualSpecimenLayerRefractiveIndex",
                       GetGLActualSpecimenLayerRefractiveIndex()));
  SetGLActualPointSourceDepthInSpecimenLayer
    (c.GetValueAsDouble(sec, "ActualPointSourceDepthInSpecimenLayer",
                       GetGLActualPointSourceDepthInSpecimenLayer()));
  SetGLDesignDistanceFromBackFocalPlaneToDetector
    (c.GetValueAsDouble(sec, "DesignDistanceFromBackFocalPlaneToDetector",
                       GetGLDesignDistanceFromBackFocalPlaneToDetector()));
  SetGLActualDistanceFromBackFocalPlaneToDetector
    (c.GetValueAsDouble(sec, "ActualDistanceFromBackFocalPlaneToDetector",
                       GetGLActualDistanceFromBackFocalPlaneToDetector()));
  SetGLBackgroundIntensity
    (c.GetValueAsDouble(sec, "BackgroundIntensity",
                       GetGLBackgroundIntensity()));
  SetGLMaximumIntensity
    (c.GetValueAsDouble(sec, "MaximumIntensity",
                       GetGLMaximumIntensity()));

  sec = std::string("ZSliceCoordinates");

  SetUseCustomZCoordinates(c.GetValueAsBool(sec, "UseCustomZCoordinates",
                                            GetUseCustomZCoordinates()));

  for (unsigned int i = 0; i < m_GibsonLanniBSFSource->GetSize()[2]; i++) {
    char name[128];
    sprintf(name, "ZCoordinate%d", i);
    SetZCoordinate(i, c.GetValueAsDouble(sec, name));
  }

  unsigned int index = 0;
  sec = std::string("Optimization");
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "VoxelSpacingX"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "VoxelSpacingY"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "VoxelSpacingZ"));

  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "BeadRadius"));

  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "BeadCenterX"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "BeadCenterY"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "BeadCenterZ"));

  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "ShearX"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "ShearY"));

  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "EmissionWavelength"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "NumericalAperture"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "Magnification"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "DesignCoverSlipRefractiveIndex"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "ActualCoverSlipRefractiveIndex"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "DesignCoverSlipThickness"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "ActualCoverSlipThickness"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "DesignImmersionOilRefractiveIndex"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "ActualImmersionOilRefractiveIndex"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "DesignImmersionOilThickness"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "DesignSpecimenLayerRefractiveIndex"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "ActualSpecimenLayerRefractiveIndex"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "ActualPointSourceDepthInSpecimenLayer"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "DesignDistanceFromBackFocalPlaneToDetector"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "ActualDistanceFromBackFocalPlaneToDetector"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "BackgroundIntensity"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "MaximumIntensity"));
}


void
DataModel
::GetConfiguration(Configuration & c) {
  // Dump the settings into the configuration structure
  std::string sec("FileInfo");
  c.SetValue(sec, "FileName", m_ImageFileName);

  sec = std::string("GibsonLanniPSFSettings");

  double vec3[3];
  GetBSFImageVoxelSpacing(vec3);
  c.SetValueFromDoubleArray(sec, "VoxelSpacing", vec3, 3);

  GetBSFPointCenter(vec3);
  c.SetValueFromDoubleArray(sec, "BeadCenter", vec3, 3);

  c.SetValueFromDouble(sec, "BeadRadius", GetBeadRadius());

  c.SetValueFromDouble(sec, "ShearX", m_GibsonLanniBSFSource->GetShearX());
  c.SetValueFromDouble(sec, "ShearY", m_GibsonLanniBSFSource->GetShearY());

  c.SetValueFromDouble(sec, "EmissionWavelength", m_GibsonLanniBSFSource->GetEmissionWavelength());

  c.SetValueFromDouble(sec, "NumericalAperture",
		      GetGLNumericalAperture());
  c.SetValueFromDouble(sec, "Magnification",
		      GetGLMagnification());
  c.SetValueFromDouble(sec, "DesignCoverSlipRefractiveIndex",
		      GetGLDesignCoverSlipRefractiveIndex());
  c.SetValueFromDouble(sec, "ActualCoverSlipRefractiveIndex",
		      GetGLActualCoverSlipRefractiveIndex());
  c.SetValueFromDouble(sec, "DesignCoverSlipThickness",
		      GetGLDesignCoverSlipThickness());
  c.SetValueFromDouble(sec, "ActualCoverSlipThickness",
		      GetGLActualCoverSlipThickness());
  c.SetValueFromDouble(sec, "DesignImmersionOilRefractiveIndex",
		      GetGLDesignImmersionOilRefractiveIndex());
  c.SetValueFromDouble(sec, "ActualImmersionOilRefractiveIndex",
		      GetGLActualImmersionOilRefractiveIndex());
  c.SetValueFromDouble(sec, "DesignImmersionOilThickness",
		      GetGLDesignImmersionOilThickness());
  c.SetValueFromDouble(sec, "DesignSpecimenLayerRefractiveIndex",
		      GetGLDesignSpecimenLayerRefractiveIndex());
  c.SetValueFromDouble(sec, "ActualSpecimenLayerRefractiveIndex",
		      GetGLActualSpecimenLayerRefractiveIndex());
  c.SetValueFromDouble(sec, "ActualPointSourceDepthInSpecimenLayer",
		      GetGLActualPointSourceDepthInSpecimenLayer());
  c.SetValueFromDouble(sec, "DesignDistanceFromBackFocalPlaneToDetector",
		      GetGLDesignDistanceFromBackFocalPlaneToDetector());
  c.SetValueFromDouble(sec, "ActualDistanceFromBackFocalPlaneToDetector",
		      GetGLActualDistanceFromBackFocalPlaneToDetector());
  c.SetValueFromDouble(sec, "BackgroundIntensity",
                      GetGLBackgroundIntensity());
  c.SetValueFromDouble(sec, "MaximumIntensity",
                      GetGLMaximumIntensity());

  sec = std::string("ZSliceCoordinates");

  c.SetValueFromBool(sec, "UseCustomZCoordinates", GetUseCustomZCoordinates());

  for (unsigned int i = 0; i < m_GibsonLanniBSFSource->GetSize()[2]; i++) {
    char name[128];
    sprintf(name, "ZCoordinate%d", i);
    c.SetValueFromDouble(sec, name, GetZCoordinate(i));
  }

  int index = 0;
  sec = std::string("Optimization");
  c.SetValueFromBool(sec, "VoxelSpacingX", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "VoxelSpacingY", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "VoxelSpacingZ", GetGLParameterEnabled(index++));

  c.SetValueFromBool(sec, "BeadRadius", GetGLParameterEnabled(index++));

  c.SetValueFromBool(sec, "BeadCenterX", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "BeadCenterY", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "BeadCenterZ", GetGLParameterEnabled(index++));

  c.SetValueFromBool(sec, "ShearX", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "ShearY", GetGLParameterEnabled(index++));

  c.SetValueFromBool(sec, "EmissionWavelength", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "NumericalAperture", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "Magnification", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "DesignCoverSlipRefractiveIndex", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "ActualCoverSlipRefractiveIndex", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "DesignCoverSlipThickness", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "ActualCoverSlipThickness", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "DesignImmersionOilRefractiveIndex", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "ActualImmersionOilRefractiveIndex", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "DesignImmersionOilThickness", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "DesignSpecimenLayerRefractiveIndex", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "ActualSpecimenLayerRefractiveIndex", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "ActualPointSourceDepthInSpecimenLayer", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "DesignDistanceFromBackFocalPlaneToDetector", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "ActualDistanceFromBackFocalPlaneToDetector", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "BackgroundIntensity", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "MaximumIntensity", GetGLParameterEnabled(index++));
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
  m_MeasuredImageData->SetSpacing(spacing);
  m_MeasuredImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetMeasuredImageVoxelSpacing(int dimension, double spacing) {
  if (!m_MeasuredImageData)
    return;

  SpacingType currentSpacing = m_MeasuredImageData->GetSpacing();
  currentSpacing[dimension] = spacing;

  SetMeasuredImageVoxelSpacing(currentSpacing.GetDataPointer());

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

  SpacingType thisSpacing = GetMeasuredImageData()->GetSpacing();
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
  SizeType size;
  for (int i = 0; i < 3; i++)
    size[i] = static_cast<SizeValueType>(dimensions[i]);

  m_GibsonLanniPSFSource->SetSize(size);
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
  SizeType size;
  for (int i = 0; i < 3; i++)
    size[i] = static_cast<SizeValueType>(dimensions[i]);

  m_GibsonLanniBSFSource->SetSize(size);
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
  SpacingType thisSpacing;
  for (int i = 0; i < 3; i++)
    thisSpacing[i] = static_cast<SpacingValueType>(spacing[i]);

  m_GibsonLanniPSFSource->SetSpacing(thisSpacing);
  m_PSFImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetPSFImageVoxelSpacing(int dimension, double spacing) {
  if (!GetMeasuredImageData())
    return;

  TImage::SpacingType currentSpacing = m_MeasuredImageData->GetSpacing();
  currentSpacing[dimension] = spacing;

  SetPSFImageVoxelSpacing(currentSpacing.GetDataPointer());
}


void
DataModel
::GetPSFImageVoxelSpacing(double spacing[3]) {
  if (!GetMeasuredImageData()) {
    for (int i = 0; i < 3; i++)
      spacing[i] = 0.0;
    return;
  }

  SpacingType thisSpacing = GetMeasuredImageData()->GetSpacing();
  for (int i = 0; i < 3; i++)
    spacing[i] = thisSpacing[i];
}


void
DataModel
::SetBSFImageVoxelSpacing(double spacing[3]) {
  SpacingType thisSpacing;
  for (int i = 0; i < 3; i++)
    thisSpacing[i] = static_cast<SpacingValueType>(spacing[i]);

  m_GibsonLanniBSFSource->SetSpacing(thisSpacing);
  m_BSFImageITKToVTKFilter->GetOutputPort()->GetProducer()->Modified();
}


void
DataModel
::SetBSFImageVoxelSpacing(int dimension, double spacing) {
  if (!GetMeasuredImageData())
    return;

  SpacingType currentSpacing = m_MeasuredImageData->GetSpacing();
  currentSpacing[dimension] = spacing;

  SetBSFImageVoxelSpacing(currentSpacing.GetDataPointer());
}


void
DataModel
::GetBSFImageVoxelSpacing(double spacing[3]) {
  if (!GetMeasuredImageData()) {
    for (int i = 0; i < 3; i++)
      spacing[i] = 0.0;
    return;
  }

  SpacingType thisSpacing = GetMeasuredImageData()->GetSpacing();
  for (int i = 0; i < 3; i++)
    spacing[i] = thisSpacing[i];
}


void
DataModel
::SetPSFImageOrigin(double origin[3]) {
  PointType thisOrigin;
  for (int i = 0; i < 3; i++)
    thisOrigin[i] = static_cast<PointValueType>(origin[i]);
  m_GibsonLanniPSFSource->SetOrigin(thisOrigin);
  m_GibsonLanniPSFSource->Modified();
}


void
DataModel
::GetPSFImageOrigin(double origin[3]) {
  PointType thisOrigin = m_GibsonLanniPSFSource->GetOrigin();
  for (int i = 0; i < 3; i++)
    origin[i] = static_cast<double>(thisOrigin[i]);
}


void
DataModel
::SetBSFImageOrigin(double origin[3]) {
  m_GibsonLanniBSFSource->SetOrigin(origin);
  m_GibsonLanniBSFSource->Modified();
}


void
DataModel
::GetBSFImageOrigin(double origin[3]) {
  PointType thisOrigin = m_GibsonLanniBSFSource->GetOrigin();
  for (int i = 0; i < 3; i++)
    origin[i] = static_cast<double>(thisOrigin[i]);
}


void
DataModel
::RecenterPSFImageOrigin() {
  double spacing[3];
  GetPSFImageVoxelSpacing(spacing);

  int size[3];
  GetPSFImageDimensions(size);

  PointType origin;
  for (int i = 0; i < 3; i++)
    origin[i] = -spacing[i]*static_cast<double>(size[i]-1)*0.5;
  m_GibsonLanniPSFSource->SetOrigin(origin);
}


void
DataModel
::SetPSFPointCenter(double center[3]) {
  PointType thisCenter;
  for (int i = 0; i < 3; i++)
    thisCenter[i] = static_cast<PointValueType>(center[i]);
  m_GibsonLanniPSFSource->SetPointCenter(thisCenter);
}


void
DataModel
::GetPSFPointCenter(double center[3]) {
  PointType thisCenter = m_GibsonLanniPSFSource->GetPointCenter();
  for (int i = 0; i < 3; i++)
    center[i] = static_cast<double>(thisCenter[i]);
}


void
DataModel
::SetBSFPointCenter(double center[3]) {
  PointType thisCenter;
  for (int i = 0; i < 3; i++)
    thisCenter[i] = static_cast<PointValueType>(center[i]);
  m_GibsonLanniBSFSource->SetBeadCenter(thisCenter);
}


void
DataModel
::GetBSFPointCenter(double center[3]) {
  PointType thisCenter = m_GibsonLanniBSFSource->GetBeadCenter();
  for (int i = 0; i < 3; i++)
    center[i] = static_cast<double>(thisCenter[i]);
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
::SetShearX(double shear) {
  m_GibsonLanniBSFSource->SetShearX(shear);
}


double
DataModel
::GetShearX() {
  return m_GibsonLanniBSFSource->GetShearX();
}


void
DataModel
::SetShearY(double shear) {
  m_GibsonLanniBSFSource->SetShearY(shear);
}


double
DataModel
::GetShearY() {
  return m_GibsonLanniBSFSource->GetShearY();
}


void
DataModel
::SetGLEmissionWavelength(double wavelength) {
  m_GibsonLanniPSFSource->SetEmissionWavelength(wavelength);
  m_GibsonLanniBSFSource->SetEmissionWavelength(wavelength);
}


double
DataModel
::GetGLEmissionWavelength() {
  return m_GibsonLanniBSFSource->GetEmissionWavelength();
}


void
DataModel
::SetGLNumericalAperture(double na) {
  m_GibsonLanniPSFSource->SetNumericalAperture(na);
  m_GibsonLanniBSFSource->SetNumericalAperture(na);
}


double
DataModel
::GetGLNumericalAperture() {
  return m_GibsonLanniBSFSource->GetNumericalAperture();
}


void
DataModel
::SetGLMagnification(double magnification) {
  m_GibsonLanniPSFSource->SetMagnification(magnification);
  m_GibsonLanniBSFSource->SetMagnification(magnification);
}


double
DataModel
::GetGLMagnification() {
  return m_GibsonLanniBSFSource->GetMagnification();
}


void
DataModel
::SetGLDesignCoverSlipRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetDesignCoverSlipRefractiveIndex(ri);
  m_GibsonLanniBSFSource->SetDesignCoverSlipRefractiveIndex(ri);
}


double
DataModel
::GetGLDesignCoverSlipRefractiveIndex() {
  return m_GibsonLanniBSFSource->GetDesignCoverSlipRefractiveIndex();
}


void
DataModel
::SetGLActualCoverSlipRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetActualCoverSlipRefractiveIndex(ri);
  m_GibsonLanniBSFSource->SetActualCoverSlipRefractiveIndex(ri);
}


double
DataModel
::GetGLActualCoverSlipRefractiveIndex() {
  return m_GibsonLanniBSFSource->GetActualCoverSlipRefractiveIndex();
}


void
DataModel
::SetGLDesignCoverSlipThickness(double thickness) {
  m_GibsonLanniPSFSource->SetDesignCoverSlipThickness(thickness);
  m_GibsonLanniBSFSource->SetDesignCoverSlipThickness(thickness);
}


double
DataModel
::GetGLDesignCoverSlipThickness() {
  return m_GibsonLanniBSFSource->GetDesignCoverSlipThickness();
}


void
DataModel
::SetGLActualCoverSlipThickness(double thickness) {
  m_GibsonLanniPSFSource->SetActualCoverSlipThickness(thickness);
  m_GibsonLanniBSFSource->SetActualCoverSlipThickness(thickness);
}


double
DataModel
::GetGLActualCoverSlipThickness() {
  return m_GibsonLanniBSFSource->GetActualCoverSlipThickness();
}


void
DataModel
::SetGLDesignImmersionOilRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetDesignImmersionOilRefractiveIndex(ri);
  m_GibsonLanniBSFSource->SetDesignImmersionOilRefractiveIndex(ri);
}


double
DataModel
::GetGLDesignImmersionOilRefractiveIndex() {
  return m_GibsonLanniBSFSource->GetDesignImmersionOilRefractiveIndex();
}


void
DataModel
::SetGLActualImmersionOilRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetActualImmersionOilRefractiveIndex(ri);
  m_GibsonLanniBSFSource->SetActualImmersionOilRefractiveIndex(ri);
}


double
DataModel
::GetGLActualImmersionOilRefractiveIndex() {
  return m_GibsonLanniBSFSource->GetActualImmersionOilRefractiveIndex();
}


void
DataModel
::SetGLDesignImmersionOilThickness(double thickness) {
  m_GibsonLanniPSFSource->SetDesignImmersionOilThickness(thickness);
  m_GibsonLanniBSFSource->SetDesignImmersionOilThickness(thickness);
}


double
DataModel
::GetGLDesignImmersionOilThickness() {
  return m_GibsonLanniBSFSource->GetDesignImmersionOilThickness();
}


void
DataModel
::SetGLDesignSpecimenLayerRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetDesignSpecimenLayerRefractiveIndex(ri);
  m_GibsonLanniBSFSource->SetDesignSpecimenLayerRefractiveIndex(ri);
}


double
DataModel
::GetGLDesignSpecimenLayerRefractiveIndex() {
  return m_GibsonLanniBSFSource->GetDesignSpecimenLayerRefractiveIndex();
}


void
DataModel
::SetGLActualSpecimenLayerRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetActualSpecimenLayerRefractiveIndex(ri);
  m_GibsonLanniBSFSource->SetActualSpecimenLayerRefractiveIndex(ri);
}


double
DataModel
::GetGLActualSpecimenLayerRefractiveIndex() {
  return m_GibsonLanniBSFSource->GetActualSpecimenLayerRefractiveIndex();
}


void
DataModel
::SetGLActualPointSourceDepthInSpecimenLayer(double depth) {
  m_GibsonLanniPSFSource->SetActualPointSourceDepthInSpecimenLayer(depth);
  m_GibsonLanniBSFSource->SetActualPointSourceDepthInSpecimenLayer(depth);
}


double
DataModel
::GetGLActualPointSourceDepthInSpecimenLayer() {
  return m_GibsonLanniBSFSource->GetActualPointSourceDepthInSpecimenLayer();
}


void
DataModel
::SetGLDesignDistanceFromBackFocalPlaneToDetector(double distance) {
  m_GibsonLanniPSFSource->SetDesignDistanceFromBackFocalPlaneToDetector(distance);
  m_GibsonLanniBSFSource->SetDesignDistanceFromBackFocalPlaneToDetector(distance);
}


double
DataModel
::GetGLDesignDistanceFromBackFocalPlaneToDetector() {
  return m_GibsonLanniBSFSource->GetDesignDistanceFromBackFocalPlaneToDetector();
}


void
DataModel
::SetGLActualDistanceFromBackFocalPlaneToDetector(double distance) {
  m_GibsonLanniPSFSource->SetActualDistanceFromBackFocalPlaneToDetector(distance);
  m_GibsonLanniBSFSource->SetActualDistanceFromBackFocalPlaneToDetector(distance);
}


double
DataModel
::GetGLActualDistanceFromBackFocalPlaneToDetector() {
  return m_GibsonLanniBSFSource->GetActualDistanceFromBackFocalPlaneToDetector();
}


void
DataModel
::SetGLBackgroundIntensity(double intensity) {
  m_GibsonLanniBSFSource->SetBackgroundIntensity(intensity);
}


double
DataModel
::GetGLBackgroundIntensity() {
  return m_GibsonLanniBSFSource->GetBackgroundIntensity();
}


void
DataModel
::SetGLMaximumIntensity(double intensity) {
  m_GibsonLanniBSFSource->SetMaximumIntensity(intensity);
}


double
DataModel
::GetGLMaximumIntensity() {
  return m_GibsonLanniBSFSource->GetMaximumIntensity();
}



void
DataModel
::SetGLParameterEnabled(unsigned int index, bool enabled) {
  try {
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
    ParametersMaskType* parametersMask = m_CostFunction->GetParametersMask();
    if (index < parametersMask->Size()) {
      enabled = parametersMask->GetElement(index) != 0;
    }
  } catch (...) {}

  return enabled;
}


void
DataModel
::SetZCoordinate(unsigned int index, double coordinate) {
  m_GibsonLanniBSFSource->SetZCoordinate(index, coordinate);
}

double
DataModel
::GetZCoordinate(unsigned int index) {
  return m_GibsonLanniBSFSource->GetZCoordinate(index);
}



void
DataModel
::SetUseCustomZCoordinates(bool use) {
  m_GibsonLanniBSFSource->SetUseCustomZCoordinates(use);
}


bool
DataModel
::GetUseCustomZCoordinates() {
  return m_GibsonLanniBSFSource->GetUseCustomZCoordinates();
}


double
DataModel
::GetImageComparisonMetricValue() {
  // Pass only the active parameter values to the cost function
  try {
    ParametersMaskType* mask = m_CostFunction->GetParametersMask();
    ParametersType activeParameters =
      ParametersType(m_CostFunction->GetNumberOfParameters());

    int activeIndex = 0;
    for (unsigned int i = 0; i < mask->Size(); i++) {
      if (mask->GetElement(i)) {
        activeParameters[activeIndex++] =
          m_GibsonLanniBSFSource->GetParameters()[i];
      }
    }

    return m_CostFunction->GetValue(activeParameters);
  } catch (...) {}

  return 0.0;
}


void
DataModel
::Optimize() {
  ParametersMaskType* mask = m_CostFunction->GetParametersMask();

  // Pluck out the active parameters
  ParametersType activeParameters(m_CostFunction->GetNumberOfParameters());
  ParametersType initialSimplexDelta(m_CostFunction->GetNumberOfParameters());
  int activeIndex = 0;
  for (unsigned int i = 0; i < mask->Size(); i++) {
    if (mask->GetElement(i)) {
      activeParameters[activeIndex] = m_GibsonLanniBSFSource->GetParameters()[i];
      initialSimplexDelta[activeIndex++] = m_InitialSimplexDelta[i];
    }
  }

  // Connect to the cost function, set the initial parameters, and optimize.
  m_ImageToImageCostFunction
    ->SetFixedImageRegion(m_GibsonLanniBSFSource->GetOutput()->GetLargestPossibleRegion());

  m_Optimizer = OptimizerType::New();
  m_Optimizer->AutomaticInitialSimplexOff();
  m_Optimizer->SetCostFunction(m_CostFunction);
  m_Optimizer->SetFunctionConvergenceTolerance(1e-1);
  m_Optimizer->SetInitialPosition(activeParameters);
  m_Optimizer->SetInitialSimplexDelta(initialSimplexDelta);

  m_Optimizer->StartOptimization();

  // Write the parameters back to the source object
  ParametersType optimizedParameters = m_Optimizer->GetCurrentPosition();
  ParametersType allParameters = m_GibsonLanniBSFSource->GetParameters();
  activeIndex = 0;
  for (unsigned int i = 0; i < mask->Size(); i++) {
    if (mask->GetElement(i)) {
      allParameters[i] = optimizedParameters[activeIndex++];
    }
  }

  m_GibsonLanniBSFSource->SetParameters(allParameters);
}

#endif // _DATA_MODEL_CXX_
