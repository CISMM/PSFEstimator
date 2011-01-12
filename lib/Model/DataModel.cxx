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

#include <itkBeadSpreadFunctionImageSource.txx>
#include <itkGibsonLanniPSFImageSource.txx>
#include <itkGridImageSource.txx>
#include <itkImageFileReader.txx>
#include <itkImageFileWriter.txx>
#include <itkImageToParametricImageSourceMetric.txx>
#include <itkMeanSquaresImageToImageMetric.txx>
#include <itkMinimumMaximumImageCalculator.txx>
//#include <itkPoissonNoiseImageToImageMetric.cxx>
#include <itkShiftScaleImageFilter.txx>
#include <itkSubtractImageFilter.h>

#include <ITKImageToVTKImage.cxx>

#include <vtkAlgorithm.h>

#include <DataModel.h>


DataModel
::DataModel() {
  // ITK will detect the number of cores on the system and set the
  // global number of threads to the number of cores by default.
  // Here we can override that setting if the proper environment
  // variable is set.
  char *var = getenv("PSFEstimator_THREADS");
  if (var) {
    int numberOfThreads = atoi(var);
    if (numberOfThreads > 0)
      SetNumberOfThreads(numberOfThreads);
  }


  m_MeasuredImageData = NULL;

  m_GibsonLanniPSFSource             = GibsonLanniPSFImageSourceType::New();
  m_BeadSpreadFunctionSource         = BeadSpreadFunctionImageSourceType::New();
  m_BeadSpreadFunctionSource->SetKernelSource(GibsonLanniPSFImageSourceType::New());
  m_BSFDifferenceImageFilter         = DifferenceFilterType::New();

  m_MeasuredImageMinMaxFilter        = MinMaxType::New();
  m_PSFImageMinMaxFilter             = MinMaxType::New();
  m_BSFImageMinMaxFilter             = MinMaxType::New();
  m_BSFDifferenceImageMinMaxFilter   = MinMaxType::New();

  m_MeasuredImageITKToVTKFilter      = new ITKImageToVTKImage<TImage>();
  m_PSFImageITKToVTKFilter           = new ITKImageToVTKImage<TImage>();
  m_BSFImageITKToVTKFilter           = new ITKImageToVTKImage<TImage>();
  m_BSFDifferenceImageITKToVTKFilter = new ITKImageToVTKImage<TImage>();

  m_ImageToImageCostFunction = ImageToImageCostFunctionType::New();
  m_CostFunction = ParametricCostFunctionType::New();
  m_CostFunction->SetInterpolator(InterpolatorType::New());
  m_CostFunction->SetDelegateMetric(m_ImageToImageCostFunction);
  m_CostFunction->SetMovingImageSource(m_BeadSpreadFunctionSource);

  SetInitialSimplexDeltas();
}


DataModel
::~DataModel() {
  delete m_MeasuredImageITKToVTKFilter;
  delete m_PSFImageITKToVTKFilter;
  delete m_BSFImageITKToVTKFilter;
  delete m_BSFDifferenceImageITKToVTKFilter;
}


bool
DataModel
::LoadSessionFile(const std::string& fileName) {
  Configuration config;
  config.Parse(fileName);

  // Read the settings from the configuration structure
  bool success = LoadImageFile(config.GetValue("FileInfo", "FileName"));
  if (!success) {
    return false;
  }

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
    SetSize(m_BeadSpreadFunctionSource->GetNumberOfParameters());

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
  m_BeadSpreadFunctionSource->SetOrigin(origin);
  m_MeasuredImageData->SetOrigin(origin);

  m_PSFImageMinMaxFilter = MinMaxType::New();
  m_PSFImageMinMaxFilter->SetImage(m_GibsonLanniPSFSource->GetOutput());
  m_PSFImageITKToVTKFilter->SetInput(m_GibsonLanniPSFSource->GetOutput());

  m_BSFImageMinMaxFilter = MinMaxType::New();
  m_BSFImageMinMaxFilter->SetImage(m_BeadSpreadFunctionSource->GetOutput());
  m_BSFImageITKToVTKFilter->SetInput(m_BeadSpreadFunctionSource->GetOutput());

  // Set up cost function
  m_CostFunction->SetFixedImage(m_MeasuredImageData);
  m_CostFunction->SetMovingImageSource(m_BeadSpreadFunctionSource);
}


bool
DataModel
::LoadImageFile(std::string fileName) {
  m_ImageFileName = fileName;
  ScalarFileReaderType::Pointer reader = ScalarFileReaderType::New();
  reader->SetFileName(fileName.c_str());

  try {
    reader->Update();
  } catch (...) {
    return false;
  }
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
  m_BeadSpreadFunctionSource->SetOrigin(origin);
  m_MeasuredImageData->SetOrigin(origin);

  m_PSFImageMinMaxFilter = MinMaxType::New();
  m_PSFImageMinMaxFilter->SetImage(m_GibsonLanniPSFSource->GetOutput());
  m_PSFImageITKToVTKFilter->SetInput(m_GibsonLanniPSFSource->GetOutput());

  m_BSFImageMinMaxFilter = MinMaxType::New();
  m_BSFImageMinMaxFilter->SetImage(m_BeadSpreadFunctionSource->GetOutput());
  m_BSFImageITKToVTKFilter->SetInput(m_BeadSpreadFunctionSource->GetOutput());

  m_BSFDifferenceImageFilter->SetInput1(m_MeasuredImageData);
  m_BSFDifferenceImageFilter->SetInput2(m_BeadSpreadFunctionSource->GetOutput());

  m_BSFDifferenceImageMinMaxFilter = MinMaxType::New();
  m_BSFDifferenceImageMinMaxFilter->SetImage(m_BSFDifferenceImageFilter->GetOutput());
  m_BSFDifferenceImageITKToVTKFilter->SetInput(m_BSFDifferenceImageFilter->GetOutput());

  // Set up cost function
  m_CostFunction->SetFixedImage(m_MeasuredImageData);

  return true;
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
  scaler->SetInput(m_BeadSpreadFunctionSource->GetOutput());
  //double min = GetBSFImageDataMinimum();
  //double max = GetBSFImageDataMaximum();
  //scaler->SetShift(-min);
  //scaler->SetScale(65535.0f / (max - min));
  scaler->SetShift(0.0);
  scaler->SetScale(1.0);

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
  m_BeadSpreadFunctionSource->SetShearX(shearX);
  double shearY = c.GetValueAsDouble(sec, "ShearY");
  m_BeadSpreadFunctionSource->SetShearY(shearY);

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
  SetGLIntensityShift
    (c.GetValueAsDouble(sec, "IntensityShift",
                       GetGLIntensityShift()));
  SetGLIntensityScale
    (c.GetValueAsDouble(sec, "IntensityScale",
                       GetGLIntensityScale()));

  sec = std::string("ZSliceCoordinates");

  SetUseCustomZCoordinates(c.GetValueAsBool(sec, "UseCustomZCoordinates",
                                            GetUseCustomZCoordinates()));

  for (unsigned int i = 0; i < m_BeadSpreadFunctionSource->GetSize()[2]; i++) {
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
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "IntensityShift"));
  SetGLParameterEnabled(index++, c.GetValueAsBool(sec, "IntensityScale"));
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

  c.SetValueFromDouble(sec, "ShearX", m_BeadSpreadFunctionSource->GetShearX());
  c.SetValueFromDouble(sec, "ShearY", m_BeadSpreadFunctionSource->GetShearY());

  c.SetValueFromDouble(sec, "EmissionWavelength",
                       m_BeadSpreadFunctionSource->GetParameter(11));

  c.SetValueFromDouble(sec, "NumericalAperture",
                       m_BeadSpreadFunctionSource->GetParameter(12));
  c.SetValueFromDouble(sec, "Magnification",
                       m_BeadSpreadFunctionSource->GetParameter(13));
  c.SetValueFromDouble(sec, "DesignCoverSlipRefractiveIndex",
                       m_BeadSpreadFunctionSource->GetParameter(14));
  c.SetValueFromDouble(sec, "ActualCoverSlipRefractiveIndex",
                       m_BeadSpreadFunctionSource->GetParameter(15));
  c.SetValueFromDouble(sec, "DesignCoverSlipThickness",
		       m_BeadSpreadFunctionSource->GetParameter(16));
  c.SetValueFromDouble(sec, "ActualCoverSlipThickness",
                       m_BeadSpreadFunctionSource->GetParameter(17));
  c.SetValueFromDouble(sec, "DesignImmersionOilRefractiveIndex",
		       m_BeadSpreadFunctionSource->GetParameter(18));
  c.SetValueFromDouble(sec, "ActualImmersionOilRefractiveIndex",
                       m_BeadSpreadFunctionSource->GetParameter(19));
  c.SetValueFromDouble(sec, "DesignImmersionOilThickness",
                       m_BeadSpreadFunctionSource->GetParameter(20));
  c.SetValueFromDouble(sec, "DesignSpecimenLayerRefractiveIndex",
                       m_BeadSpreadFunctionSource->GetParameter(21));
  c.SetValueFromDouble(sec, "ActualSpecimenLayerRefractiveIndex",
                       m_BeadSpreadFunctionSource->GetParameter(22));
  c.SetValueFromDouble(sec, "ActualPointSourceDepthInSpecimenLayer",
                       m_BeadSpreadFunctionSource->GetParameter(23));
  c.SetValueFromDouble(sec, "DesignDistanceFromBackFocalPlaneToDetector",
                       m_BeadSpreadFunctionSource->GetParameter(24));
  c.SetValueFromDouble(sec, "ActualDistanceFromBackFocalPlaneToDetector",
                       m_BeadSpreadFunctionSource->GetParameter(25));
  c.SetValueFromDouble(sec, "IntensityShift",
                       m_BeadSpreadFunctionSource->GetParameter(9));
  c.SetValueFromDouble(sec, "IntensityScale",
                       m_BeadSpreadFunctionSource->GetParameter(10));

  sec = std::string("ZSliceCoordinates");

  c.SetValueFromBool(sec, "UseCustomZCoordinates", GetUseCustomZCoordinates());

  for (unsigned int i = 0; i < m_BeadSpreadFunctionSource->GetSize()[2]; i++) {
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
  c.SetValueFromBool(sec, "IntensityShift", GetGLParameterEnabled(index++));
  c.SetValueFromBool(sec, "IntensityScale", GetGLParameterEnabled(index++));
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


vtkAlgorithmOutput*
DataModel
::GetBSFDifferenceImageOutputPort() {
  return m_BSFDifferenceImageITKToVTKFilter->GetOutputPort();
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
  return m_BeadSpreadFunctionSource->GetNumberOfParameters();
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

  m_BeadSpreadFunctionSource->UpdateLargestPossibleRegion();

  m_BSFImageMinMaxFilter->Compute();
  return m_BSFImageMinMaxFilter->GetMinimum();
}


double
DataModel
::GetBSFImageDataMaximum() {
  if (!GetMeasuredImageData()) {
    return 0.0;
  }

  m_BeadSpreadFunctionSource->UpdateLargestPossibleRegion();

  m_BSFImageMinMaxFilter->Compute();
  return m_BSFImageMinMaxFilter->GetMaximum();
}


double
DataModel
::GetBSFDifferenceImageDataMinimum() {
  if (!GetMeasuredImageData()) {
    return 0.0;
  }

  m_BSFDifferenceImageFilter->Modified();
  m_BSFDifferenceImageFilter->UpdateLargestPossibleRegion();

  m_BSFDifferenceImageMinMaxFilter->Compute();
  return m_BSFDifferenceImageMinMaxFilter->GetMinimum();
}


double
DataModel
::GetBSFDifferenceImageDataMaximum() {
  if (!GetMeasuredImageData()) {
    return 0.0;
  }

  m_BSFDifferenceImageFilter->Modified();
  m_BSFDifferenceImageFilter->UpdateLargestPossibleRegion();

  m_BSFDifferenceImageMinMaxFilter->Compute();
  return m_BSFDifferenceImageMinMaxFilter->GetMaximum();
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

  m_BeadSpreadFunctionSource->SetSize(size);
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

  m_BeadSpreadFunctionSource->SetSpacing(thisSpacing);
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
  m_BeadSpreadFunctionSource->SetOrigin(origin);
  m_BeadSpreadFunctionSource->Modified();
}


void
DataModel
::GetBSFImageOrigin(double origin[3]) {
  PointType thisOrigin = m_BeadSpreadFunctionSource->GetOrigin();
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
  m_BeadSpreadFunctionSource->SetBeadCenter(thisCenter);
}


void
DataModel
::GetBSFPointCenter(double center[3]) {
  PointType thisCenter = m_BeadSpreadFunctionSource->GetBeadCenter();
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
  m_BeadSpreadFunctionSource->UpdateLargestPossibleRegion();
  m_BSFImageMinMaxFilter->Compute();
}


void
DataModel
::UpdateBSFDifferenceImage() {
  m_BSFDifferenceImageFilter->UpdateLargestPossibleRegion();
  m_BSFDifferenceImageMinMaxFilter->Compute();
}


void
DataModel
::SetBeadRadius(double radius) {
  m_BeadSpreadFunctionSource->SetBeadRadius(radius);
}


double
DataModel
::GetBeadRadius() {
  return m_BeadSpreadFunctionSource->GetBeadRadius();
}


void
DataModel
::SetShearX(double shear) {
  m_BeadSpreadFunctionSource->SetShearX(shear);
}


double
DataModel
::GetShearX() {
  return m_BeadSpreadFunctionSource->GetShearX();
}


void
DataModel
::SetShearY(double shear) {
  m_BeadSpreadFunctionSource->SetShearY(shear);
}


double
DataModel
::GetShearY() {
  return m_BeadSpreadFunctionSource->GetShearY();
}


void
DataModel
::SetGLEmissionWavelength(double wavelength) {
  m_GibsonLanniPSFSource->SetEmissionWavelength(wavelength);
  m_BeadSpreadFunctionSource->SetParameter(11, wavelength);
}


double
DataModel
::GetGLEmissionWavelength() {
  return m_BeadSpreadFunctionSource->GetParameter(11);
}


void
DataModel
::SetGLNumericalAperture(double na) {
  m_GibsonLanniPSFSource->SetNumericalAperture(na);
  m_BeadSpreadFunctionSource->SetParameter(12, na);
}


double
DataModel
::GetGLNumericalAperture() {
  return m_BeadSpreadFunctionSource->GetParameter(12);
}


void
DataModel
::SetGLMagnification(double magnification) {
  m_GibsonLanniPSFSource->SetMagnification(magnification);
  m_BeadSpreadFunctionSource->SetParameter(13, magnification);
}


double
DataModel
::GetGLMagnification() {
  return m_BeadSpreadFunctionSource->GetParameter(13);
}


void
DataModel
::SetGLDesignCoverSlipRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetDesignCoverSlipRefractiveIndex(ri);
  m_BeadSpreadFunctionSource->SetParameter(14, ri);
}


double
DataModel
::GetGLDesignCoverSlipRefractiveIndex() {
  return m_BeadSpreadFunctionSource->GetParameter(14);
}


void
DataModel
::SetGLActualCoverSlipRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetActualCoverSlipRefractiveIndex(ri);
  m_BeadSpreadFunctionSource->SetParameter(15, ri);
}


double
DataModel
::GetGLActualCoverSlipRefractiveIndex() {
  return m_BeadSpreadFunctionSource->GetParameter(15);
}


void
DataModel
::SetGLDesignCoverSlipThickness(double thickness) {
  m_GibsonLanniPSFSource->SetDesignCoverSlipThickness(thickness);
  m_BeadSpreadFunctionSource->SetParameter(16, thickness);
}


double
DataModel
::GetGLDesignCoverSlipThickness() {
  return m_BeadSpreadFunctionSource->GetParameter(16);
}


void
DataModel
::SetGLActualCoverSlipThickness(double thickness) {
  m_GibsonLanniPSFSource->SetActualCoverSlipThickness(thickness);
  m_BeadSpreadFunctionSource->SetParameter(17, thickness);
}


double
DataModel
::GetGLActualCoverSlipThickness() {
  return m_BeadSpreadFunctionSource->GetParameter(17);
}


void
DataModel
::SetGLDesignImmersionOilRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetDesignImmersionOilRefractiveIndex(ri);
  m_BeadSpreadFunctionSource->SetParameter(18, ri);
}


double
DataModel
::GetGLDesignImmersionOilRefractiveIndex() {
  return m_BeadSpreadFunctionSource->GetParameter(18);
}


void
DataModel
::SetGLActualImmersionOilRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetActualImmersionOilRefractiveIndex(ri);
  m_BeadSpreadFunctionSource->SetParameter(19, ri);
}


double
DataModel
::GetGLActualImmersionOilRefractiveIndex() {
  return m_BeadSpreadFunctionSource->GetParameter(19);
}


void
DataModel
::SetGLDesignImmersionOilThickness(double thickness) {
  m_GibsonLanniPSFSource->SetDesignImmersionOilThickness(thickness);
  m_BeadSpreadFunctionSource->SetParameter(20, thickness);
}


double
DataModel
::GetGLDesignImmersionOilThickness() {
  return m_BeadSpreadFunctionSource->GetParameter(20);
}


void
DataModel
::SetGLDesignSpecimenLayerRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetDesignSpecimenLayerRefractiveIndex(ri);
  m_BeadSpreadFunctionSource->SetParameter(21, ri);
}


double
DataModel
::GetGLDesignSpecimenLayerRefractiveIndex() {
  return m_BeadSpreadFunctionSource->GetParameter(21);
}


void
DataModel
::SetGLActualSpecimenLayerRefractiveIndex(double ri) {
  m_GibsonLanniPSFSource->SetActualSpecimenLayerRefractiveIndex(ri);
  m_BeadSpreadFunctionSource->SetParameter(22, ri);
}


double
DataModel
::GetGLActualSpecimenLayerRefractiveIndex() {
  return m_BeadSpreadFunctionSource->GetParameter(22);
}


void
DataModel
::SetGLActualPointSourceDepthInSpecimenLayer(double depth) {
  m_GibsonLanniPSFSource->SetActualPointSourceDepthInSpecimenLayer(depth);
  m_BeadSpreadFunctionSource->SetParameter(23, depth);
}


double
DataModel
::GetGLActualPointSourceDepthInSpecimenLayer() {
  return m_BeadSpreadFunctionSource->GetParameter(23);
}


void
DataModel
::SetGLDesignDistanceFromBackFocalPlaneToDetector(double distance) {
  m_GibsonLanniPSFSource->SetDesignDistanceFromBackFocalPlaneToDetector(distance);
  m_BeadSpreadFunctionSource->SetParameter(24, distance);
}


double
DataModel
::GetGLDesignDistanceFromBackFocalPlaneToDetector() {
  return m_BeadSpreadFunctionSource->GetParameter(24);
}


void
DataModel
::SetGLActualDistanceFromBackFocalPlaneToDetector(double distance) {
  m_GibsonLanniPSFSource->SetActualDistanceFromBackFocalPlaneToDetector(distance);
  m_BeadSpreadFunctionSource->SetParameter(25, distance);
}


double
DataModel
::GetGLActualDistanceFromBackFocalPlaneToDetector() {
  return m_BeadSpreadFunctionSource->GetParameter(25);
}


void
DataModel
::SetGLIntensityShift(double intensity) {
  m_BeadSpreadFunctionSource->SetParameter(10, intensity);
}


double
DataModel
::GetGLIntensityShift() {
  return m_BeadSpreadFunctionSource->GetParameter(10);
}


void
DataModel
::SetGLIntensityScale(double scale) {
  m_BeadSpreadFunctionSource->SetParameter(9, scale);
}


double
DataModel
::GetGLIntensityScale() {
  return m_BeadSpreadFunctionSource->GetParameter(9);
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
  m_BeadSpreadFunctionSource->SetZCoordinate(index, coordinate);
}

double
DataModel
::GetZCoordinate(unsigned int index) {
  return m_BeadSpreadFunctionSource->GetZCoordinate(index);
}



void
DataModel
::SetUseCustomZCoordinates(bool use) {
  m_BeadSpreadFunctionSource->SetUseCustomZCoordinates(use);
}


bool
DataModel
::GetUseCustomZCoordinates() {
  return m_BeadSpreadFunctionSource->GetUseCustomZCoordinates();
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
          m_BeadSpreadFunctionSource->GetParameters()[i];
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
      activeParameters[activeIndex] = m_BeadSpreadFunctionSource->GetParameters()[i];
      initialSimplexDelta[activeIndex++] = m_InitialSimplexDelta[i];
    }
  }

  // Connect to the cost function, set the initial parameters, and optimize.
  m_ImageToImageCostFunction
    ->SetFixedImageRegion(m_BeadSpreadFunctionSource->GetOutput()->GetLargestPossibleRegion());

  m_Optimizer = OptimizerType::New();
  m_Optimizer->AutomaticInitialSimplexOff();
  m_Optimizer->SetCostFunction(m_CostFunction);
  m_Optimizer->SetFunctionConvergenceTolerance(1e-1);
  m_Optimizer->SetInitialPosition(activeParameters);
  m_Optimizer->SetInitialSimplexDelta(initialSimplexDelta);

  m_Optimizer->StartOptimization();

  // Write the parameters back to the source object
  ParametersType optimizedParameters = m_Optimizer->GetCurrentPosition();
  ParametersType allParameters = m_BeadSpreadFunctionSource->GetParameters();
  activeIndex = 0;
  for (unsigned int i = 0; i < mask->Size(); i++) {
    if (mask->GetElement(i)) {
      allParameters[i] = optimizedParameters[activeIndex++];
    }
  }

  m_BeadSpreadFunctionSource->SetParameters(allParameters);
}

#endif // _DATA_MODEL_CXX_
