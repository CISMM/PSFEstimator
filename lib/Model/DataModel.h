#ifndef _DATA_MODEL_H_
#define _DATA_MODEL_H_

#include <string>

#include <Configuration.h>

#define ITK_MANUAL_INSTANTIATION
#include <itkGridImageSource.h>
#include <itkGibsonLanniBSFImageSource.h>
#include <itkGibsonLanniPSFImageSource.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkShiftScaleImageFilter.h>

#include <itkAmoebaOptimizer.h>
#include <itkMeanSquaresImageToImageMetric.h>
//#include <itkPoissonNoiseImageToImageMetric.h>
#include <itkImageToParametricImageSourceMetric.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <ITKImageToVTKImage.h>
#undef ITK_MANUAL_INSTANTIATION

#include <vtkAlgorithmOutput.h>


// This is the data model for the PSF Optimizer library.
class DataModel {

public:
  typedef float                                   FloatPixelType;
  static const unsigned int                       Dimension3 = 3;
  typedef itk::Image<FloatPixelType, Dimension3>  Float3DImageType;
  typedef itk::ImageFileReader<Float3DImageType>  Float3DImageReaderType;
  typedef itk::VTKImageExport<Float3DImageType>   Float3DExporterType;
  typedef Float3DImageType::PointType             Float3DPointType;

  typedef Float3DImageType::SpacingType      SpacingType;
  typedef Float3DImageType::SpacingValueType SpacingValueType;
  typedef Float3DImageType::SizeType         SizeType;
  typedef Float3DImageType::SizeValueType    SizeValueType;
  typedef Float3DImageType::PointType        PointType;
  typedef Float3DImageType::PointValueType   PointValueType;

  typedef itk::GridImageSource<Float3DImageType>
    DummyImageSourceType;
  typedef DummyImageSourceType::Pointer
    DummyImageSourcePointer;
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
  typedef itk::ImageToParametricImageSourceMetric<TImage, GibsonLanniBSFImageSourceType>
    ParametricCostFunctionType;
  typedef ParametricCostFunctionType::ParametersMaskType
    ParametersMaskType;
  typedef ParametricCostFunctionType::ParametersType
    ParametersType;

  typedef itk::NearestNeighborInterpolateImageFunction<TImage, double>
    InterpolatorType;
  //  typedef itk::PoissonNoiseImageToImageMetric<TImage, TImage>
  typedef itk::MeanSquaresImageToImageMetric<TImage, TImage>
    ImageToImageCostFunctionType;
  typedef itk::AmoebaOptimizer
    OptimizerType;

  DataModel();
  virtual ~DataModel();

  bool LoadSessionFile(const std::string& fileName);
  bool SaveSessionFile(const std::string& fileName);

  void SetInitialSimplexDeltas();

  void CreateImageFile(int xSize, int ySize, int zSize,
                       double xSpacing, double ySpacing, double zSpacing);
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

  int GetNumberOfProperties();

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

  // Sets the X and Y shear
  void SetShearX(double shear) {
    m_GibsonLanniBSFSource->SetShearX(shear);
  }
  double GetShearX() {
    return m_GibsonLanniBSFSource->GetShearX();
  }

  void SetShearY(double shear) {
    m_GibsonLanniBSFSource->SetShearY(shear);
  }
  double GetShearY() {
    return m_GibsonLanniBSFSource->GetShearY();
  }

  void UpdateGibsonLanniPSFImage();
  void UpdateGibsonLanniBSFImage();

  void SetBeadRadius(double radius);
  double GetBeadRadius();

  void  SetGLEmissionWavelength(double wavelength) {
    m_GibsonLanniPSFSource->SetEmissionWavelength(wavelength);
    m_GibsonLanniBSFSource->SetEmissionWavelength(wavelength);
  }
  double GetGLEmissionWavelength() {
    return m_GibsonLanniBSFSource->GetEmissionWavelength(); }

  void  SetGLNumericalAperture(double na) {
    m_GibsonLanniPSFSource->SetNumericalAperture(na);
    m_GibsonLanniBSFSource->SetNumericalAperture(na);
  }
  double GetGLNumericalAperture() {
    return m_GibsonLanniBSFSource->GetNumericalAperture(); }

  void  SetGLMagnification(double magnification) {
    m_GibsonLanniPSFSource->SetMagnification(magnification);
    m_GibsonLanniBSFSource->SetMagnification(magnification);
  }
  double GetGLMagnification() {
    return m_GibsonLanniBSFSource->GetMagnification(); }

  void  SetGLDesignCoverSlipRefractiveIndex(double ri) {
    m_GibsonLanniPSFSource->SetDesignCoverSlipRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetDesignCoverSlipRefractiveIndex(ri);
  }
  double GetGLDesignCoverSlipRefractiveIndex() {
    return m_GibsonLanniBSFSource->GetDesignCoverSlipRefractiveIndex(); }
  void  SetGLActualCoverSlipRefractiveIndex(double ri) {
    m_GibsonLanniPSFSource->SetActualCoverSlipRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetActualCoverSlipRefractiveIndex(ri);
  }
  double GetGLActualCoverSlipRefractiveIndex() {
    return m_GibsonLanniBSFSource->GetActualCoverSlipRefractiveIndex(); }

  void  SetGLDesignCoverSlipThickness(double thickness) {
    m_GibsonLanniPSFSource->SetDesignCoverSlipThickness(thickness);
    m_GibsonLanniBSFSource->SetDesignCoverSlipThickness(thickness);
  }
  double GetGLDesignCoverSlipThickness() {
    return m_GibsonLanniBSFSource->GetDesignCoverSlipThickness(); }
  void  SetGLActualCoverSlipThickness(double thickness) {
    m_GibsonLanniPSFSource->SetActualCoverSlipThickness(thickness);
    m_GibsonLanniBSFSource->SetActualCoverSlipThickness(thickness);
  }
  double GetGLActualCoverSlipThickness() {
    return m_GibsonLanniBSFSource->GetActualCoverSlipThickness(); }

  void  SetGLDesignImmersionOilRefractiveIndex(double ri) {
    m_GibsonLanniPSFSource->SetDesignImmersionOilRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetDesignImmersionOilRefractiveIndex(ri);
  }
  double GetGLDesignImmersionOilRefractiveIndex() {
    return m_GibsonLanniBSFSource->GetDesignImmersionOilRefractiveIndex(); }
  void  SetGLActualImmersionOilRefractiveIndex(double ri) {
    m_GibsonLanniPSFSource->SetActualImmersionOilRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetActualImmersionOilRefractiveIndex(ri);
  }
  double GetGLActualImmersionOilRefractiveIndex() {
    return m_GibsonLanniBSFSource->GetActualImmersionOilRefractiveIndex(); }

  void  SetGLDesignImmersionOilThickness(double thickness) {
    m_GibsonLanniPSFSource->SetDesignImmersionOilThickness(thickness);
    m_GibsonLanniBSFSource->SetDesignImmersionOilThickness(thickness);
  }
  double GetGLDesignImmersionOilThickness() {
    return m_GibsonLanniBSFSource->GetDesignImmersionOilThickness(); }

  void  SetGLDesignSpecimenLayerRefractiveIndex(double ri) {
    m_GibsonLanniPSFSource->SetDesignSpecimenLayerRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetDesignSpecimenLayerRefractiveIndex(ri);
  }
  double GetGLDesignSpecimenLayerRefractiveIndex() {
    return m_GibsonLanniBSFSource->GetDesignSpecimenLayerRefractiveIndex(); }
  void  SetGLActualSpecimenLayerRefractiveIndex(double ri) {
    m_GibsonLanniPSFSource->SetActualSpecimenLayerRefractiveIndex(ri);
    m_GibsonLanniBSFSource->SetActualSpecimenLayerRefractiveIndex(ri);
  }
  double GetGLActualSpecimenLayerRefractiveIndex() {
    return m_GibsonLanniBSFSource->GetActualSpecimenLayerRefractiveIndex(); }

  void  SetGLActualPointSourceDepthInSpecimenLayer(double depth) {
    m_GibsonLanniPSFSource->SetActualPointSourceDepthInSpecimenLayer(depth);
    m_GibsonLanniBSFSource->SetActualPointSourceDepthInSpecimenLayer(depth);
  }
  double GetGLActualPointSourceDepthInSpecimenLayer() {
    return m_GibsonLanniBSFSource->GetActualPointSourceDepthInSpecimenLayer(); }

  void  SetGLDesignDistanceFromBackFocalPlaneToDetector(double distance) {
    m_GibsonLanniPSFSource->SetDesignDistanceFromBackFocalPlaneToDetector(distance);
    m_GibsonLanniBSFSource->SetDesignDistanceFromBackFocalPlaneToDetector(distance);
  }
  double GetGLDesignDistanceFromBackFocalPlaneToDetector() {
    return m_GibsonLanniBSFSource->GetDesignDistanceFromBackFocalPlaneToDetector(); }
  void  SetGLActualDistanceFromBackFocalPlaneToDetector(double distance) {
    m_GibsonLanniPSFSource->SetActualDistanceFromBackFocalPlaneToDetector(distance);
    m_GibsonLanniBSFSource->SetActualDistanceFromBackFocalPlaneToDetector(distance);
  }
  double GetGLActualDistanceFromBackFocalPlaneToDetector() {
    return m_GibsonLanniBSFSource->GetActualDistanceFromBackFocalPlaneToDetector(); }

  void  SetGLBackgroundIntensity(double intensity) {
    m_GibsonLanniBSFSource->SetBackgroundIntensity(intensity); }
  double GetGLBackgroundIntensity() {
    return m_GibsonLanniBSFSource->GetBackgroundIntensity(); }

  void  SetGLMaximumIntensity(double intensity) {
    m_GibsonLanniBSFSource->SetMaximumIntensity(intensity); }
  double GetGLMaximumIntensity() {
    return m_GibsonLanniBSFSource->GetMaximumIntensity(); }

  void SetGLParameterEnabled(unsigned int index, bool enabled);
  bool GetGLParameterEnabled(unsigned int index);

  void SetZCoordinate(unsigned int index, double coordinate) {
    m_GibsonLanniBSFSource->SetZCoordinate(index, coordinate);
  }
  double GetZCoordinate(unsigned int index) {
    return m_GibsonLanniBSFSource->GetZCoordinate(index);
  }

  void SetUseCustomZCoordinates(bool use) {
    m_GibsonLanniBSFSource->SetUseCustomZCoordinates(use);
  }
  bool GetUseCustomZCoordinates() {
    return m_GibsonLanniBSFSource->GetUseCustomZCoordinates();
  }

  double GetImageComparisonMetricValue();

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
  ParametricCostFunctionType::Pointer m_CostFunction;

  // The delegate cost function used by m_CostFunction
  ImageToImageCostFunctionType::Pointer  m_ImageToImageCostFunction;

  // The optimizer
  OptimizerType::Pointer m_Optimizer;

  // Optimizer initial simplex deltas. These specify the starting size of the
  // search simplex in the Nelder-Mead (Amoeba) optimization algorithm. It is
  // important to set these appropriately because small changes to some
  // parameters (e.g. actual refractive index) produce huge changes in the
  // shape of the PSF.
  ParametersType m_InitialSimplexDelta;
};

// _DATA_MODEL_H_
#endif
