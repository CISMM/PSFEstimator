#ifndef _DATA_MODEL_H_
#define _DATA_MODEL_H_

#include <string>

#include <Configuration.h>

#define ITK_MANUAL_INSTANTIATION
#include <itkBeadSpreadFunctionImageSource.h>
#include <itkGridImageSource.h>
#include <itkGibsonLanniPSFImageSource.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkShiftScaleImageFilter.h>
#include <itkSubtractImageFilter.h>

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
  typedef itk::BeadSpreadFunctionImageSource< Float3DImageType >
    BeadSpreadFunctionImageSourceType;
  typedef BeadSpreadFunctionImageSourceType::Pointer
    BeadSpreadFunctionImageSourcePointer;

  typedef itk::SubtractImageFilter<Float3DImageType, Float3DImageType, Float3DImageType>
    DifferenceFilterType;
  typedef DifferenceFilterType::Pointer
    DifferenceFilterPointer;

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
  typedef itk::ImageToParametricImageSourceMetric<TImage, BeadSpreadFunctionImageSourceType>
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
  bool LoadImageFile(std::string fileName);
  void SavePSFImageFile(std::string fileName);
  void SaveBSFImageFile(std::string fileName);

  void SetConfiguration(Configuration & configuration);
  void GetConfiguration(Configuration & configuration);

  std::string GetMeasuredImageFileName();

  // Number of threads to run for all the multithreaded ITK algorithms
  void SetNumberOfThreads(int threads);
  int  GetNumberOfThreads();

  void            SetMeasuredImageData(TImage::Pointer image);
  TImage::Pointer GetMeasuredImageData();

  // Returns the VTK output port for the original scalar image data.
  vtkAlgorithmOutput* GetMeasuredImageOutputPort();

  // Returns the VTK output port for the generated point-spread function
  // image data.
  vtkAlgorithmOutput* GetPSFImageOutputPort();

  // Returns the VTK output port for the generated bead-spread function
  // image data.
  vtkAlgorithmOutput* GetBSFImageOutputPort();

  // Returns the VTK output port for the difference between the
  // measured and calculated bead-spread function
  vtkAlgorithmOutput* GetBSFDifferenceImageOutputPort();

  double GetMeasuredImageDataMinimum();
  double GetMeasuredImageDataMaximum();
  Float3DPointType GetMeasuredImageDataMaximumCoordinates();

  void   GetMeasuredImageDimensions(int dimensions[3]);

  int    GetNumberOfProperties();

  void   SetMeasuredImageVoxelSpacing(double spacing[3]);
  void   SetMeasuredImageVoxelSpacing(int dimension, double spacing);
  void   GetMeasuredImageVoxelSpacing(double spacing[3]);

  void   SetMeasuredImageOrigin(double origin[3]);
  void   GetMeasuredImageOrigin(double origin[3]);

  double GetPSFImageDataMinimum();
  double GetPSFImageDataMaximum();

  double GetBSFImageDataMinimum();
  double GetBSFImageDataMaximum();

  double GetBSFDifferenceImageDataMinimum();
  double GetBSFDifferenceImageDataMaximum();

  void   SetPSFImageDimensions(int dimensions[3]);
  void   SetPSFImageDimension(int index, int dimension);
  void   GetPSFImageDimensions(int dimensions[3]);

  void   SetBSFImageDimensions(int dimensions[3]);
  void   SetBSFImageDimension(int index, int dimension);
  void   GetBSFImageDimensions(int dimensions[3]);

  void   SetPSFImageVoxelSpacing(double spacing[3]);
  void   SetPSFImageVoxelSpacing(int dimension, double spacing);
  void   GetPSFImageVoxelSpacing(double spacing[3]);

  void   SetBSFImageVoxelSpacing(double spacing[3]);
  void   SetBSFImageVoxelSpacing(int dimension, double spacing);
  void   GetBSFImageVoxelSpacing(double spacing[3]);

  void   SetPSFImageOrigin(double origin[3]);
  void   GetPSFImageOrigin(double origin[3]);

  void   SetBSFImageOrigin(double origin[3]);
  void   GetBSFImageOrigin(double origin[3]);

  // Recenters the PSF image origin to the center of the image bounds.
  void   RecenterPSFImageOrigin();

  // Sets the PSF center
  void   SetPSFPointCenter(double center[3]);
  void   GetPSFPointCenter(double center[3]);
  void   SetBSFPointCenter(double center[3]);
  void   GetBSFPointCenter(double center[3]);

  // Sets the X and Y shear
  void   SetShearX(double shear);
  double GetShearX();

  void   SetShearY(double shear);
  double GetShearY();

  void   UpdateGibsonLanniPSFImage();
  void   UpdateGibsonLanniBSFImage();
  void   UpdateBSFDifferenceImage();

  void   SetBeadRadius(double radius);
  double GetBeadRadius();

  void   SetGLEmissionWavelength(double wavelength);
  double GetGLEmissionWavelength();

  void   SetGLNumericalAperture(double na);
  double GetGLNumericalAperture();

  void   SetGLMagnification(double magnification);
  double GetGLMagnification();

  void   SetGLDesignCoverSlipRefractiveIndex(double ri);
  double GetGLDesignCoverSlipRefractiveIndex();

  void   SetGLActualCoverSlipRefractiveIndex(double ri);
  double GetGLActualCoverSlipRefractiveIndex();

  void   SetGLDesignCoverSlipThickness(double thickness);
  double GetGLDesignCoverSlipThickness();

  void   SetGLActualCoverSlipThickness(double thickness);
  double GetGLActualCoverSlipThickness();

  void   SetGLDesignImmersionOilRefractiveIndex(double ri);
  double GetGLDesignImmersionOilRefractiveIndex();

  void   SetGLActualImmersionOilRefractiveIndex(double ri);
  double GetGLActualImmersionOilRefractiveIndex();

  void   SetGLDesignImmersionOilThickness(double thickness);
  double GetGLDesignImmersionOilThickness();

  void   SetGLDesignSpecimenLayerRefractiveIndex(double ri);
  double GetGLDesignSpecimenLayerRefractiveIndex();
  void   SetGLActualSpecimenLayerRefractiveIndex(double ri);
  double GetGLActualSpecimenLayerRefractiveIndex();

  void   SetGLActualPointSourceDepthInSpecimenLayer(double depth);
  double GetGLActualPointSourceDepthInSpecimenLayer();

  void   SetGLDesignDistanceFromBackFocalPlaneToDetector(double distance);
  double GetGLDesignDistanceFromBackFocalPlaneToDetector();
  void   SetGLActualDistanceFromBackFocalPlaneToDetector(double distance);
  double GetGLActualDistanceFromBackFocalPlaneToDetector();

  void   SetGLIntensityShift(double intensity);
  double GetGLIntensityShift();

  void   SetGLIntensityScale(double intensity);
  double GetGLIntensityScale();

  void   SetGLParameterEnabled(unsigned int index, bool enabled);
  bool   GetGLParameterEnabled(unsigned int index);

  void   SetZCoordinate(unsigned int index, double coordinate);
  double GetZCoordinate(unsigned int index);

  void   SetUseCustomZCoordinates(bool use);
  bool   GetUseCustomZCoordinates();

  double GetImageComparisonMetricValue();

  void Optimize();

protected:
  Configuration m_Configuration;

  std::string m_ImageFileName;

  TImage::Pointer m_MeasuredImageData;

  GibsonLanniPSFImageSourcePointer m_GibsonLanniPSFSource;
  BeadSpreadFunctionImageSourcePointer m_BeadSpreadFunctionSource;
  DifferenceFilterPointer          m_BSFDifferenceImageFilter;

  MinMaxType::Pointer m_MeasuredImageMinMaxFilter;
  MinMaxType::Pointer m_PSFImageMinMaxFilter;
  MinMaxType::Pointer m_BSFImageMinMaxFilter;
  MinMaxType::Pointer m_BSFDifferenceImageMinMaxFilter;

  ITKImageToVTKImage<TImage>* m_MeasuredImageITKToVTKFilter;
  ITKImageToVTKImage<TImage>* m_PSFImageITKToVTKFilter;
  ITKImageToVTKImage<TImage>* m_BSFImageITKToVTKFilter;
  ITKImageToVTKImage<TImage>* m_BSFDifferenceImageITKToVTKFilter;

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
