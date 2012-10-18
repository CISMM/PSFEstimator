#ifndef _DATA_MODEL_H_
#define _DATA_MODEL_H_

#include <string>

#include "Configuration.h"

#include "Validation.h"

#define ITK_MANUAL_INSTANTIATION

#ifdef VALIDATE_CONVOLUTION
#include <itkBeadSpreadFunctionImageSource2.h>
#else
#include <itkBeadSpreadFunctionImageSource.h>
#endif

// IO
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// PSF sources
#include <itkGaussianImageSource.h>
#include <itkMaskedParametricImageSource.h>
#include <itkGibsonLanniCOSMOSPointSpreadFunctionImageSource.h>
#include <itkHaeberleCOSMOSPointSpreadFunctionImageSource.h>

// Optimizers
#include <itkAmoebaOptimizer.h>
#include <itkConjugateGradientOptimizer.h>
#include <itkGradientDescentOptimizer.h>
#include <itkLBFGSBOptimizer.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkPowellOptimizer.h>

// Metrics
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
//#include <itkPoissonNoiseImageToImageMetric.h>
#include <itkImageToParametricImageSourceMetric.h>

// Misc
#include <itkGridImageSource.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkShiftScaleImageFilter.h>
#include <itkSubtractImageFilter.h>

#include <ITKImageToVTKImage.h>
#undef ITK_MANUAL_INSTANTIATION

#include <vtkAlgorithmOutput.h>


// This is the data model for the PSF Optimizer library.
class DataModel {

public:
  typedef enum {
    GAUSSIAN_PSF = 0,
    GIBSON_LANNI_PSF,
    HAEBERLE_PSF,
    NUM_PSFS
  } PointSpreadFunctionType;

  typedef enum {
    MEAN_SQUARED_ERROR = 0,
    NORMALIZED_CORRELATION,
    NUM_OBJECTIVE_FUNCTIONS
  } ObjectiveFunctionType;

  typedef enum {
    AMOEBA_OPTIMIZER = 0,
    GRADIENT_DESCENT_OPTIMIZER,
    ONE_PLUS_ONE_EVOLUTIONARY_OPTIMIZER,
    NUM_OPTIMIZERS
  } OptimizerType;

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

  typedef itk::ParametricImageSource<Float3DImageType>
    ParametricImageSourceType;
  typedef ParametricImageSourceType::Pointer
    ParametricImageSourcePointer;

  typedef itk::GaussianImageSource<Float3DImageType>
    GaussianPSFImageSourceType;
  typedef GaussianPSFImageSourceType::Pointer
    GaussianPSFImageSourcePointer;

  typedef itk::MaskedParametricImageSource<Float3DImageType>
    MaskedGaussianPSFImageSourceType;
  typedef MaskedGaussianPSFImageSourceType::Pointer
    MaskedGaussianPSFImageSourcePointer;

  typedef itk::GibsonLanniCOSMOSPointSpreadFunctionImageSource<Float3DImageType>
    GibsonLanniPSFImageSourceType;
  typedef GibsonLanniPSFImageSourceType::Pointer
    GibsonLanniPSFImageSourcePointer;

  typedef itk::HaeberleCOSMOSPointSpreadFunctionImageSource<Float3DImageType>
    HaeberlePSFImageSourceType;
  typedef HaeberlePSFImageSourceType::Pointer
    HaeberlePSFImageSourcePointer;

#ifdef VALIDATE_CONVOLUTION
  typedef itk::BeadSpreadFunctionImageSource2< Float3DImageType >
    BeadSpreadFunctionImageSourceType;
#else
  typedef itk::BeadSpreadFunctionImageSource< Float3DImageType >
    BeadSpreadFunctionImageSourceType;
#endif

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
  typedef ParametricCostFunctionType::ParametersMaskType ParametersMaskType;
  typedef ParametricCostFunctionType::ParametersType     ParametersType;

  typedef itk::SingleValuedNonLinearOptimizer  OptimizerBaseType;
  typedef itk::AmoebaOptimizer                 AmoebaOptimizerType;
  typedef itk::GradientDescentOptimizer        GradientDescentOptimizerType;
  typedef itk::OnePlusOneEvolutionaryOptimizer OnePlusOneEvolutionaryOptimizerType;

  typedef itk::NearestNeighborInterpolateImageFunction<TImage, double>
    InterpolatorType;
  //  typedef itk::PoissonNoiseImageToImageMetric<TImage, TImage>
  typedef itk::ImageToImageMetric< TImage, TImage >
    ImageToImageCostFunctionBaseType;
  typedef itk::MeanSquaresImageToImageMetric< TImage, TImage >
    MeanSquaredCostFunctionType;
  typedef itk::NormalizedCorrelationImageToImageMetric< TImage, TImage >
    NormalizedCorrelationCostFunctionType;

  DataModel();
  virtual ~DataModel();

  void SetPointSpreadFunctionType(PointSpreadFunctionType psfType);
  PointSpreadFunctionType GetPointSpreadFunctionType() const;

  void                  SetObjectiveFunctionType(ObjectiveFunctionType ofType);
  ObjectiveFunctionType GetObjectiveFunctionType() const;

  void SetOptimizerType(OptimizerType optimizerType);
  OptimizerType GetOptimizerType() const;

  bool LoadSessionFile(const std::string& fileName);
  bool SaveSessionFile(const std::string& fileName);

  void Initialize();

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

  // Recenters the origins of the images to the center of the image bounds.
  void   RecenterImageOrigin();

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

  void   UpdatePSFImage();
  void   UpdateBSFImage();
  void   UpdateBSFDifferenceImage();

  void   SetBeadRadius(double radius);
  double GetBeadRadius();

  void   SetIntensityShift(double intensity);
  double GetIntensityShift();

  void   SetIntensityScale(double intensity);
  double GetIntensityScale();

  std::string GetParameterName(unsigned int index) const;

  void   SetParameterValue(unsigned int index, double value);
  double GetParameterValue(unsigned int index);

  std::string GetParameterUnit(unsigned int index) const;

  void   SetParameterEnabled(unsigned int index, bool enabled);
  bool   GetParameterEnabled(unsigned int index);

  void   SetParameterScale(unsigned int index, double scale);
  double GetParameterScale(unsigned int index);

  void   SetZCoordinate(unsigned int index, double coordinate);
  double GetZCoordinate(unsigned int index);

  void   SetUseCustomZCoordinates(bool use);
  bool   GetUseCustomZCoordinates();

  void   UpdateMetricParameterMask();

  double GetImageComparisonMetricValue();

  void SetUpOptimizer(const ParametersType & parameterScales);
  void Optimize();

protected:
  PointSpreadFunctionType m_PointSpreadFunctionType;
  ObjectiveFunctionType   m_ObjectiveFunctionType;
  OptimizerType           m_OptimizerType;

  Configuration m_Configuration;

  std::string m_ImageFileName;

  TImage::Pointer m_MeasuredImageData;

  // The different point-spread function types
  ParametricImageSourcePointer             m_PointSpreadFunctionSource;
  MaskedGaussianPSFImageSourcePointer      m_GaussianPSFSource;
  MaskedGaussianPSFImageSourcePointer      m_GaussianPSFKernelSource;
  GibsonLanniPSFImageSourcePointer         m_GibsonLanniPSFSource;
  GibsonLanniPSFImageSourcePointer         m_GibsonLanniPSFKernelSource;
  HaeberlePSFImageSourcePointer            m_HaeberlePSFSource;
  HaeberlePSFImageSourcePointer            m_HaeberlePSFKernelSource;

  // Lists of parameter names for the BSF and the different PSFs
  std::vector<std::string> m_BSFParameterNames;
  std::vector<std::string> m_GaussianPSFParameterNames;
  std::vector<std::string> m_GibsonLanniPSFParameterNames;
  std::vector<std::string> m_HaeberlePSFParameterNames;

  // Lists of units for the BSF and the different PSFs
  std::vector<std::string> m_BSFParameterUnits;
  std::vector<std::string> m_GaussianPSFParameterUnits;
  std::vector<std::string> m_GibsonLanniPSFParameterUnits;
  std::vector<std::string> m_HaeberlePSFParameterUnits;

  // Parameter optimization masks
  std::vector<bool> m_BSFParameterMask;
  std::vector<bool> m_GaussianPSFParameterMask;
  std::vector<bool> m_GibsonLanniPSFParameterMask;
  std::vector<bool> m_HaeberlePSFParameterMask;

  // Parameter scales. It is important to set these appropriately
  // because small changes to some parameters (e.g. actual refractive
  // index) produce huge changes in the shape of the PSF.
  std::vector<double> m_BSFParameterScales;
  std::vector<double> m_GaussianPSFParameterScales;
  std::vector<double> m_GibsonLanniPSFParameterScales;
  std::vector<double> m_HaeberlePSFParameterScales;

  // The bead-spread function
  BeadSpreadFunctionImageSourcePointer m_BeadSpreadFunctionSource;
  DifferenceFilterPointer              m_BSFDifferenceImageFilter;

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
  ImageToImageCostFunctionBaseType::Pointer  m_ImageToImageCostFunction;

  // The optimizer
  OptimizerBaseType::Pointer m_Optimizer;

  // Objective function and optimizer names
  std::vector< std::string > m_ObjectiveFunctionNames;
  std::vector< std::string > m_OptimizerNames;

};

// _DATA_MODEL_H_
#endif
