#ifndef _ITK_IMAGE_TO_VTK_IMAGE_H_
#define _ITK_IMAGE_TO_VTK_IMAGE_H_

#include <itkImage.h>
#include <itkVTKImageExport.h>

class vtkImageImport;
class vtkAlgorithmOutput;

template <class TImage>
class ITKImageToVTKImage {
	
public:
  ITKImageToVTKImage();
  virtual ~ITKImageToVTKImage();

  void SetInput(typename TImage::Pointer input);

  vtkAlgorithmOutput* GetOutputPort();

  void Modified();

  void Update();
	
protected:
  typename itk::VTKImageExport<TImage>::Pointer exporter;

  vtkImageImport* importer;
	
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "ITKImageToVTKImage.cxx"
#endif

#endif // _ITK_IMAGE_TO_VTK_IMAGE_H_
