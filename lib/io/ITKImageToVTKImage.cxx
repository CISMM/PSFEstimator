#ifndef _ITK_IMAGE_TO_VTK_IMAGE_CXX_
#define _ITK_IMAGE_TO_VTK_IMAGE_CXX_

#include "ITKImageToVTKImage.h"

#include <vtkImageImport.h>

template <class TImage>
ITKImageToVTKImage<TImage>
::ITKImageToVTKImage() {
  this->exporter = itk::VTKImageExport<TImage>::New();

  this->importer = vtkImageImport::New();
  this->importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
  this->importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
  this->importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
  this->importer->SetSpacingCallback(exporter->GetSpacingCallback());
  this->importer->SetOriginCallback(exporter->GetOriginCallback());
  this->importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
  this->importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
  this->importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
  this->importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
  this->importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
  this->importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
  this->importer->SetCallbackUserData(exporter->GetCallbackUserData());
}


template <class TImage>
ITKImageToVTKImage<TImage>
::~ITKImageToVTKImage() {
  this->importer->Delete();
}


template <class TImage>
void
ITKImageToVTKImage<TImage>
::SetInput(typename TImage::Pointer input) {
  this->exporter->SetInput(input);
}


template <class TImage>
vtkAlgorithmOutput*
ITKImageToVTKImage<TImage>
::GetOutputPort() {
  return this->importer->GetOutputPort();
}


template <class TImage>
void
ITKImageToVTKImage<TImage>
::Modified() {
  this->importer->Modified();
}


template <class TImage>
void
ITKImageToVTKImage<TImage>
::Update() {
  this->importer->Update();
}

#endif // _ITK_IMAGE_TO_VTK_IMAGE_CXX_ 
