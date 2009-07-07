#include "IsosurfaceVisualizationPipeline.h"

#include <vtkActor.h>
#include <vtkAlgorithmOutput.h>
#include <vtkContourFilter.h>
#include <vtkImageClip.h>
#include <vtkImageData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>


IsosurfaceVisualizationPipeline
::IsosurfaceVisualizationPipeline() {
  this->isoValue = 0.0;
  this->zPlane = 0;
  this->deltaZ = 2;
  
  // Set up pipeline.
  this->imageClip = vtkImageClip::New();
  this->imageClip->ClipDataOff();

  this->isoContourer = vtkContourFilter::New();
  this->isoContourer->SetNumberOfContours(1);
  this->isoContourer->SetValue(0, this->isoValue);
  this->isoContourer->SetInputConnection(this->imageClip->GetOutputPort());

  this->isoMapper = vtkPolyDataMapper::New();
  this->isoMapper->ScalarVisibilityOff();
  this->isoMapper->ImmediateModeRenderingOff();
  this->isoMapper->SetInputConnection(this->isoContourer->GetOutputPort()); //

  this->isoActor = vtkActor::New();
  this->isoActor->SetMapper(this->isoMapper);

  // It is essential to set the inputAlgorithm
  this->SetInputAlgorithm(imageClip);

  // Connect the input to the input algorithm.
  //this->isoMapper->SetInputConnection(this->inputAlgorithm->GetOutputPort());
  //this->imageClip->SetInputConnection(this->inputAlgorithm->GetOutputPort());
}


IsosurfaceVisualizationPipeline
::~IsosurfaceVisualizationPipeline() {
}


void
IsosurfaceVisualizationPipeline
::AddToRenderer(vtkRenderer* renderer) {
  this->isoContourer->Modified();
  renderer->AddActor(this->isoActor);
}


void
IsosurfaceVisualizationPipeline
::SetColor(double r, double g, double b) {
  this->isoActor->GetProperty()->SetColor(r, g, b);
}


void
IsosurfaceVisualizationPipeline
::SetVisible(bool visible) {
  this->isoActor->SetVisibility(visible ? 1 : 0);
}


bool
IsosurfaceVisualizationPipeline
::GetVisible() {
  return this->isoActor->GetVisibility() == 1;
}


void
IsosurfaceVisualizationPipeline
::ClipDataOn() {
  this->imageClip->ClipDataOn();
}


void
IsosurfaceVisualizationPipeline
::ClipDataOff() {
  this->imageClip->ClipDataOff();
}


void
IsosurfaceVisualizationPipeline
::SetClipData(bool clip) {
  this->imageClip->SetClipData(clip ? 1 : 0);
}


bool
IsosurfaceVisualizationPipeline
::GetClipData() {
  return this->imageClip->GetClipData() == 1;
}


void
IsosurfaceVisualizationPipeline
::SetZPlane(int zPlane) {
  this->zPlane = zPlane;
  int dims[3];
  this->imageClip->GetImageDataInput(0)->GetDimensions(dims);

  if (this->imageClip->GetClipData()) {
    int minZ = zPlane - this->deltaZ;
    if (minZ < 0) 
      minZ = 0;
    int maxZ = zPlane + this->deltaZ;
    if (maxZ >= dims[2]) 
      maxZ = dims[2]-1;
    this->imageClip->SetOutputWholeExtent(0, dims[0]-1, 0, dims[1]-1, minZ, maxZ);
    std::cout << deltaZ << ", " << minZ << ", " << maxZ << std::endl;
  } else {
    this->imageClip->SetOutputWholeExtent(0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);
  }

}


int
IsosurfaceVisualizationPipeline
::GetZPlane() {
  return this->zPlane;
}


void
IsosurfaceVisualizationPipeline
::SetDeltaZ(int deltaZ) {
  this->deltaZ = deltaZ;
}


int
IsosurfaceVisualizationPipeline
::GetDeltaZ() {
  return this->deltaZ;
}


void
IsosurfaceVisualizationPipeline
::SetIsoValue(double isoValue) {
  this->isoValue = isoValue;
  this->isoContourer->SetValue(0, this->isoValue);
  this->isoContourer->Modified();
}


double
IsosurfaceVisualizationPipeline
::GetIsoValue() {
  return this->isoValue;
}


vtkAlgorithmOutput*
IsosurfaceVisualizationPipeline
::GetIsosurfaceOutputPort() {
  return this->isoContourer->GetOutputPort();
}


void
IsosurfaceVisualizationPipeline
::FastRenderingOn() {
  this->isoMapper->ImmediateModeRenderingOff();
}


void
IsosurfaceVisualizationPipeline
::FastRenderingOff() {
  this->isoMapper->ImmediateModeRenderingOn();
}


bool
IsosurfaceVisualizationPipeline
::GetFastRenderingOn() {
  return this->isoMapper->GetImmediateModeRendering() == 0;
}

