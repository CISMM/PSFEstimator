#include "OutlineVisualizationPipeline.h"


OutlineVisualizationPipeline
::OutlineVisualizationPipeline() {
  this->outlineFilter = vtkOutlineFilter::New();

  this->SetInputAlgorithm(this->outlineFilter);

  this->outlineMapper = vtkPolyDataMapper::New();
  this->outlineMapper->SetInputConnection(this->outlineFilter->GetOutputPort());
  this->outlineActor = vtkActor::New();
  this->outlineActor->SetMapper(this->outlineMapper);
}


OutlineVisualizationPipeline
::~OutlineVisualizationPipeline() {
  this->outlineFilter->Delete();
  this->outlineMapper->Delete();
  this->outlineActor->Delete();
}


void
OutlineVisualizationPipeline
::AddToRenderer(vtkRenderer* renderer) {
  renderer->AddActor(this->outlineActor);
}


void
OutlineVisualizationPipeline
::SetVisible(bool visible) {
  this->outlineActor->SetVisibility(static_cast<int>(visible));
}


bool
OutlineVisualizationPipeline
::GetVisible() {
  return this->outlineActor->GetVisibility() == 1;
}

