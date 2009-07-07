#include "VisualizationPipeline.h"

#include <vtkActor.h>
#include <vtkAlgorithm.h>
#include <vtkPolyDataMapper.h>
#include <vtkOutlineFilter.h>
#include <vtkRenderer.h>


VisualizationPipeline
::VisualizationPipeline() {
}


VisualizationPipeline
::~VisualizationPipeline() {
}


void
VisualizationPipeline
::SetInputConnection(vtkAlgorithmOutput* input) {
  this->input = input;
  this->inputAlgorithm->SetInputConnection(input);
  
  // Need this here for visualizations to update properly.
  this->inputAlgorithm->Update();
}


void
VisualizationPipeline
::SetInputAlgorithm(vtkAlgorithm* algorithm) {
  this->inputAlgorithm = algorithm;
}


vtkAlgorithmOutput*
VisualizationPipeline
::GetInputConnection() {
  return this->input;
}
