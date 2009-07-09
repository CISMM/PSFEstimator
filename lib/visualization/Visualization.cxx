#include "Visualization.h"

#include <vtkCamera.h>


Visualization
::Visualization() {
  m_OutlineVisualization = new OutlineVisualizationPipeline();
  m_XPlane = new ImagePlaneVisualizationPipeline();
  m_XPlane->SetToXPlane();
  m_XPlane->SetSliceNumber(0);
  m_XPlane->SetAutoScalingOff();

  m_YPlane = new ImagePlaneVisualizationPipeline();
  m_YPlane->SetToYPlane();
  m_YPlane->SetSliceNumber(0);
  m_YPlane->SetAutoScalingOff();

  m_ZPlane = new ImagePlaneVisualizationPipeline();
  m_ZPlane->SetToZPlane();
  m_ZPlane->SetSliceNumber(0);
  m_ZPlane->SetAutoScalingOff();
}


Visualization
::~Visualization() {
  delete m_OutlineVisualization;
  delete m_XPlane;
  delete m_YPlane;
  delete m_ZPlane;
}


void
Visualization
::SetImageInputConnection(vtkAlgorithmOutput* input) {
  input->GetProducer()->Modified();
  input->GetProducer()->Update();
  m_OutlineVisualization->SetInputConnection(input);
  m_XPlane->SetInputConnection(input);
  m_YPlane->SetInputConnection(input);
  m_ZPlane->SetInputConnection(input);
}


void
Visualization
::AddToRenderer() {
  m_OutlineVisualization->AddToRenderer(m_Renderer);
  m_XPlane->AddToRenderer(m_Renderer);
  m_YPlane->AddToRenderer(m_Renderer);
  m_ZPlane->AddToRenderer(m_Renderer);
}


void
Visualization
::SetShowOutline(bool show) {
  m_OutlineVisualization->SetVisible(show);
}


bool
Visualization
::GetShowOutline() {
  return m_OutlineVisualization->GetVisible();
}


void
Visualization
::ShowOutlineOn() {
  SetShowOutline(true);
}


void
Visualization
::ShowOutlineOff() {
  SetShowOutline(false);
}


void
Visualization
::SetShowXPlane(bool show) {
  m_XPlane->SetVisible(show);
}


bool
Visualization
::GetShowXPlane() {
  return m_XPlane->GetVisible();
}


void
Visualization
::SetShowYPlane(bool show) {
  m_YPlane->SetVisible(show);
}


bool
Visualization
::GetShowYPlane() {
  return m_YPlane->GetVisible();
}


void
Visualization
::SetShowZPlane(bool show) {
  m_ZPlane->SetVisible(show);
}


bool
Visualization
::GetShowZPlane() {
  return m_ZPlane->GetVisible();
}


void
Visualization
::SetXPlane(int plane) {
  m_XPlane->SetSliceNumber(plane);
}


int
Visualization
::GetXPlane() {
  return m_XPlane->GetSliceNumber();
}


void
Visualization
::SetYPlane(int plane) {
  m_YPlane->SetSliceNumber(plane);
}


int
Visualization
::GetYPlane() {
  return m_YPlane->GetSliceNumber();
}


void
Visualization
::SetZPlane(int plane) {
  m_ZPlane->SetSliceNumber(plane);
}


int
Visualization
::GetZPlane() {
  return m_ZPlane->GetSliceNumber();
}


void
Visualization
::SetImagePlanesBlackValue(double value) {
  m_XPlane->SetMapsToBlack(value);
  m_YPlane->SetMapsToBlack(value);
  m_ZPlane->SetMapsToBlack(value);
}


double
Visualization
::GetImagePlanesBlackValue() {
  return m_XPlane->GetMapsToBlack();
}


void
Visualization
::SetImagePlanesWhiteValue(double value) {
  m_XPlane->SetMapsToWhite(value);
  m_YPlane->SetMapsToWhite(value);
  m_ZPlane->SetMapsToWhite(value);
}


double
Visualization
::GetImagePlanesWhiteValue() {
  return m_XPlane->GetMapsToWhite();
}


void
Visualization
::ResetView() {
  vtkCamera* camera = m_Renderer->GetActiveCamera();
  camera->SetFocalPoint(0, 0, 0);
  camera->SetPosition(0, 0, 1);
  camera->SetViewUp(0, 1, 0);
  m_Renderer->ResetCamera();
}


void 
Visualization
::SetViewToXPlus() {
  ResetView();
  vtkCamera* camera = m_Renderer->GetActiveCamera();
  camera->Azimuth(-90.0);
}


void
Visualization
::SetViewToXMinus() {
  ResetView();
  vtkCamera* camera = m_Renderer->GetActiveCamera();
  camera->Azimuth(90.0);
}


void
Visualization
::SetViewToYPlus() {
  ResetView();
  vtkCamera* camera = m_Renderer->GetActiveCamera();
  camera->Elevation(-90.0);
}


void
Visualization
::SetViewToYMinus() {
  ResetView();
  vtkCamera* camera = m_Renderer->GetActiveCamera();
  camera->Elevation(90.0);
}


void
Visualization
::SetViewToZPlus() {
  ResetView();
  vtkCamera* camera = m_Renderer->GetActiveCamera();
  camera->Azimuth(180.0);
}


void
Visualization
::SetViewToZMinus() {
  ResetView();
  // No rotation needed.
}


void
Visualization
::Update() {
  m_XPlane->Update();
  m_YPlane->Update();
  m_ZPlane->Update();
}
