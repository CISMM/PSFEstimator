#include "ImagePlaneVisualizationPipeline.h"

#include <vtkAlgorithmOutput.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageShiftScale.h>
#include <vtkRenderer.h>


const int ImagePlaneVisualizationPipeline::X_PLANE = 1;
const int ImagePlaneVisualizationPipeline::Y_PLANE = 2;
const int ImagePlaneVisualizationPipeline::Z_PLANE = 3;


ImagePlaneVisualizationPipeline
::ImagePlaneVisualizationPipeline() {
  m_PlaneDimension = X_PLANE;
  m_AutoScalingOn = true;
  m_MapsToWhite = 0.0;
  m_MapsToBlack = 1.0;

  // Set up pipeline.
  m_ShiftScaler = vtkSmartPointer<vtkImageShiftScale>::New();
  m_ShiftScaler->SetOutputScalarTypeToUnsignedChar();
  m_ShiftScaler->ClampOverflowOn();

  // It is essential to set the input algorithm.
  SetInputAlgorithm(m_ShiftScaler);

  m_ImageActor = vtkSmartPointer<vtkImageActor>::New();
  m_ImageActor->InterpolateOff();
  m_ImageActor->SetInput(m_ShiftScaler->GetOutput());
}


ImagePlaneVisualizationPipeline
::~ImagePlaneVisualizationPipeline() {
}


void
ImagePlaneVisualizationPipeline
::SetInputConnection(vtkAlgorithmOutput* input) {
  input = input;
  inputAlgorithm->SetInputConnection(input);

    // Update shift/scale filter automatically if auto-rescaling is on.
  double scalarRange[2];
  if (m_AutoScalingOn) {
    // Gotta be a better way to do this.
    input->GetProducer()->Update();
    vtkImageData* originalImage 
      = vtkImageData::SafeDownCast(input->GetProducer()->GetOutputDataObject(0));
    originalImage->GetScalarRange(scalarRange);
  } else {
    scalarRange[0] = m_MapsToBlack;
    scalarRange[1] = m_MapsToWhite;
  }

  m_ShiftScaler->SetShift(-scalarRange[0]);
  m_ShiftScaler->SetScale(255.0 / (scalarRange[1] - scalarRange[0]));
  m_ShiftScaler->Update();
}


void
ImagePlaneVisualizationPipeline
::AddToRenderer(vtkRenderer* renderer) {
  renderer->AddActor(m_ImageActor);
}


void
ImagePlaneVisualizationPipeline
::SetToXPlane() {
  m_PlaneDimension = X_PLANE;
  SetSliceNumber(0);
}


void
ImagePlaneVisualizationPipeline
::SetToYPlane() {
  m_PlaneDimension = Y_PLANE;
  SetSliceNumber(0);
}


void
ImagePlaneVisualizationPipeline
::SetToZPlane() {
  m_PlaneDimension = Z_PLANE;
  SetSliceNumber(0);
}


void
ImagePlaneVisualizationPipeline
::SetVisible(bool visible) {
  m_ImageActor->SetVisibility(visible ? 1 : 0);
}


bool
ImagePlaneVisualizationPipeline
::GetVisible() {
  return m_ImageActor->GetVisibility() == 1;
}


void
ImagePlaneVisualizationPipeline
::SetSliceNumber(int sliceNumber) {
  int* extent = m_ImageActor->GetInput()->GetWholeExtent();
  
  if (m_PlaneDimension == X_PLANE) {
    m_ImageActor->SetDisplayExtent(sliceNumber, sliceNumber,
                 extent[2], extent[3],
                 extent[4], extent[5]);
  
  } else if (m_PlaneDimension == Y_PLANE) {
    m_ImageActor->SetDisplayExtent(extent[0], extent[1],
                 sliceNumber, sliceNumber,
                 extent[4], extent[5]);
  } else {
    m_ImageActor->SetDisplayExtent(extent[0], extent[1],
                 extent[2], extent[3],
                 sliceNumber, sliceNumber);
  }
  m_ImageActor->Modified();
}


int
ImagePlaneVisualizationPipeline
::GetSliceNumber() {
  int slice = 0;
  if (m_PlaneDimension == X_PLANE) {
    slice = m_ImageActor->GetDisplayExtent()[0];
  } else if (m_PlaneDimension == Y_PLANE) {
    slice = m_ImageActor->GetDisplayExtent()[2];
  } else {
    slice = m_ImageActor->GetDisplayExtent()[4];
  }

  return slice;
}


void
ImagePlaneVisualizationPipeline
::SetAutoScalingOn() {
  m_AutoScalingOn = true;
}


void
ImagePlaneVisualizationPipeline
::SetAutoScalingOff() {
  m_AutoScalingOn = false;
}


void
ImagePlaneVisualizationPipeline
::SetMapsToBlack(double value) {
  m_MapsToBlack = value;

  m_ShiftScaler->SetShift(-m_MapsToBlack);
  m_ShiftScaler->SetScale(255.0 / (m_MapsToWhite - m_MapsToBlack));
  m_ShiftScaler->Update();
}


double
ImagePlaneVisualizationPipeline
::GetMapsToBlack() {
  return m_MapsToBlack;
}


void
ImagePlaneVisualizationPipeline
::SetMapsToWhite(double value) {
  m_MapsToWhite = value;

  m_ShiftScaler->SetShift(-m_MapsToBlack);
  m_ShiftScaler->SetScale(255.0 / (m_MapsToWhite - m_MapsToBlack));
  m_ShiftScaler->Update();
}


double
ImagePlaneVisualizationPipeline
::GetMapsToWhite() {
  return m_MapsToWhite;
}


void
ImagePlaneVisualizationPipeline
::Update() {
  m_ShiftScaler->Modified();
}
