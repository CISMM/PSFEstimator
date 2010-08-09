#ifndef _IMAGE_PLANE_VISUALIZATION_PIPELINE_H_
#define _IMAGE_PLANE_VISUALIZATION_PIPELINE_H_

#include "VisualizationPipeline.h"

#include <vtkSmartPointer.h>

class vtkImageActor;
class vtkImageShiftScale;
class vtkRenderer;

class ImagePlaneVisualizationPipeline : public VisualizationPipeline {

public:
  ImagePlaneVisualizationPipeline();
  ~ImagePlaneVisualizationPipeline();
  
  const static int X_PLANE;
  const static int Y_PLANE;
  const static int Z_PLANE;

  void SetInputConnection(vtkAlgorithmOutput* input);

  void AddToRenderer(vtkRenderer* renderer);
  
  void SetToXPlane();
  void SetToYPlane();
  void SetToZPlane();

  void SetVisible(bool visible);
  bool GetVisible();

  void SetSliceNumber(int sliceNumber);
  int  GetSliceNumber();

  void SetAutoScalingOn();
  void SetAutoScalingOff();

  void   SetMapsToBlack(double value);
  double GetMapsToBlack();
  void   SetMapsToWhite(double value);
  double GetMapsToWhite();

  void Update();

protected:
  int    m_PlaneDimension;
  bool   m_AutoScalingOn;
  double m_MapsToWhite;
  double m_MapsToBlack;

  vtkSmartPointer<vtkImageShiftScale> m_ShiftScaler;

  vtkSmartPointer<vtkImageActor> m_ImageActor;

};

// _IMAGE_PLANE_VISUALIZATION_PIPELINE_H_
#endif
