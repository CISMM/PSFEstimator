#ifndef _VISUALIZATION_H_
#define _VISUALIZATION_H_

#include "OutlineVisualizationPipeline.h"
#include "ImagePlaneVisualizationPipeline.h"

#include <vtkAlgorithmOutput.h>
#include <vtkRenderer.h>

class Visualization {

public:
  Visualization();
  virtual ~Visualization();

  void SetImageInputConnection(vtkAlgorithmOutput* input);
  void AddToRenderer(vtkRenderer* renderer);

  void SetShowOutline(bool show);
  bool GetShowOutline();
  void ShowOutlineOn();
  void ShowOutlineOff();
  
  void SetShowXPlane(bool show);
  bool GetShowXPlane();
  void SetShowYPlane(bool show);
  bool GetShowYPlane();
  void SetShowZPlane(bool show);
  bool GetShowZPlane();
  
  void SetXPlane(int plane);
  int  GetXPlane();
  void SetYPlane(int plane);
  int  GetYPlane();
  void SetZPlane(int plane);
  int  GetZPlane();

  void   SetImagePlanesBlackValue(double value);
  double GetImagePlanesBlackValue();
  void   SetImagePlanesWhiteValue(double value);
  double GetImagePlanesWhiteValue();

protected:

  OutlineVisualizationPipeline* m_OutlineVisualization;
  ImagePlaneVisualizationPipeline* m_XPlane;
  ImagePlaneVisualizationPipeline* m_YPlane;
  ImagePlaneVisualizationPipeline* m_ZPlane;

};

#endif // _VISUALIZATION_H_
