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

  void SetRenderer(vtkRenderer* renderer) {
    m_Renderer = renderer;
  }
  vtkRenderer* GetRenderer() {
    return m_Renderer;
  }

  void SetImageInputConnection(vtkAlgorithmOutput* input);
  void AddToRenderer();

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

  void ResetView();
  void SetViewToXPlus();
  void SetViewToXMinus();
  void SetViewToYPlus();
  void SetViewToYMinus();
  void SetViewToZPlus();
  void SetViewToZMinus();

  void Update();

protected:

  vtkRenderer* m_Renderer;
  OutlineVisualizationPipeline* m_OutlineVisualization;
  ImagePlaneVisualizationPipeline* m_XPlane;
  ImagePlaneVisualizationPipeline* m_YPlane;
  ImagePlaneVisualizationPipeline* m_ZPlane;

};

#endif // _VISUALIZATION_H_
