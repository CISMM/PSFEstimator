#ifndef _VISUALIZATION_PIPELINE_H_
#define _VISUALIZATION_PIPELINE_H_

#include <vtkAlgorithm.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>


class VisualizationPipeline {

public:
  VisualizationPipeline();
  virtual ~VisualizationPipeline();

  void SetInputConnection(vtkAlgorithmOutput* input);
  vtkAlgorithmOutput* GetInputConnection();

  virtual void AddToRenderer(vtkRenderer* renderer) = 0;
  
protected:
  vtkAlgorithmOutput* input;

  vtkAlgorithm* inputAlgorithm;
  void SetInputAlgorithm(vtkAlgorithm* algorithm);
  
};


#endif // _VISUALIZATION_PIPELINE_H_
