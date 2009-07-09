#ifndef _QTVTKITK_GENERIC_APPLICATION_H_
#define _QTVTKITK_GENERIC_APPLICATION_H_

#include <qerrormessage.h>
#include <qmainwindow.h>
#include <qstandarditemmodel.h>

#include "ui_VisualPSFOptimizer.h"

#include "DataModel.h"
#include "Visualization.h"

// Forward class declarations
class vtkRenderer;

class VisualPSFOptimizer : public QMainWindow, private Ui_MainWindow {
  Q_OBJECT

public:

  // Constructor/destructor
  VisualPSFOptimizer(QWidget* parent=0);
  virtual ~VisualPSFOptimizer();
  
protected:

  void writeProgramSettings();
  void readProgramSettings();
  
public slots:

  // Use Qt's auto-connect magic to tie GUI widgets to slots.
  // Names of the methods must follow the naming convention
  // on_<widget name>_<signal name>(<signal parameters>).
  virtual void on_actionOpenImage_triggered();
  virtual void on_actionSavePSFImage_triggered();
  virtual void on_actionExit_triggered();
  
  virtual void on_actionAboutApplication_triggered();
  
  virtual void on_showXPlaneCheckBox_toggled(bool show);
  virtual void on_xPlaneSlider_sliderMoved(int value);
  virtual void on_xPlaneEdit_textEdited(QString text);
  
  virtual void on_showYPlaneCheckBox_toggled(bool show);
  virtual void on_yPlaneSlider_sliderMoved(int value);
  virtual void on_yPlaneEdit_textEdited(QString text);
  
  virtual void on_showZPlaneCheckBox_toggled(bool show);
  virtual void on_zPlaneSlider_sliderMoved(int value);
  virtual void on_zPlaneEdit_textEdited(QString text);
  
  virtual void on_mapsToBlackSlider_sliderMoved(int value);
  virtual void on_mapsToWhiteSlider_sliderMoved(int value);

  virtual void on_showDataOutlineCheckBox_toggled(bool show);

  virtual void on_xPlusButton_clicked();
  virtual void on_xMinusButton_clicked();
  virtual void on_yPlusButton_clicked();
  virtual void on_yMinusButton_clicked();
  virtual void on_zPlusButton_clicked();
  virtual void on_zMinusButton_clicked();

  virtual void on_applyButton_clicked();
  
  virtual void handle_imageInformationTableModel_dataChanged(const QModelIndex& topLeft,
    const QModelIndex& bottomRight);


protected:
  DataModel* m_DataModel;
  
  Visualization* m_Visualization;
  
  QStandardItemModel* m_ImageInformationTableModel;

  QStandardItemModel* m_GibsonLanniPSFSettingsTableModel;
  
  void OpenFile(std::string fileName);
  
  void RefreshUI();

protected slots:

private:
  vtkRenderer* m_Renderer;
  
  QErrorMessage m_ErrorDialog;
  
};

#endif // _QTVTKITK_GENERIC_APPLICATION_H_
