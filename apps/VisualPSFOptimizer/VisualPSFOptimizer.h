#ifndef _QTVTKITK_GENERIC_APPLICATION_H_
#define _QTVTKITK_GENERIC_APPLICATION_H_

#include <qerrormessage.h>
#include <qmainwindow.h>
#include <qstandarditemmodel.h>

#include "ui_VisualPSFOptimizer.h"

// Forward class declarations
class DataModel;
class Visualization;
class vtkRenderer;

class VisualPSFOptimizer : public QMainWindow, private Ui_MainWindow {
  Q_OBJECT

public:

  // Constructor/destructor
  VisualPSFOptimizer(QWidget* parent=0);
  virtual ~VisualPSFOptimizer();
  
public slots:

  // Use Qt's auto-connect magic to tie GUI widgets to slots.
  // Names of the methods must follow the naming convention
  // on_<widget name>_<signal name>(<signal parameters>).
  virtual void on_actionOpenImage_triggered();
  virtual void on_actionSavePSFImage_triggered();
  virtual void on_actionLoadSettings_triggered();
  virtual void on_actionSaveSettings_triggered();
  virtual void on_actionExit_triggered();
  
  virtual void on_actionAboutApplication_triggered();

  virtual void on_measuredPSFRadioButton_clicked(bool state);
  virtual void on_calculatedPSFRadioButton_clicked(bool state);
  virtual void on_calculatedBSFRadioButton_clicked(bool state);
  
  virtual void on_showXPlaneCheckBox_toggled(bool show);
  virtual void on_xPlaneSlider_valueChanged(int value);
  virtual void on_xPlaneEdit_textEdited(QString text);
  
  virtual void on_showYPlaneCheckBox_toggled(bool show);
  virtual void on_yPlaneSlider_valueChanged(int value);
  virtual void on_yPlaneEdit_textEdited(QString text);
  
  virtual void on_showZPlaneCheckBox_toggled(bool show);
  virtual void on_zPlaneSlider_valueChanged(int value);
  virtual void on_zPlaneEdit_textEdited(QString text);
  
  virtual void on_mapsToBlackSlider_valueChanged(int value);
  virtual void on_mapsToWhiteSlider_valueChanged(int value);

  virtual void on_showDataOutlineCheckBox_toggled(bool show);

  virtual void on_xPlusButton_clicked();
  virtual void on_xMinusButton_clicked();
  virtual void on_yPlusButton_clicked();
  virtual void on_yMinusButton_clicked();
  virtual void on_zPlusButton_clicked();
  virtual void on_zMinusButton_clicked();

  virtual void on_applyButton_clicked();
  virtual void on_estimatePSFCenterButton_clicked();
  virtual void on_optimizePSFParametersButton_clicked();

  virtual void handle_imageInformationTableModel_dataChanged(const QModelIndex& topLeft,
    const QModelIndex& bottomRight);

  /** Mark the session state as changed. */
  void Sully();


protected:
  /** Dirty bit on the session. */
  bool m_Dirty;

  typedef enum { 
    MEASURED_PSF_IMAGE,
    CALCULATED_PSF_IMAGE,
    CALCULATED_BSF_IMAGE
  } DisplayImageType;

  DataModel* m_DataModel;
  
  Visualization* m_Visualization;
  
  QStandardItemModel* m_ImageInformationTableModel;

  QStandardItemModel* m_GibsonLanniPSFSettingsTableModel;
  
  DisplayImageType m_DisplayedImage;

  void OpenFile(std::string fileName);
  
  void SetDisplayedImageToMeasuredPSF();
  void SetDisplayedImageToCalculatedPSF();
  void SetDisplayedImageToCalculatedBSF();

  void RefreshUI();

  double GetDisplayedImageDataMinimum();
  double GetDisplayedImageDataMaximum();

  void SetMapsToBlackValueFromSliderPosition(int position);
  void SetMapsToWhiteValueFromSliderPosition(int position);

  void writeProgramSettings();
  void readProgramSettings();

protected slots:

private:
  vtkRenderer* m_Renderer;
  
  QErrorMessage m_ErrorDialog;
  
};

#endif // _QTVTKITK_GENERIC_APPLICATION_H_
