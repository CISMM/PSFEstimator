#ifndef _PSF_ESTIMATOR_H_
#define _PSF_ESTIMATOR_H_

#include <QDialog>
#include <QErrorMessage>
#include <QMainWindow>
#include <QStandardItemModel>

#include <QPointSpreadFunctionPropertyTableModel.h>

#include "ui_PSFEstimator.h"
#include "ui_NewImageDialog.h"

// Forward class declarations
class DataModel;
class Visualization;
class vtkRenderer;
class QCloseEvent;

class PSFEstimator : public QMainWindow {
  Q_OBJECT

public:

  // Constructor/destructor
  PSFEstimator(QWidget* parent=0);
  virtual ~PSFEstimator();

public slots:

  // Use Qt's auto-connect magic to tie GUI widgets to slots.
  // Names of the methods must follow the naming convention
  // on_<widget name>_<signal name>(<signal parameters>).
  virtual void on_actionNewImage_triggered();
  virtual void on_actionOpenImage_triggered();
  virtual void on_actionSavePSFImage_triggered();
  virtual void on_actionSaveBSFImage_triggered();
  virtual void on_actionLoadSession_triggered();
  virtual void on_actionSaveSession_triggered();
  virtual void on_actionExit_triggered();

  virtual void on_actionCopy_triggered();
  virtual void on_actionPaste_triggered();

  virtual void on_actionAboutApplication_triggered();

  virtual void on_measuredBSFRadioButton_clicked(bool state);
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

  virtual void on_useCustomZSlicePositions_toggled(bool use);

  virtual void on_resetCustomSlicePositionsButton_clicked();
  virtual void on_estimatePSFCenterButton_clicked();
  virtual void on_applyButton_clicked();
  virtual void on_optimizePSFParametersButton_clicked();
  virtual void on_submitOptimizationJobToQueueButton_clicked();

  virtual void handle_imageInformationTableModel_dataChanged
    (const QModelIndex& topLeft, const QModelIndex& bottomRight);

  virtual void handle_PSFPropertyTableModel_dataChanged
    (const QModelIndex& topLeft, const QModelIndex& bottomRight);

  /** Mark the session state as changed. */
  void Sully();


protected:
  Ui_MainWindow* gui;

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

  QPointSpreadFunctionPropertyTableModel* m_PSFPropertyTableModel;

  DisplayImageType m_DisplayedImage;

  void SetupInterface(bool hasMeasuredImage);
  void SetupRenderer();

  void SetDisplayedImageToMeasuredPSF();
  void SetDisplayedImageToCalculatedPSF();
  void SetDisplayedImageToCalculatedBSF();

  void RefreshUI();

  double GetDisplayedImageDataMinimum();
  double GetDisplayedImageDataMaximum();

  void SetMapsToBlackValueFromSliderPosition(int position);
  void SetMapsToWhiteValueFromSliderPosition(int position);

  void Exit();
  void WriteProgramSettings();
  void ReadProgramSettings();

  void    SaveFileChooserDirectory(const QString& path);
  QString GetFileChooserDirectory();


  // Override the closeEvent handler.
  void closeEvent(QCloseEvent* event);

protected slots:

private:
  vtkRenderer*  m_Renderer;

  QDialog       m_NewFileDialog;
  Ui::Dialog    m_NewFileDialogUI;

  QErrorMessage m_ErrorDialog;

};

#endif // _PSF_ESTIMATOR_H_
