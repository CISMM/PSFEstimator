#include "VisualPSFOptimizer.h"
#include "Version.h"

#include <qmessagebox.h>

#if defined(_WIN32) // Turn off deprecation warnings in Visual Studio
#pragma warning( disable : 4996 )
#endif

#include "Configuration.h"
#include "DataModel.h"
#include "Visualization.h"

#include <qapplication.h>
#include <qfiledialog.h>
#include <qsettings.h>
#include <qvariant.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#define CLAMP(value, min, max) (value < min ? min : (value > max ? max : value))

// Constructor
VisualPSFOptimizer
::VisualPSFOptimizer(QWidget* p)
 : QMainWindow(p) {
  
  gui = new Ui_MainWindow();
  gui->setupUi(this);

  // Mark as initially clean
  m_Dirty = false;
  m_DisplayedImage = MEASURED_PSF_IMAGE;
  
  // QT/VTK interaction
  m_Renderer = vtkRenderer::New();
  gui->qvtkWidget->GetRenderWindow()->AddRenderer(m_Renderer);
  
  // Instantiate data model.
  m_DataModel = new DataModel();
  
  // Instantiate m_Visualization pipelines.
  m_Visualization = new Visualization();
  m_Visualization->SetRenderer(m_Renderer);
  
  // Set application information
  QCoreApplication::setOrganizationName("CISMM");
  QCoreApplication::setOrganizationDomain("cismm.org");
  QCoreApplication::setApplicationName("VisualPSFOptimizer");
  
  // Restore inter-session GUI settings.
  readProgramSettings();
  
  // Set up dialog boxes.
  m_NewFileDialogUI.setupUi(&m_NewFileDialog);

  m_ErrorDialog.setModal(true);
  
  // Create and populate image information table model.
  int LEFT_COLUMN = 0;
  int RIGHT_COLUMN = 1;
  m_ImageInformationTableModel = new QStandardItemModel(5, 2, this);
  m_ImageInformationTableModel->setHeaderData(LEFT_COLUMN,  Qt::Horizontal, tr("Property"));
  m_ImageInformationTableModel->setHeaderData(RIGHT_COLUMN, Qt::Horizontal, tr("Value"));
  
  QStandardItem* labelItems[5];
  labelItems[ 0] = new QStandardItem(tr("Intensity minimum"));
  labelItems[ 1] = new QStandardItem(tr("Intensity maximum"));
  labelItems[ 2] = new QStandardItem(tr("X dimension (pixels)"));
  labelItems[ 3] = new QStandardItem(tr("Y dimension (pixels)"));
  labelItems[ 4] = new QStandardItem(tr("Z dimension (slices)"));
 
  for (unsigned int i = 0; i < sizeof(labelItems) / sizeof(QStandardItem*); i++) {
    labelItems[i]->setEditable(false);
    m_ImageInformationTableModel->setItem(i, LEFT_COLUMN, labelItems[i]);

    QStandardItem* item = new QStandardItem(tr(""));
    item->setEditable(false);

    m_ImageInformationTableModel->setItem(i, RIGHT_COLUMN, item);
  }
  gui->imageDataView->setModel(m_ImageInformationTableModel);
  
  connect(m_ImageInformationTableModel, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(handle_imageInformationTableModel_dataChanged(const QModelIndex&, const QModelIndex&)));
  
  int LEFT_COLUMN_WIDTH = 160;
  gui->imageDataView->setColumnWidth(LEFT_COLUMN, LEFT_COLUMN_WIDTH);

  // Create and popuplate Gibson-Lanni PSF settings table model.
  int column = 0;
  m_GibsonLanniPSFSettingsTableModel = new QStandardItemModel(24, 4, this);
  m_GibsonLanniPSFSettingsTableModel->
    setHeaderData(column++, Qt::Horizontal, tr("Property"));
  m_GibsonLanniPSFSettingsTableModel->
    setHeaderData(column++, Qt::Horizontal, tr("Value"));
  m_GibsonLanniPSFSettingsTableModel->
    setHeaderData(column++, Qt::Horizontal, tr("Units"));
  m_GibsonLanniPSFSettingsTableModel->
    setHeaderData(column++, Qt::Horizontal, tr("Optimize?"));

  QStandardItem* psfPropertyItems[26];
  QStandardItem* unitItems[26];

  int item = 0;
  psfPropertyItems[item] = new QStandardItem(tr("X pixel size"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Y pixel size"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Z slice spacing"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("CCD border X"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("CCD border Y"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Bead radius"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Bead center X"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Bead center Y"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Bead center Z"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Shear X"));
  unitItems[item++] = new QStandardItem(tr("nanometers in X vs nanometers in Z"));

  psfPropertyItems[item] = new QStandardItem(tr("Shear Y"));
  unitItems[item++] = new QStandardItem(tr("nanometers in Y vs nanometers in Z"));

  psfPropertyItems[item] = new QStandardItem(tr("Emission Wavelength"));
  unitItems[item++] = new QStandardItem(tr("nanometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Numerical Aperture"));
  unitItems[item++] = new QStandardItem(tr("unitless"));

  psfPropertyItems[item] = new QStandardItem(tr("Magnification"));
  unitItems[item++] = new QStandardItem(tr("unitless"));

  psfPropertyItems[item] = new QStandardItem(tr("Design Cover Slip Refractive Index"));
  unitItems[item++] = new QStandardItem(tr("unitless"));

  psfPropertyItems[item] = new QStandardItem(tr("Actual Cover Slip Refractive Index"));
  unitItems[item++] = new QStandardItem(tr("unitless"));

  psfPropertyItems[item] = new QStandardItem(tr("Design Cover Slip Thickness"));
  unitItems[item++] = new QStandardItem(tr("micrometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Actual Cover Slip Thickness"));
  unitItems[item++] = new QStandardItem(tr("micrometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Design Immersion Oil Refractive Index"));
  unitItems[item++] = new QStandardItem(tr("unitless"));

  psfPropertyItems[item] = new QStandardItem(tr("Actual Immersion Oil Refractive Index"));
  unitItems[item++] = new QStandardItem(tr("unitless"));

  psfPropertyItems[item] = new QStandardItem(tr("Design Immersion Oil Thickness"));
  unitItems[item++] = new QStandardItem(tr("micrometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Design Specimen Layer Refractive Index"));
  unitItems[item++] = new QStandardItem(tr("unitless"));

  psfPropertyItems[item] = new QStandardItem(tr("Actual Specimen Layer Refractive Index"));
  unitItems[item++] = new QStandardItem(tr("unitless"));

  psfPropertyItems[item] = new QStandardItem(tr("Actual Point Source Depth in Specimen Layer"));
  unitItems[item++] = new QStandardItem(tr("micrometers"));

  psfPropertyItems[item] = new QStandardItem(tr("Design Distance from Back Focal Plane to Detector"));
  unitItems[item++] = new QStandardItem(tr("millimeters"));

  psfPropertyItems[item] = new QStandardItem(tr("Actual Distance from Back Focal Plane to Detector"));
  unitItems[item++] = new QStandardItem(tr("millimeters"));

  for (unsigned int i = 0; i < sizeof(psfPropertyItems) / sizeof(QStandardItem*); i++) {

    // Property column
    psfPropertyItems[i]->setEditable(false);
    m_GibsonLanniPSFSettingsTableModel->setItem(i, 0, psfPropertyItems[i]);

    // Value column
    QStandardItem* item = new QStandardItem(tr(""));
    m_GibsonLanniPSFSettingsTableModel->setItem(i, 1, item);

    // Units column
    unitItems[i]->setEditable(false);
    m_GibsonLanniPSFSettingsTableModel->setItem(i, 2, unitItems[i]);

    // Optimize checkbox column
    item = new QStandardItem();
    item->setCheckable(true);
    item->setEditable(false);
    m_GibsonLanniPSFSettingsTableModel->setItem(i, 3, item);
  }

  gui->psfSettingsTableView->setModel(m_GibsonLanniPSFSettingsTableModel);
  gui->psfSettingsTableView->setColumnWidth(0, 300);

  // Refresh the UI
  RefreshUI();

  // Reset camera
  m_Renderer->ResetCamera();
  
  // Render
  gui->qvtkWidget->GetRenderWindow()->Render();
}


// Destructor
VisualPSFOptimizer
::~VisualPSFOptimizer() {
  delete m_DataModel;
  delete m_Visualization;
}


void
VisualPSFOptimizer
::on_actionNewImage_triggered() {
  // Create dialog box and show it.
  if (m_NewFileDialog.exec() == QDialog::Accepted) {

    // Read out settings from the interface and create a new image.
    int xSize = m_NewFileDialogUI.xSizeEdit->text().toInt();
    int ySize = m_NewFileDialogUI.ySizeEdit->text().toInt();
    int zSize = m_NewFileDialogUI.zSizeEdit->text().toInt();
    float xSpacing = m_NewFileDialogUI.xSpacingEdit->text().toFloat();
    float ySpacing = m_NewFileDialogUI.ySpacingEdit->text().toFloat();
    float zSpacing = m_NewFileDialogUI.zSpacingEdit->text().toFloat();
    CreateFile(xSize, ySize, zSize, xSpacing, ySpacing, zSpacing);
  }

}


void
VisualPSFOptimizer
::on_actionOpenImage_triggered() {

  // Locate file.
  QString fileName = QFileDialog::getOpenFileName(this, "Open Image Data", "", "TIF Images (*.tif);;VTK Images (*.vtk);;LSM Images (*.lsm)");

  // Now read the file
  if (fileName == "") {
    return;
  }

  // Should probably report if opening the image failed.
  OpenFile(fileName.toStdString());
}


void
VisualPSFOptimizer
::CreateFile(int xSize, int ySize, int zSize,
             float xSpacing, float ySpacing, float zSpacing) {
  // Create new image in data model
  m_DataModel->CreateImageFile(xSize, ySize, zSize,
                               xSpacing, ySpacing, zSpacing);

  // Set status bar.
  QString imageInfo("Created new image.");
  gui->statusbar->showMessage(imageInfo);

  SetupRenderer();

  // Enable/disable appropriate GUI widgets
  gui->measuredPSFRadioButton->setEnabled(false);
  gui->calculatedPSFRadioButton->setEnabled(true);
  gui->calculatedBSFRadioButton->setEnabled(true);

  gui->estimatePSFCenterButton->setEnabled(false);
  gui->optimizePSFParametersButton->setEnabled(false);

  gui->calculatedPSFRadioButton->click();
}


void
VisualPSFOptimizer
::OpenFile(std::string fileName) {
  m_DataModel->LoadImageFile(fileName);
  
  // Set status bar with info about the file.
  QString imageInfo("Loaded image '");
  imageInfo.append(fileName.c_str()); imageInfo.append("'.");
  gui->statusbar->showMessage(imageInfo);

  SetupRenderer();

  gui->measuredPSFRadioButton->setEnabled(true);
  gui->calculatedPSFRadioButton->setEnabled(true);
  gui->calculatedBSFRadioButton->setEnabled(true);

  gui->estimatePSFCenterButton->setEnabled(true);
  gui->optimizePSFParametersButton->setEnabled(true);

  gui->measuredPSFRadioButton->click();
}


void
VisualPSFOptimizer
::SetupRenderer() {
  // Set up m_Visualization pipeline.
  m_Visualization->SetImageInputConnection(m_DataModel->GetMeasuredImageOutputPort());

  // Should clamp this to the valid range for the newly-opened file.
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);

  m_Visualization->SetXPlane(CLAMP(gui->xPlaneEdit->text().toInt()-1,0,dims[0]-1));
  m_Visualization->SetYPlane(CLAMP(gui->yPlaneEdit->text().toInt()-1,0,dims[1]-1));
  m_Visualization->SetZPlane(CLAMP(gui->zPlaneEdit->text().toInt()-1,0,dims[2]-1));

  gui->measuredPSFRadioButton->setEnabled(true);
  gui->calculatedPSFRadioButton->setEnabled(true);
  gui->calculatedBSFRadioButton->setEnabled(true);

  // Refresh the UI
  RefreshUI();

  // Reset camera
  m_Renderer->ResetCamera();
  
  // Render
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_actionSavePSFImage_triggered() {

  // Locate file.
  QString fileName = QFileDialog::getSaveFileName(this, "Save PSF Image Data", "", "TIF Images (*.tif);;VTK Images (*.vtk);;LSM Images (*.lsm)");

  // Now read the file
  if (fileName == "") {
    return;
  }

  m_DataModel->SavePSFImageFile(fileName.toStdString());
  
}


void
VisualPSFOptimizer
::on_actionSaveBSFImage_triggered() {

  // Locate file.
  QString fileName = QFileDialog::getSaveFileName(this, "Save BSF Image Data", "", "TIF Images (*.tif);;VTK Images (*.vtk);;LSM Images (*.lsm)");

  // Now read the file
  if (fileName == "") {
    return;
  }

  m_DataModel->SaveBSFImageFile(fileName.toStdString());
  
}


void
VisualPSFOptimizer
::on_actionLoadSession_triggered() {
  // Locate file.
  QString fileName = QFileDialog::getOpenFileName(this, "Load Settings", "", "VisualPSFOptimizer Settings Files (*.vpo);;All Files (*)");

  if (fileName == "") {
    return;
  }

  Configuration config;
  config.Parse(fileName.toStdString());

  OpenFile(config.GetValue("FileInfo", "FileName"));
  m_DataModel->SetConfiguration(config);

  RefreshUI();
  on_applyButton_clicked();
}


void
VisualPSFOptimizer
::on_actionSaveSession_triggered() {
  // Locate file.
  QString fileName = QFileDialog::getSaveFileName(this, "Save Settings", "", "VisualPSFOptimizer Settings Files (*.vpo);;All Files (*)");

  if (fileName == "") {
    return;
  }

  Configuration config;
  m_DataModel->GetConfiguration(config);
  std::ofstream os(fileName.toStdString().c_str());
  config.Write(os);
}


void
VisualPSFOptimizer
::on_actionExit_triggered() {
  // Ask if user really wants to quit.
  QMessageBox messageBox;
  messageBox.setText("Do you really want to exit?");
  messageBox.setInformativeText("If you exit now, all unsaved settings will be lost.");
  messageBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
  messageBox.setDefaultButton(QMessageBox::Cancel);
  int selected = messageBox.exec();

  if (selected == QMessageBox::Ok) {
    writeProgramSettings();
    qApp->exit();
  }

}


void
VisualPSFOptimizer
::on_actionAboutApplication_triggered() {
  QString version = QString().sprintf("%d.%d.%d", 
				      VisualPSFOptimizer_MAJOR_NUMBER,
				      VisualPSFOptimizer_MINOR_NUMBER,
				      VisualPSFOptimizer_REVISION_NUMBER);
  QChar copyright(169);
  QString title = QString("About VisualPSFOptimizer ").append(version);
  QString text  = QString("VisualPSFOptimizer ").append(version).append("\n");
  text.append(copyright).append(" 2010, UNC CISMM\n\n");
  text.append("Developed by:\n");
  text.append("Cory Quammen");
  QMessageBox::about(this, title, text);
}


void
VisualPSFOptimizer
::on_measuredPSFRadioButton_clicked(bool state) {
  SetDisplayedImageToMeasuredPSF();
}


void
VisualPSFOptimizer
::on_calculatedPSFRadioButton_clicked(bool state) {
  SetDisplayedImageToCalculatedPSF();
}


void
VisualPSFOptimizer
::on_calculatedBSFRadioButton_clicked(bool state) {
  SetDisplayedImageToCalculatedBSF();
}


void
VisualPSFOptimizer
::SetDisplayedImageToMeasuredPSF() {
  m_DisplayedImage = MEASURED_PSF_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetMeasuredImageOutputPort());
  SetMapsToBlackValueFromSliderPosition(gui->mapsToBlackSlider->sliderPosition());
  SetMapsToWhiteValueFromSliderPosition(gui->mapsToWhiteSlider->sliderPosition());

  RefreshUI();
}


void
VisualPSFOptimizer
::SetDisplayedImageToCalculatedPSF() {
  m_DisplayedImage = CALCULATED_PSF_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetPSFImageOutputPort());
  SetMapsToBlackValueFromSliderPosition(gui->mapsToBlackSlider->sliderPosition());
  SetMapsToWhiteValueFromSliderPosition(gui->mapsToWhiteSlider->sliderPosition());

  RefreshUI();
}


void
VisualPSFOptimizer
::SetDisplayedImageToCalculatedBSF() {
  m_DisplayedImage = CALCULATED_BSF_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetBSFImageOutputPort());
  SetMapsToBlackValueFromSliderPosition(gui->mapsToBlackSlider->sliderPosition());
  SetMapsToWhiteValueFromSliderPosition(gui->mapsToWhiteSlider->sliderPosition());

  RefreshUI();
}


void
VisualPSFOptimizer
::on_showXPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowXPlane(show);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_xPlaneSlider_valueChanged(int plane) {
  m_Visualization->SetXPlane(plane-1);
  gui->xPlaneEdit->setText(QString().sprintf("%d", plane));
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_xPlaneEdit_textEdited(QString text) {
  int value = text.toInt();
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);
  int plane = value-1;
  if (plane >= 0 && plane < dims[0]) {
    gui->xPlaneSlider->setValue(value);
    m_Visualization->SetXPlane(plane);
    gui->qvtkWidget->GetRenderWindow()->Render();
  }
}

  
void
VisualPSFOptimizer
::on_showYPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowYPlane(show);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_yPlaneSlider_valueChanged(int plane) {
  m_Visualization->SetYPlane(plane-1);
  gui->yPlaneEdit->setText(QString().sprintf("%d", plane));
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_yPlaneEdit_textEdited(QString text) {
  int value = text.toInt();
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);
  int plane = value-1;
  if (plane >= 0 && plane < dims[1]) {
    gui->yPlaneSlider->setValue(value);
    m_Visualization->SetYPlane(plane);
    gui->qvtkWidget->GetRenderWindow()->Render();
  }
}


void
VisualPSFOptimizer
::on_showZPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowZPlane(show);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_zPlaneSlider_valueChanged(int plane) {
  m_Visualization->SetZPlane(plane-1);
  gui->zPlaneEdit->setText(QString().sprintf("%d", plane));
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_zPlaneEdit_textEdited(QString text) {
  int value = text.toInt();
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);
  int plane = value-1;
  if (plane >= 0 && plane < dims[2]) {
    gui->zPlaneSlider->setValue(value);
    m_Visualization->SetZPlane(plane);
    gui->qvtkWidget->GetRenderWindow()->Render();
  }
}


void
VisualPSFOptimizer
::on_mapsToBlackSlider_valueChanged(int value) {
  SetMapsToBlackValueFromSliderPosition(value);
  gui->qvtkWidget->GetRenderWindow()->Render();  
}


void
VisualPSFOptimizer
::on_mapsToWhiteSlider_valueChanged(int value) {
  SetMapsToWhiteValueFromSliderPosition(value);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_showDataOutlineCheckBox_toggled(bool show) {
  m_Visualization->SetShowOutline(show);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_xPlusButton_clicked() {
  m_Visualization->SetViewToXPlus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_xMinusButton_clicked() {
  m_Visualization->SetViewToXMinus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_yPlusButton_clicked() {
  m_Visualization->SetViewToYPlus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_yMinusButton_clicked() {
  m_Visualization->SetViewToYMinus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_zPlusButton_clicked() {
  m_Visualization->SetViewToZPlus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_zMinusButton_clicked() {
  m_Visualization->SetViewToZMinus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_applyButton_clicked() {
  m_Dirty = false;

  float value = 0.0f;
  QStandardItem* item = NULL;
  int COLUMN = 1;
  int itemRow = 0;

  double vec3[3];
  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  vec3[0] = item->text().toDouble();

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  vec3[1] = item->text().toDouble();

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  vec3[2] = item->text().toDouble();
  m_DataModel->SetMeasuredImageVoxelSpacing(vec3);
  m_DataModel->SetPSFImageVoxelSpacing(vec3);
  m_DataModel->SetBSFImageVoxelSpacing(vec3);

  double ccdBorder[2];
  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  ccdBorder[0] = item->text().toDouble();
  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  ccdBorder[1] = item->text().toDouble();
  m_DataModel->SetCCDBorderWidth(ccdBorder);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  double beadRadius = item->text().toDouble();
  m_DataModel->SetBeadRadius(beadRadius);

  double pointCenter[3];
  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  pointCenter[0] = item->text().toDouble();
  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  pointCenter[1] = item->text().toDouble();
  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  pointCenter[2] = item->text().toDouble();
  m_DataModel->SetPSFPointCenter(pointCenter);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  m_DataModel->SetShearX(item->text().toFloat());

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  m_DataModel->SetShearY(item->text().toFloat());

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLEmissionWavelength(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLNumericalAperture(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLMagnification(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLDesignCoverSlipRefractiveIndex(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLActualCoverSlipRefractiveIndex(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLDesignCoverSlipThickness(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLActualCoverSlipThickness(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLDesignImmersionOilRefractiveIndex(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLActualImmersionOilRefractiveIndex(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLDesignImmersionOilThickness(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLDesignSpecimenLayerRefractiveIndex(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLActualSpecimenLayerRefractiveIndex(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLActualPointSourceDepthInSpecimenLayer(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLDesignDistanceFromBackFocalPlaneToDetector(value);

  item = m_GibsonLanniPSFSettingsTableModel->item(itemRow++, COLUMN);
  value = item->text().toFloat();
  m_DataModel->SetGLActualDistanceFromBackFocalPlaneToDetector(value);

  // Now update which parameters should be optimized
  for (int i = 0; i < m_GibsonLanniPSFSettingsTableModel->rowCount(); i++) {
    item = m_GibsonLanniPSFSettingsTableModel->item(i, 3);
    m_DataModel->SetGLParameterEnabled(i, item->checkState() == Qt::Checked);
  }
    
  // General image parameters
  COLUMN = 1;
  itemRow = 5;

  // Set up origin so that (0, 0, 0) is centered in the image volume.
  int dimensions[3];
  double spacing[3], origin[3];
  m_DataModel->GetPSFImageDimensions(dimensions);
  m_DataModel->GetPSFImageVoxelSpacing(spacing);
  for (int i = 0; i < 3; i++) {
    origin[i] = -0.5*static_cast<double>(dimensions[i]-1)*spacing[i];
  }

  m_DataModel->SetMeasuredImageOrigin(origin);
  m_DataModel->SetPSFImageOrigin(origin);
  m_DataModel->SetBSFImageOrigin(origin);

  // Now update
  if (gui->calculatedPSFRadioButton->isChecked()) {
    m_DataModel->UpdateGibsonLanniPSFImage();
  } else if (gui->calculatedBSFRadioButton->isChecked()) {
    m_DataModel->UpdateGibsonLanniBSFImage();
  }
  
  SetMapsToBlackValueFromSliderPosition(gui->mapsToBlackSlider->sliderPosition());
  SetMapsToWhiteValueFromSliderPosition(gui->mapsToWhiteSlider->sliderPosition());

  RefreshUI();
  
}


void
VisualPSFOptimizer
::on_estimatePSFCenterButton_clicked() {
  DataModel::Float3DPointType center = m_DataModel->GetMeasuredImageDataMaximumCoordinates();

  char decimalFormat[] = "%.3f";
  int item = 6;
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, center[0]));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, center[1]));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, center[2]));

  Sully();
}


void
VisualPSFOptimizer
::on_optimizePSFParametersButton_clicked() {
  m_DataModel->Optimize();

  RefreshUI();
}


void
VisualPSFOptimizer
::handle_imageInformationTableModel_dataChanged(const QModelIndex& topLeft,
    const QModelIndex& bottomRight) {

  if (topLeft != bottomRight) {
    return;
  }

  Sully();
}


void
VisualPSFOptimizer
::Sully() {
  m_Dirty = true;
}


void
VisualPSFOptimizer
::RefreshUI() {

  ///////////////// Update window title /////////////////
  QFileInfo fileInfo(m_DataModel->GetMeasuredImageFileName().c_str());
  QString windowTitle("VisualPSFOptimizer");
  if (fileInfo.fileName() != "")
    windowTitle.append(tr(" - '").append(fileInfo.fileName()).append("'"));
  setWindowTitle(windowTitle);
  
  const char* decimalFormat = "%.3f";
  const char* intFormat = "%d";

  gui->showDataOutlineCheckBox->setChecked(m_Visualization->GetShowOutline());
  
  ///////////////// Image planes stuff /////////////////
  int dim[3];
  m_DataModel->GetMeasuredImageDimensions(dim);
  
  gui->showXPlaneCheckBox->setChecked(m_Visualization->GetShowXPlane());
  gui->xPlaneSlider->setMinimum(1);
  gui->xPlaneSlider->setMaximum(dim[0]);
  gui->xPlaneEdit->setText(QString().sprintf(intFormat, m_Visualization->GetXPlane()+1));
  
  gui->showYPlaneCheckBox->setChecked(m_Visualization->GetShowYPlane());
  gui->yPlaneSlider->setMinimum(1);
  gui->yPlaneSlider->setMaximum(dim[1]);
  gui->yPlaneEdit->setText(QString().sprintf(intFormat, m_Visualization->GetYPlane()+1));
  
  gui->showZPlaneCheckBox->setChecked(m_Visualization->GetShowZPlane());
  gui->zPlaneSlider->setMinimum(1);
  gui->zPlaneSlider->setMaximum(dim[2]);
  gui->zPlaneEdit->setText(QString().sprintf(intFormat, m_Visualization->GetZPlane()+1));
  
  ///////////////// Image information update /////////////////
  int item = 0;
  QString dataMin = QString().sprintf(decimalFormat, GetDisplayedImageDataMinimum());
  m_ImageInformationTableModel->item(item++, 1)->setText(dataMin);
  QString dataMax = QString().sprintf(decimalFormat, GetDisplayedImageDataMaximum());
  m_ImageInformationTableModel->item(item++, 1)->setText(dataMax);
  
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);
  QString xDim = QString().sprintf(intFormat, dims[0]);
  m_ImageInformationTableModel->item(item++, 1)->setText(xDim);
  QString yDim = QString().sprintf(intFormat, dims[1]);
  m_ImageInformationTableModel->item(item++, 1)->setText(yDim);
  QString zDim = QString().sprintf(intFormat, dims[2]);
  m_ImageInformationTableModel->item(item++, 1)->setText(zDim);

  ///////////////// PSF settings update /////////////////
  item = 0;

  double spacing[3];
  m_DataModel->GetMeasuredImageVoxelSpacing(spacing);
  QString xSpacing = QString().sprintf(decimalFormat, spacing[0]);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(xSpacing);

  QString ySpacing = QString().sprintf(decimalFormat, spacing[1]);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(ySpacing);

  QString zSpacing = QString().sprintf(decimalFormat, spacing[2]);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(zSpacing);

  double ccdBorder[2];
  m_DataModel->GetCCDBorderWidth(ccdBorder);
  QString ccdBorderX = QString().sprintf(decimalFormat, ccdBorder[0]);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(ccdBorderX);

  QString ccdBorderY = QString().sprintf(decimalFormat, ccdBorder[1]);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(ccdBorderY);

  QString beadRadius = QString().sprintf(decimalFormat, m_DataModel->GetBeadRadius());
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(beadRadius);

  double pointCenter[3];
  m_DataModel->GetPSFPointCenter(pointCenter);
  QString xPointCenter = QString().sprintf(decimalFormat, pointCenter[0]);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(xPointCenter);

  QString yPointCenter = QString().sprintf(decimalFormat, pointCenter[1]);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(yPointCenter);

  QString zPointCenter = QString().sprintf(decimalFormat, pointCenter[2]);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(zPointCenter);

  float xShear = m_DataModel->GetShearX();
  QString xShearStr = QString().sprintf(decimalFormat, xShear);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(xShearStr);

  float yShear = m_DataModel->GetShearY();
  QString yShearStr = QString().sprintf(decimalFormat, yShear);
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->setText(yShearStr);

  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLEmissionWavelength()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLNumericalAperture()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLMagnification()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLDesignCoverSlipRefractiveIndex()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLActualCoverSlipRefractiveIndex()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLDesignCoverSlipThickness()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLActualCoverSlipThickness()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLDesignImmersionOilRefractiveIndex()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLActualImmersionOilRefractiveIndex()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLDesignImmersionOilThickness()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLDesignSpecimenLayerRefractiveIndex()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLActualSpecimenLayerRefractiveIndex()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLActualPointSourceDepthInSpecimenLayer()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLDesignDistanceFromBackFocalPlaneToDetector()));
  m_GibsonLanniPSFSettingsTableModel->item(item++, 1)->
    setText(QString().sprintf(decimalFormat, m_DataModel->GetGLActualDistanceFromBackFocalPlaneToDetector()));
  
  for (int i = 0; i < m_GibsonLanniPSFSettingsTableModel->rowCount(); i++) {
    m_GibsonLanniPSFSettingsTableModel->item(i, 3)->setCheckState(m_DataModel->GetGLParameterEnabled(i) ? Qt::Checked : Qt::Unchecked);
  }

  ///////////////// Update visualization stuff /////////////////
  m_Renderer->RemoveAllViewProps();

  if (m_DataModel->GetMeasuredImageData()) {
    m_Visualization->AddToRenderer();
  }

  gui->qvtkWidget->GetRenderWindow()->Render();
}


double
VisualPSFOptimizer
::GetDisplayedImageDataMinimum() {
  if (m_DisplayedImage == MEASURED_PSF_IMAGE) {
    return m_DataModel->GetMeasuredImageDataMinimum();
  } else if (m_DisplayedImage == CALCULATED_PSF_IMAGE) {
    return m_DataModel->GetPSFImageDataMinimum();
  } else if (m_DisplayedImage == CALCULATED_BSF_IMAGE) {
    return m_DataModel->GetBSFImageDataMinimum();
  }
  return 0.0;
}


double
VisualPSFOptimizer
::GetDisplayedImageDataMaximum() {
  if (m_DisplayedImage == MEASURED_PSF_IMAGE) {
    return m_DataModel->GetMeasuredImageDataMaximum();
  } else if (m_DisplayedImage == CALCULATED_PSF_IMAGE) {
    return m_DataModel->GetPSFImageDataMaximum();
  } else if (m_DisplayedImage == CALCULATED_BSF_IMAGE) {
    return m_DataModel->GetBSFImageDataMaximum();
  }
  return 0.0;
}


void
VisualPSFOptimizer
::SetMapsToBlackValueFromSliderPosition(int position) {
  double dataMin = GetDisplayedImageDataMinimum();
  double dataMax = GetDisplayedImageDataMaximum();
  double dd = dataMax - dataMin;
  double sliderMax = static_cast<double>(gui->mapsToBlackSlider->maximum());
  double dvalue = static_cast<double>(position);
  double mapped = (dvalue/sliderMax) * dd + dataMin;
  m_Visualization->SetImagePlanesBlackValue(mapped);
}


void
VisualPSFOptimizer
::SetMapsToWhiteValueFromSliderPosition(int position) {
  double dataMin = GetDisplayedImageDataMinimum();
  double dataMax = GetDisplayedImageDataMaximum();
  double dd = dataMax - dataMin;
  double sliderMax = static_cast<double>(gui->mapsToWhiteSlider->maximum());
  double dvalue = static_cast<double>(position);
  double mapped = (dvalue/sliderMax) * dd + dataMin;
  m_Visualization->SetImagePlanesWhiteValue(mapped);
}


void
VisualPSFOptimizer
::writeProgramSettings() {
  QSettings settings;

  settings.beginGroup("MainWindow");
  settings.setValue("size", size());
  settings.setValue("pos", pos());
  settings.endGroup();

  QList<QDockWidget*> widgets = findChildren<QDockWidget*>();
  QListIterator<QDockWidget*> iterator(widgets);
  while (iterator.hasNext()) {
    QDockWidget* dockWidget = iterator.next();
    settings.beginGroup(dockWidget->objectName());
    settings.setValue("size", dockWidget->size());
    settings.setValue("pos", dockWidget->pos());
    settings.setValue("visible", dockWidget->isVisible());
    settings.setValue("floating", dockWidget->isFloating());
    settings.setValue("dockArea", dockWidgetArea(dockWidget));
    settings.endGroup();
  }

}


void
VisualPSFOptimizer
::readProgramSettings() {
  QSettings settings;

  settings.beginGroup("MainWindow");
  resize(settings.value("size", QSize(1000, 743)).toSize());
  move(settings.value("pos", QPoint(0, 20)).toPoint());
  settings.endGroup();

  QList<QDockWidget*> widgets = findChildren<QDockWidget*>();
  QListIterator<QDockWidget*> iterator(widgets);
  while (iterator.hasNext()) {
    QDockWidget* dockWidget = iterator.next();
    settings.beginGroup(dockWidget->objectName());
    dockWidget->resize(settings.value("size", QSize(340, 200)).toSize());
    dockWidget->move(settings.value("pos", QPoint(0, 0)).toPoint());
    dockWidget->setVisible(settings.value("visible", true).toBool());
    dockWidget->setFloating(settings.value("floating", false).toBool());
    addDockWidget(static_cast<Qt::DockWidgetArea>(settings.value("dockArea", Qt::LeftDockWidgetArea).toUInt()), dockWidget);
    settings.endGroup();
  }

}
