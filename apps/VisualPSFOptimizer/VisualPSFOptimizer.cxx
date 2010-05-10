#include <VisualPSFOptimizer.h>
#include <Version.h>

#include <QMessageBox>

#if defined(_WIN32) // Turn off deprecation warnings in Visual Studio
#pragma warning( disable : 4996 )
#endif

#include <Configuration.h>
#include <DataModel.h>
#include <Visualization.h>

#include <QApplication>
#include <QFileDialog>
#include <QItemEditorFactory>
#include <QSettings>
#include <QStandardItemEditorCreator>
#include <QVariant>


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

  // Change the double item editor to QLineEdit
  QItemEditorFactory* factory = new QItemEditorFactory();
  factory->registerEditor(QVariant::Double, new QStandardItemEditorCreator<QLineEdit>());
  QItemEditorFactory::setDefaultFactory(factory);
  
  // Mark as initially clean
  m_Dirty = false;
  m_DisplayedImage = MEASURED_PSF_IMAGE;
  
  // QT/VTK interaction
  m_Renderer = vtkRenderer::New();
  m_Renderer->SetBackground(0.2, 0.2, 0.2);
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
  
  connect(m_ImageInformationTableModel, 
          SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), 
          this, SLOT(handle_imageInformationTableModel_dataChanged(const QModelIndex&, const QModelIndex&)));

  int LEFT_COLUMN_WIDTH = 160;
  gui->imageDataView->setColumnWidth(LEFT_COLUMN, LEFT_COLUMN_WIDTH);

  m_PSFPropertyTableModel = new QPointSpreadFunctionPropertyTableModel();
  m_PSFPropertyTableModel->SetDataModel(m_DataModel);

  connect(m_PSFPropertyTableModel,
          SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)),
          this, SLOT(handle_PSFPropertyTableModel_dataChanged(const QModelIndex&, const QModelIndex&)));
  
  gui->psfSettingsTableView->setModel(m_PSFPropertyTableModel);
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
    
    on_applyButton_clicked();
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

  m_PSFPropertyTableModel->InitializeSettingsCache();
  m_PSFPropertyTableModel->Refresh();
}


void
VisualPSFOptimizer
::OpenFile(std::string fileName) {

  m_DataModel->LoadImageFile(fileName);
  
  m_PSFPropertyTableModel->InitializeSettingsCache();
  m_PSFPropertyTableModel->Refresh();

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

  m_DataModel->LoadSessionFile(fileName.toStdString());

  m_PSFPropertyTableModel->InitializeSettingsCache();
  gui->measuredPSFRadioButton->setEnabled(true);
  gui->calculatedPSFRadioButton->setEnabled(true);
  gui->calculatedBSFRadioButton->setEnabled(true);

  gui->estimatePSFCenterButton->setEnabled(true);
  gui->optimizePSFParametersButton->setEnabled(true);

  gui->measuredPSFRadioButton->click();
  SetupRenderer();

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

  m_DataModel->SaveSessionFile(fileName.toStdString());
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
::on_useCustomZSlicePositions_toggled(bool use) {
  m_DataModel->SetUseCustomZCoordinates(use);
  gui->resetCustomSlicePositionsButton->setEnabled(use);

  m_PSFPropertyTableModel->Refresh();

  on_applyButton_clicked();
}


void
VisualPSFOptimizer
::on_resetCustomSlicePositionsButton_clicked() {
  
  // Reset the individual slice z positions with even increments centered
  // about z = 0
  int dims[3];
  m_DataModel->GetBSFImageDimensions(dims);
  double spacing[3];
  m_DataModel->GetBSFImageVoxelSpacing(spacing);

  double zMax = 0.5*(dims[2]-1)*spacing[2];
  for (unsigned int i = 0; i < static_cast<unsigned int>(dims[2]); i++) {
    double z = zMax - static_cast<double>(i)*spacing[2];
    m_DataModel->SetZCoordinate(i, z);
  }

  m_PSFPropertyTableModel->Refresh();

  gui->applyButton->setEnabled(true);
}


void
VisualPSFOptimizer
::on_estimatePSFCenterButton_clicked() {
  DataModel::Float3DPointType center = m_DataModel->GetMeasuredImageDataMaximumCoordinates();
  double triplet[3];

  for (int i = 0; i < 3; i++)
    triplet[i] = static_cast<double>(center[i]);
  m_DataModel->SetPSFPointCenter(triplet);
  m_DataModel->SetBSFPointCenter(triplet);

  m_PSFPropertyTableModel->InitializeSettingsCache();
  m_PSFPropertyTableModel->Refresh();

  gui->applyButton->setEnabled(true);

  Sully();
}


void
VisualPSFOptimizer
::on_applyButton_clicked() {
  m_Dirty = false;

  gui->applyButton->setEnabled(false);

  m_PSFPropertyTableModel->SaveSettingsCache();

  // Now update the image
  if (gui->calculatedPSFRadioButton->isChecked()) {
    m_DataModel->UpdateGibsonLanniPSFImage();
  } else if (gui->calculatedBSFRadioButton->isChecked()) {
    m_DataModel->UpdateGibsonLanniBSFImage();
  }
  
  SetMapsToBlackValueFromSliderPosition(gui->mapsToBlackSlider->sliderPosition());
  SetMapsToWhiteValueFromSliderPosition(gui->mapsToWhiteSlider->sliderPosition());

  RefreshUI();

  if (gui->measuredPSFRadioButton->isEnabled()) {
    double value = m_DataModel->GetImageComparisonMetricValue();
    gui->comparisonMetricLineEdit->setText(QString().sprintf("%.3f", value));
  } else {
    gui->comparisonMetricLineEdit->setText(QString("-"));
  }

}


void
VisualPSFOptimizer
::on_optimizePSFParametersButton_clicked() {
  m_DataModel->Optimize();

  Sully();
  RefreshUI();

  // Load settings from the data model.
  m_PSFPropertyTableModel->InitializeSettingsCache();
  m_PSFPropertyTableModel->Refresh();

  on_applyButton_clicked();
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
::handle_PSFPropertyTableModel_dataChanged(const QModelIndex& topLeft,
                                           const QModelIndex& bottomRight) {

  if (topLeft != bottomRight) {
    return;
  }

  gui->applyButton->setEnabled(true);

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

  ///////////////// Other widgets //////////////////////////////
  gui->useCustomZSlicePositions->
    setCheckState(m_DataModel->GetUseCustomZCoordinates() ? Qt::Checked : Qt::Unchecked);

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
