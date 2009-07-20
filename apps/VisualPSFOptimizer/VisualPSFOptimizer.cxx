#include "VisualPSFOptimizer.h"
#include "Version.h"

#include <qmessagebox.h>

#if defined(_WIN32) // Turn off deprecation warnings in Visual Studio
#pragma warning( disable : 4996 )
#endif

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


// Constructor
VisualPSFOptimizer
::VisualPSFOptimizer(QWidget* p)
 : QMainWindow(p) {
  setupUi(this);
  
  // QT/VTK interaction
  m_Renderer = vtkRenderer::New();
  qvtkWidget->GetRenderWindow()->AddRenderer(m_Renderer);
  
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
  
  // Set up error dialog box.
  m_ErrorDialog.setModal(true);
  
  // Create and populate image information table model.
  int LEFT_COLUMN = 0;
  int RIGHT_COLUMN = 1;
  m_ImageInformationTableModel = new QStandardItemModel(8, 2, this);
  m_ImageInformationTableModel->setHeaderData(LEFT_COLUMN,  Qt::Horizontal, tr("Property"));
  m_ImageInformationTableModel->setHeaderData(RIGHT_COLUMN, Qt::Horizontal, tr("Value"));
  
  QStandardItem* labelItems[11];
  labelItems[ 0] = new QStandardItem(tr("Intensity minimum"));
  labelItems[ 1] = new QStandardItem(tr("Intensity maximum"));
  labelItems[ 2] = new QStandardItem(tr("X dimension (pixels)"));
  labelItems[ 3] = new QStandardItem(tr("Y dimension (pixels)"));
  labelItems[ 4] = new QStandardItem(tr("Z dimension (slices)"));
  labelItems[ 5] = new QStandardItem(tr("X pixel size (nm)"));
  labelItems[ 6] = new QStandardItem(tr("Y pixel size (nm)"));
  labelItems[ 7] = new QStandardItem(tr("Z slice spacing (nm)"));
  labelItems[ 8] = new QStandardItem(tr("X origin (nm)"));
  labelItems[ 9] = new QStandardItem(tr("Y origin (nm)"));
  labelItems[10] = new QStandardItem(tr("Z origin (nm)"));
 
  for (unsigned int i = 0; i < sizeof(labelItems) / sizeof(QStandardItem*); i++) {
    labelItems[i]->setEditable(false);
    m_ImageInformationTableModel->setItem(i, LEFT_COLUMN, labelItems[i]);

    QStandardItem* item = new QStandardItem(tr(""));

    // Allow editing of voxel spacing and origin.
    if (i < 5)
      item->setEditable(false);
    m_ImageInformationTableModel->setItem(i, RIGHT_COLUMN, item);
  }
  imageDataView->setModel(m_ImageInformationTableModel);
  
  connect(m_ImageInformationTableModel, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(handle_imageInformationTableModel_dataChanged(const QModelIndex&, const QModelIndex&)));
  
  int LEFT_COLUMN_WIDTH = 160;
  imageDataView->setColumnWidth(LEFT_COLUMN, LEFT_COLUMN_WIDTH);

  // Create and popuplate Gibson-Lanni PSF settings table model.
  m_GibsonLanniPSFSettingsTableModel = new QStandardItemModel(16, 3, this);
  m_GibsonLanniPSFSettingsTableModel->setHeaderData(0, Qt::Horizontal, tr("Property"));
  m_GibsonLanniPSFSettingsTableModel->setHeaderData(1, Qt::Horizontal, tr("Value"));
  m_GibsonLanniPSFSettingsTableModel->setHeaderData(2, Qt::Horizontal, tr("Units"));

  QStandardItem* psfPropertyItems[15];
  QStandardItem* unitItems[15];

  int item = 0;
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

  psfPropertyItems[item] = new QStandardItem(tr("Design Point Source Distance from Back Focal Plane"));
  unitItems[item++] = new QStandardItem(tr("millimeters"));

  psfPropertyItems[item] = new QStandardItem(tr("Actual Point Source Distance from Back Focal Plane"));
  unitItems[item++] = new QStandardItem(tr("millimeters"));

  for (unsigned int i = 0; i < sizeof(psfPropertyItems) / sizeof(QStandardItem*); i++) {

    // Left column
    psfPropertyItems[i]->setEditable(false);
    m_GibsonLanniPSFSettingsTableModel->setItem(i, 0, psfPropertyItems[i]);

    // Middle column
    QStandardItem* item = new QStandardItem(tr(""));
    m_GibsonLanniPSFSettingsTableModel->setItem(i, 1, item);

    // Right column
    unitItems[i]->setEditable(false);
    m_GibsonLanniPSFSettingsTableModel->setItem(i, 2, unitItems[i]);
  }
  psfSettingsTableView->setModel(m_GibsonLanniPSFSettingsTableModel);
  psfSettingsTableView->setColumnWidth(0, 300);

  connect(m_GibsonLanniPSFSettingsTableModel, SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)), this, SLOT(handle_GibsonLanniPSFSettingsTableModel_dataChanged(const QModelIndex&, const QModelIndex&)));

  // Refresh the UI
  RefreshUI();

  // Reset camera
  m_Renderer->ResetCamera();
  
  // Render
  qvtkWidget->GetRenderWindow()->Render();
}


// Destructor
VisualPSFOptimizer
::~VisualPSFOptimizer() {
  delete m_DataModel;
  delete m_Visualization;
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
  
  // Set up m_Visualization pipeline.
  m_Visualization->SetImageInputConnection(m_DataModel->GetMeasuredImageOutputPort());
  m_Visualization->SetXPlane(0);
  m_Visualization->SetYPlane(0);
  m_Visualization->SetZPlane(0);

  SetDisplayedImageToMeasuredPSF();

  measuredPSFRadioButton->setEnabled(true);
  calculatedPSFRadioButton->setEnabled(true);

  // Refresh the UI
  RefreshUI();

  // Reset camera
  m_Renderer->ResetCamera();
  
  // Render
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::OpenFile(std::string fileName) {
  std::cout << "Loading file '" << fileName << "'" << std::endl;
  m_DataModel->LoadImageFile(fileName);
  
  // Set status bar with info about the file.
  QString imageInfo("Loaded image '");
  imageInfo.append(fileName.c_str()); imageInfo.append("'.");
  statusbar->showMessage(imageInfo);
}


void
VisualPSFOptimizer
::on_actionSavePSFImage_triggered() {

  // Locate file.
  QString fileName = QFileDialog::getSaveFileName(this, "Save Image Data", "", "TIF Images (*.tif);;VTK Images (*.vtk);;LSM Images (*.lsm)");

  // Now read the file
  if (fileName == "") {
    return;
  }

  m_DataModel->SavePSFImageFile(fileName.toStdString());
  
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
  text.append(copyright).append(" 2009, UNC CISMM\n\n");
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
::SetDisplayedImageToMeasuredPSF() {
  m_DisplayedImage = MEASURED_PSF_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetMeasuredImageOutputPort());
  on_mapsToBlackSlider_sliderMoved(mapsToBlackSlider->sliderPosition());
  on_mapsToWhiteSlider_sliderMoved(mapsToWhiteSlider->sliderPosition());
  m_Visualization->Update();

  RefreshUI();
}


void
VisualPSFOptimizer
::SetDisplayedImageToCalculatedPSF() {
  m_DisplayedImage = CALCULATED_PSF_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetPSFImageOutputPort());
  on_mapsToBlackSlider_sliderMoved(mapsToBlackSlider->sliderPosition());
  on_mapsToWhiteSlider_sliderMoved(mapsToWhiteSlider->sliderPosition());
  m_Visualization->Update();

  RefreshUI();
}


void
VisualPSFOptimizer
::on_showXPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowXPlane(show);
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_xPlaneSlider_sliderMoved(int plane) {
  m_Visualization->SetXPlane(plane-1);
  xPlaneEdit->setText(QString().sprintf("%d", plane));
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_xPlaneEdit_textEdited(QString text) {
  int value = text.toInt();
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);
  int plane = value-1;
  if (plane >= 0 && plane < dims[0]) {
    xPlaneSlider->setValue(value);
    m_Visualization->SetXPlane(plane);
    qvtkWidget->GetRenderWindow()->Render();
  }
}

  
void
VisualPSFOptimizer
::on_showYPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowYPlane(show);
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_yPlaneSlider_sliderMoved(int plane) {
  m_Visualization->SetYPlane(plane-1);
  yPlaneEdit->setText(QString().sprintf("%d", plane));
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_yPlaneEdit_textEdited(QString text) {
  int value = text.toInt();
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);
  int plane = value-1;
  if (plane >= 0 && plane < dims[1]) {
    yPlaneSlider->setValue(value);
    m_Visualization->SetYPlane(plane);
    qvtkWidget->GetRenderWindow()->Render();
  }
}


void
VisualPSFOptimizer
::on_showZPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowZPlane(show);
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_zPlaneSlider_sliderMoved(int plane) {
  m_Visualization->SetZPlane(plane-1);
  zPlaneEdit->setText(QString().sprintf("%d", plane));
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_zPlaneEdit_textEdited(QString text) {
  int value = text.toInt();
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);
  int plane = value-1;
  if (plane >= 0 && plane < dims[2]) {
    zPlaneSlider->setValue(value);
    m_Visualization->SetZPlane(plane);
    qvtkWidget->GetRenderWindow()->Render();
  }
}


void
VisualPSFOptimizer
::on_mapsToBlackSlider_sliderMoved(int value) {
  SetMapsToBlackValueFromSliderPosition(value);
  qvtkWidget->GetRenderWindow()->Render();  
}


void
VisualPSFOptimizer
::on_mapsToWhiteSlider_sliderMoved(int value) {
  SetMapsToWhiteValueFromSliderPosition(value);
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_showDataOutlineCheckBox_toggled(bool show) {
  m_Visualization->SetShowOutline(show);
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_xPlusButton_clicked() {
  m_Visualization->SetViewToXPlus();
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_xMinusButton_clicked() {
  m_Visualization->SetViewToXMinus();
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_yPlusButton_clicked() {
  m_Visualization->SetViewToYPlus();
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_yMinusButton_clicked() {
  m_Visualization->SetViewToYMinus();
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_zPlusButton_clicked() {
  m_Visualization->SetViewToZPlus();
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_zMinusButton_clicked() {
  m_Visualization->SetViewToZMinus();
  qvtkWidget->GetRenderWindow()->Render();
}


void
VisualPSFOptimizer
::on_applyButton_clicked() {
  float value = 0.0f;
  int itemRow = 0;
  QStandardItem* item = NULL;
  int COLUMN = 1;

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

  // General image parameters
  COLUMN = 1;
  itemRow = 5;
  double dValue = 0.0;

  item = m_ImageInformationTableModel->item(itemRow++, COLUMN);
  dValue = item->text().toDouble();
  m_DataModel->SetMeasuredImageVoxelXSpacing(dValue);
  m_DataModel->SetPSFImageVoxelXSpacing(dValue);

  item = m_ImageInformationTableModel->item(itemRow++, COLUMN);
  dValue = item->text().toDouble();
  m_DataModel->SetMeasuredImageVoxelYSpacing(dValue);
  m_DataModel->SetPSFImageVoxelYSpacing(dValue);

  item = m_ImageInformationTableModel->item(itemRow++, COLUMN);
  dValue = item->text().toDouble();
  m_DataModel->SetMeasuredImageVoxelZSpacing(dValue);
  m_DataModel->SetPSFImageVoxelZSpacing(dValue);

  double origin[3];
  item = m_ImageInformationTableModel->item(itemRow++, COLUMN);
  origin[0] = item->text().toDouble();
  item = m_ImageInformationTableModel->item(itemRow++, COLUMN);
  origin[1] = item->text().toDouble();
  item = m_ImageInformationTableModel->item(itemRow++, COLUMN);
  origin[2] = item->text().toDouble();
  m_DataModel->SetMeasuredImageOrigin(origin);
  m_DataModel->SetPSFImageOrigin(origin);

  // Now update
  m_DataModel->UpdateGibsonLanniPSFImage();
  
  SetMapsToBlackValueFromSliderPosition(mapsToBlackSlider->sliderPosition());
  SetMapsToWhiteValueFromSliderPosition(mapsToWhiteSlider->sliderPosition());

  RefreshUI();
}


void
VisualPSFOptimizer
::on_optimizePSFParametersButton_clicked() {
  m_DataModel->Optimize();
  std::cout << "Metric: " << m_DataModel->GetImageComparisonMetric()
	    << std::endl;
}


void
VisualPSFOptimizer
::on_recenterPSFOriginButton_clicked() {
  // Recenter the PSF origin parameters so that it is in the middle of the
  // image whose voxel spacings are specified in the GUI. Usage model:
  // user specifies voxel spacing and updates the PSF to recenter it
  // according to that new spacing.
  double spacing[3];
  spacing[0] = m_ImageInformationTableModel->item(5, 1)->text().toDouble();
  spacing[1] = m_ImageInformationTableModel->item(6, 1)->text().toDouble();
  spacing[2] = m_ImageInformationTableModel->item(7, 1)->text().toDouble();

  // Go to the data model for this info (it is immutable).
  int size[3];
  m_DataModel->GetPSFImageDimensions(size);
  
  double origin[3];
  for (int i = 0; i < 3; i++)
    origin[i] = -spacing[i]*static_cast<double>(size[i]-1)*0.5;

  char *decimalFormat = "%.3f";
  QString xOrigin = QString().sprintf(decimalFormat, origin[0]);
  m_ImageInformationTableModel->item(8, 1)->
    setText(QString().sprintf(decimalFormat, origin[0]));
  m_ImageInformationTableModel->item(9, 1)->
    setText(QString().sprintf(decimalFormat, origin[1]));
  m_ImageInformationTableModel->item(10, 1)->
    setText(QString().sprintf(decimalFormat, origin[2]));

}


void
VisualPSFOptimizer
::handle_imageInformationTableModel_dataChanged(const QModelIndex& topLeft,
    const QModelIndex& bottomRight) {

  if (topLeft != bottomRight) {
    return;
  }
  
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

  showDataOutlineCheckBox->setChecked(m_Visualization->GetShowOutline());
  
  ///////////////// Image planes stuff /////////////////
  int dim[3];
  m_DataModel->GetMeasuredImageDimensions(dim);
  
  showXPlaneCheckBox->setChecked(m_Visualization->GetShowXPlane());
  xPlaneSlider->setMinimum(1);
  xPlaneSlider->setMaximum(dim[0]);
  xPlaneEdit->setText(QString().sprintf(intFormat, m_Visualization->GetXPlane()+1));
  
  showYPlaneCheckBox->setChecked(m_Visualization->GetShowYPlane());
  yPlaneSlider->setMinimum(1);
  yPlaneSlider->setMaximum(dim[1]);
  yPlaneEdit->setText(QString().sprintf(intFormat, m_Visualization->GetYPlane()+1));
  
  showZPlaneCheckBox->setChecked(m_Visualization->GetShowZPlane());
  zPlaneSlider->setMinimum(1);
  zPlaneSlider->setMaximum(dim[2]);
  zPlaneEdit->setText(QString().sprintf(intFormat, m_Visualization->GetZPlane()+1));
  
  ///////////////// Image information update /////////////////
  QString dataMin = QString().sprintf(decimalFormat, GetDisplayedImageDataMinimum());
  m_ImageInformationTableModel->item(0, 1)->setText(dataMin);
  QString dataMax = QString().sprintf(decimalFormat, GetDisplayedImageDataMaximum());
  m_ImageInformationTableModel->item(1, 1)->setText(dataMax);
  
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);
  QString xDim = QString().sprintf(intFormat, dims[0]);
  m_ImageInformationTableModel->item(2, 1)->setText(xDim);
  QString yDim = QString().sprintf(intFormat, dims[1]);
  m_ImageInformationTableModel->item(3, 1)->setText(yDim);
  QString zDim = QString().sprintf(intFormat, dims[2]);
  m_ImageInformationTableModel->item(4, 1)->setText(zDim);

  double spacing[3];
  m_DataModel->GetMeasuredImageVoxelSpacing(spacing);
  QString xSpacing = QString().sprintf(decimalFormat, spacing[0]);
  m_ImageInformationTableModel->item(5, 1)->setText(xSpacing);

  QString ySpacing = QString().sprintf(decimalFormat, spacing[1]);
  m_ImageInformationTableModel->item(6, 1)->setText(ySpacing);

  QString zSpacing = QString().sprintf(decimalFormat, spacing[2]);
  m_ImageInformationTableModel->item(7, 1)->setText(zSpacing);

  double origin[3];
  m_DataModel->GetPSFImageOrigin(origin);
  QString xOrigin = QString().sprintf(decimalFormat, origin[0]);
  m_ImageInformationTableModel->item(8, 1)->setText(xOrigin);

  QString yOrigin = QString().sprintf(decimalFormat, origin[1]);
  m_ImageInformationTableModel->item(9, 1)->setText(yOrigin);

  QString zOrigin = QString().sprintf(decimalFormat, origin[2]);
  m_ImageInformationTableModel->item(10, 1)->setText(zOrigin);

  ///////////////// PSF settings update /////////////////
  int item = 0;
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
  
  ///////////////// Update visualization stuff /////////////////
  m_Renderer->RemoveAllViewProps();

  if (m_DataModel->GetMeasuredImageData()) {
    m_Visualization->AddToRenderer();
  }

  qvtkWidget->GetRenderWindow()->Render();
}


double
VisualPSFOptimizer
::GetDisplayedImageDataMinimum() {
  if (m_DisplayedImage == MEASURED_PSF_IMAGE) {
    return m_DataModel->GetMeasuredImageDataMinimum();
  } else if (m_DisplayedImage == CALCULATED_PSF_IMAGE) {
    return m_DataModel->GetPSFImageDataMinimum();
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
  }
  return 0.0;
}


void
VisualPSFOptimizer
::SetMapsToBlackValueFromSliderPosition(int position) {
  double dataMin = GetDisplayedImageDataMinimum();
  double dataMax = GetDisplayedImageDataMaximum();
  double dd = dataMax - dataMin;
  double sliderMax = static_cast<double>(mapsToBlackSlider->maximum());
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
  double sliderMax = static_cast<double>(mapsToWhiteSlider->maximum());
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
