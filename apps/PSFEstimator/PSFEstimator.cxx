#include <PSFEstimator.h>
#include <Version.h>

#include <QMessageBox>

#if defined(_WIN32) // Turn off deprecation warnings in Visual Studio
#pragma warning( disable : 4996 )
#endif

#include <Configuration.h>
#include <DataModel.h>
#include <Visualization.h>

#include <QApplication>
#include <QClipboard>
#include <QCloseEvent>
#include <QDateTime>
#include <QFileDialog>
#include <QFileInfo>
#include <QItemEditorFactory>
#include <QList>
#include <QMimeData>
#include <QProcess>
#include <QRegExp>
#include <QSettings>
#include <QStandardItemEditorCreator>
#include <QVariant>
#include <QtNetwork/QHostInfo>


#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#define CLAMP(value, min, max) (value < min ? min : (value > max ? max : value))


// Constructor
PSFEstimator
::PSFEstimator(QWidget* p)
 : QMainWindow(p) {

  gui = new Ui_MainWindow();
  gui->setupUi(this);

  // Disable these widgets for now
  gui->noiseTypeWidget->setVisible(false);

  // Hide queue submission button if domain is not bass.cs.unc.edu
  QString hostName = QHostInfo::localHostName();
  QRegExp regExp("bass-comp\\d*\\.cs\\.unc\\.edu");
  if (!regExp.exactMatch(hostName)) {
    gui->submitOptimizationJobToQueueButton->setVisible(false);
  }

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
  QCoreApplication::setApplicationName("PSFEstimator");

  // Set up dialog boxes.
  m_NewFileDialogUI.setupUi(&m_NewFileDialog);

  m_ErrorDialog.setModal(true);

  // Create and populate image information table model.
  int LEFT_COLUMN  = 0;
  int RIGHT_COLUMN = 1;
  m_ImageInformationTableModel = new QStandardItemModel(5, 2, this);
  m_ImageInformationTableModel->setHeaderData(LEFT_COLUMN,  Qt::Horizontal, tr("Property"));
  m_ImageInformationTableModel->setHeaderData(RIGHT_COLUMN, Qt::Horizontal, tr("Value"));

  QStandardItem* labelItems[5];
  labelItems[ 0] = new QStandardItem(tr("Intensity Minimum"));
  labelItems[ 1] = new QStandardItem(tr("Intensity Maximum"));
  labelItems[ 2] = new QStandardItem(tr("X Dimension (pixels)"));
  labelItems[ 3] = new QStandardItem(tr("Y Dimension (pixels)"));
  labelItems[ 4] = new QStandardItem(tr("Z Dimension (slices)"));

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

  m_BSFPropertyTableModel = new QBeadSpreadFunctionPropertyTableModel();
  m_BSFPropertyTableModel->SetDataModel(m_DataModel);

  connect(m_BSFPropertyTableModel,
          SIGNAL(dataChanged(const QModelIndex&, const QModelIndex&)),
          this, SLOT(handle_BSFPropertyTableModel_dataChanged(const QModelIndex&, const QModelIndex&)));

  gui->psfSettingsTableView->setModel(m_BSFPropertyTableModel);
  gui->psfSettingsTableView->setColumnWidth(0, 300);

  // Refresh the UI
  RefreshUI();

  // Reset camera
  m_Renderer->ResetCamera();

  // Render
  on_yPlusButton_clicked();
  gui->qvtkWidget->GetRenderWindow()->Render();

  // Restore inter-session GUI settings.
  ReadProgramSettings();

}


// Destructor
PSFEstimator
::~PSFEstimator() {
  delete m_DataModel;
  delete m_Visualization;
}


void
PSFEstimator
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
    m_DataModel->CreateImageFile(xSize, ySize, zSize,
                                 xSpacing, ySpacing, zSpacing);

    m_BSFPropertyTableModel->InitializeSettingsCache();
    m_BSFPropertyTableModel->Refresh();

    SetupInterface(false);
    SetupRenderer();

    on_applyButton_clicked();
    gui->calculatedPSFRadioButton->click();

    // Set status bar.
    QString imageInfo("Created new image.");
    gui->statusbar->showMessage(imageInfo);
  }
}


void
PSFEstimator
::on_actionOpenImage_triggered() {

  // Locate file.
  QString fileName =
    QFileDialog::getOpenFileName
    (this, "Open Image Data", GetFileChooserDirectory(), "TIF Images (*.tif);;VTK Images (*.vtk);;LSM Images (*.lsm)");
  if (fileName == "") {
    return;
  }
  SaveFileChooserDirectory(fileName);

  // Should probably report if opening the image failed.
  m_DataModel->LoadImageFile(fileName.toStdString());

  m_BSFPropertyTableModel->InitializeSettingsCache();
  m_BSFPropertyTableModel->Refresh();

  SetupInterface(true);
  SetupRenderer();

  // Set status bar with info about the file.
  QString imageInfo("Loaded image '");
  imageInfo.append(fileName); imageInfo.append("'.");
  gui->statusbar->showMessage(imageInfo);
}


void
PSFEstimator
::SetupInterface(bool hasMeasuredImage) {
  // Enable the use of custom Z slices positions in any case
  gui->useCustomZSlicePositions->setEnabled(true);

  gui->measuredBSFRadioButton->setEnabled(hasMeasuredImage);
  gui->calculatedPSFRadioButton->setEnabled(true);
  gui->calculatedBSFRadioButton->setEnabled(true);
  gui->measuredMinusCalculatedBSFRadioButton->setEnabled(hasMeasuredImage);

  gui->psfModelComboBox->setEnabled(true);
  gui->noiseTypeWidget->setEnabled(true);

  gui->estimatePSFCenterButton->setEnabled(hasMeasuredImage);
  gui->optimizePSFParametersButton->setEnabled(hasMeasuredImage);
  gui->submitOptimizationJobToQueueButton->setEnabled(hasMeasuredImage);

  if (hasMeasuredImage)
    gui->measuredBSFRadioButton->click();
  else
    gui->calculatedPSFRadioButton->click();
}


void
PSFEstimator
::SetupRenderer() {
  // Set up m_Visualization pipeline.
  m_Visualization->SetImageInputConnection(m_DataModel->GetMeasuredImageOutputPort());

  // Should clamp this to the valid range for the newly-opened file.
  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);

  m_Visualization->SetXPlane(CLAMP(dims[0]/2, 0, dims[0]-1));
  m_Visualization->SetYPlane(CLAMP(dims[1]/2, 0, dims[1]-1));
  m_Visualization->SetZPlane(CLAMP(dims[2]/2, 0, dims[2]-1));
  gui->xPlaneSlider->setValue(dims[0]/2);
  gui->yPlaneSlider->setValue(dims[1]/2);
  gui->zPlaneSlider->setValue(dims[2]/2);

  // Refresh the UI
  RefreshUI();

  // Reset camera
  m_Renderer->ResetCamera();

  // Render
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_actionSavePSFImage_triggered() {

  // Locate file.
  QString fileName =
    QFileDialog::getSaveFileName
    (this, "Save PSF Image Data", GetFileChooserDirectory(),
     "TIF Images (*.tif);;VTK Images (*.vtk);;LSM Images (*.lsm)");

  // Now read the file
  if (fileName == "") {
    return;
  }
  SaveFileChooserDirectory(fileName);

  m_DataModel->SavePSFImageFile(fileName.toStdString());
}


void
PSFEstimator
::on_actionSaveBSFImage_triggered() {

  // Locate file.
  QString fileName =
    QFileDialog::getSaveFileName
    (this, "Save BSF Image Data", GetFileChooserDirectory(),
     "TIF Images (*.tif);;VTK Images (*.vtk);;LSM Images (*.lsm)");

  // Now read the file
  if (fileName == "") {
    return;
  }
  SaveFileChooserDirectory(fileName);

  m_DataModel->SaveBSFImageFile(fileName.toStdString());

}


void
PSFEstimator
::on_actionLoadSession_triggered() {
  // Locate file.
  QString fileName =
    QFileDialog::getOpenFileName
    (this, "Load Settings", GetFileChooserDirectory(),
     "PSF Estimator Settings Files (*.psfe);;All Files (*)");
  if (fileName == "") {
    return;
  }
  SaveFileChooserDirectory(fileName);

  bool success = m_DataModel->LoadSessionFile(fileName.toStdString());
  if (!success) {
    m_ErrorDialog.showMessage(tr("Error loading session file. Please check that the image file path is valid."));
    return;
  }

  m_BSFPropertyTableModel->InitializeSettingsCache();
  m_BSFPropertyTableModel->Refresh();

  SetupInterface(true);
  SetupRenderer();

  on_applyButton_clicked();
  gui->measuredBSFRadioButton->click();

  // Set status bar with info about the file.
  QString sessionFileInfo("Loaded session file '");
  sessionFileInfo.append(fileName); sessionFileInfo.append("'.");
  gui->statusbar->showMessage(sessionFileInfo);
}


void
PSFEstimator
::on_actionSaveSession_triggered() {
  // Locate file.
  QString fileName =
    QFileDialog::getSaveFileName
    (this, "Save Settings", GetFileChooserDirectory(),
     "PSF Estimator Settings Files (*.psfe);;All Files (*)");
  if (fileName == "") {
    return;
  }
  SaveFileChooserDirectory(fileName);

  if (!fileName.toLower().endsWith(tr(".psfe"))) {
    fileName.append(tr(".psfe"));
  }

  m_DataModel->SaveSessionFile(fileName.toStdString());
  QString message("Saved session file '");
  message.append(fileName); message.append("'.");
  gui->statusbar->showMessage(message);
}


void
PSFEstimator
::on_actionExit_triggered() {
  Exit();
}


void
PSFEstimator
::on_actionCopy_triggered() {
  QItemSelectionModel* selection = gui->psfSettingsTableView->selectionModel();
  QModelIndexList indices = selection->selectedIndexes();
  QStringList valueStrings;
  for (int i = 0; i < indices.length(); i++) {
    QVariant value = m_BSFPropertyTableModel->data(indices[i], Qt::DisplayRole);
    valueStrings.append(value.toString());
  }

  QClipboard* clipboard = QApplication::clipboard();
  clipboard->setText(valueStrings.join("\n"));
}


void
PSFEstimator
::on_actionPaste_triggered() {
  // Split strings pasted in. Assumes new rows are separated by newlines.
  QClipboard* clipboard = QApplication::clipboard();
  QStringList valueStrings = clipboard->text().split("\n");

  QModelIndex startIndex = gui->psfSettingsTableView->currentIndex();
  int numModelValues = m_BSFPropertyTableModel->rowCount() - startIndex.row();
  int numStringValues = valueStrings.length();

  int count = (numModelValues < numStringValues) ? numModelValues : numStringValues;
  for (int i = 0; i < count; i++) {
    QVariant value(valueStrings[i]);
    QModelIndex index = m_BSFPropertyTableModel->index(i+startIndex.row(), 1);
    m_BSFPropertyTableModel->setData(index, value);
  }

}


void
PSFEstimator
::on_actionAboutApplication_triggered() {
  QString version = QString().sprintf("%d.%d.%d",
				      PSFEstimator_MAJOR_NUMBER,
				      PSFEstimator_MINOR_NUMBER,
				      PSFEstimator_REVISION_NUMBER);
  QChar copyright(169);
  QString title = QString("About PSF Estimator ").append(version);
  QString text  = QString("PSF Estimator ").append(version).append("\n");
  text.append(copyright).append(" 2010, UNC CISMM\n\n");
  text.append("Developed by:\n");
  text.append("Cory Quammen");
  QMessageBox::about(this, title, text);
}


void
PSFEstimator
::on_measuredBSFRadioButton_clicked(bool state) {
  SetDisplayedImageToMeasuredPSF();
}


void
PSFEstimator
::on_calculatedPSFRadioButton_clicked(bool state) {
  SetDisplayedImageToCalculatedPSF();
}


void
PSFEstimator
::on_calculatedBSFRadioButton_clicked(bool state) {
  SetDisplayedImageToCalculatedBSF();
}


void
PSFEstimator
::on_measuredMinusCalculatedBSFRadioButton_clicked(bool state) {
  SetDisplayedImageToBSFDifference();
}


void
PSFEstimator
::SetDisplayedImageToMeasuredPSF() {
  m_DisplayedImage = MEASURED_PSF_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetMeasuredImageOutputPort());

  RefreshUI();
}


void
PSFEstimator
::SetDisplayedImageToCalculatedPSF() {
  m_DisplayedImage = CALCULATED_PSF_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetPSFImageOutputPort());

  RefreshUI();
}


void
PSFEstimator
::SetDisplayedImageToCalculatedBSF() {
  m_DisplayedImage = CALCULATED_BSF_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetBSFImageOutputPort());

  RefreshUI();
}


void
PSFEstimator
::SetDisplayedImageToBSFDifference() {
  m_DisplayedImage = BSF_DIFFERENCE_IMAGE;
  m_Visualization->SetImageInputConnection(m_DataModel->GetBSFDifferenceImageOutputPort());

  RefreshUI();
}


void
PSFEstimator
::on_showXPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowXPlane(show);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_xPlaneSlider_valueChanged(int plane) {
  m_Visualization->SetXPlane(plane-1);
  gui->xPlaneEdit->setText(QString().sprintf("%d", plane));
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
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
PSFEstimator
::on_showYPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowYPlane(show);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_yPlaneSlider_valueChanged(int plane) {
  m_Visualization->SetYPlane(plane-1);
  gui->yPlaneEdit->setText(QString().sprintf("%d", plane));
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
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
PSFEstimator
::on_showZPlaneCheckBox_toggled(bool show) {
  m_Visualization->SetShowZPlane(show);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_zPlaneSlider_valueChanged(int plane) {
  m_Visualization->SetZPlane(plane-1);
  gui->zPlaneEdit->setText(QString().sprintf("%d", plane));
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
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
PSFEstimator
::on_mapsToBlackSlider_valueChanged(int value) {
  SetMapsToBlackValueFromSliderPosition(value);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_mapsToWhiteSlider_valueChanged(int value) {
  SetMapsToWhiteValueFromSliderPosition(value);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_showDataOutlineCheckBox_toggled(bool show) {
  m_Visualization->SetShowOutline(show);
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_xPlusButton_clicked() {
  m_Visualization->SetViewToXPlus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_xMinusButton_clicked() {
  m_Visualization->SetViewToXMinus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_yPlusButton_clicked() {
  m_Visualization->SetViewToYPlus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_yMinusButton_clicked() {
  m_Visualization->SetViewToYMinus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_zPlusButton_clicked() {
  m_Visualization->SetViewToZPlus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_zMinusButton_clicked() {
  m_Visualization->SetViewToZMinus();
  gui->qvtkWidget->GetRenderWindow()->Render();
}


void
PSFEstimator
::on_psfModelComboBox_currentIndexChanged(int index) {
  switch (index) {

  case DataModel::GAUSSIAN_PSF:
    m_DataModel->SetPointSpreadFunctionType(DataModel::GAUSSIAN_PSF);
    break;

  case DataModel::GIBSON_LANNI_PSF:
    m_DataModel->SetPointSpreadFunctionType(DataModel::GIBSON_LANNI_PSF);
    break;

  case DataModel::HAEBERLE_PSF:
    m_DataModel->SetPointSpreadFunctionType(DataModel::HAEBERLE_PSF);
    break;

  }

  m_BSFPropertyTableModel->SetDataModel(m_DataModel);
  RefreshUI();
}


void
PSFEstimator
::on_useCustomZSlicePositions_toggled(bool use) {
  m_DataModel->SetUseCustomZCoordinates(use);
  gui->resetCustomSlicePositionsButton->setEnabled(use);

  m_BSFPropertyTableModel->Refresh();

  on_applyButton_clicked();
}


void
PSFEstimator
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

  m_BSFPropertyTableModel->Refresh();

  gui->applyButton->setEnabled(true);
}


void
PSFEstimator
::on_estimatePSFCenterButton_clicked() {
  DataModel::Float3DPointType center = m_DataModel->GetMeasuredImageDataMaximumCoordinates();
  double triplet[3];

  for (int i = 0; i < 3; i++)
    triplet[i] = static_cast<double>(center[i]);
  m_DataModel->SetPSFPointCenter(triplet);
  m_DataModel->SetBSFPointCenter(triplet);

  m_BSFPropertyTableModel->InitializeSettingsCache();
  m_BSFPropertyTableModel->Refresh();

  gui->applyButton->setEnabled(true);

  Sully();
}


void
PSFEstimator
::on_applyButton_clicked() {
  m_Dirty = false;

  gui->applyButton->setEnabled(false);

  m_BSFPropertyTableModel->SaveSettingsCache();

  // Now update the image
  if (gui->calculatedPSFRadioButton->isChecked()) {
    m_DataModel->UpdatePSFImage();
  } else if (gui->calculatedBSFRadioButton->isChecked()) {
    m_DataModel->UpdateBSFImage();
  } else if (gui->measuredMinusCalculatedBSFRadioButton->isChecked()) {
    m_DataModel->UpdateBSFDifferenceImage();
  }

  RefreshUI();

  if (gui->measuredBSFRadioButton->isEnabled()) {
    double value = m_DataModel->GetImageComparisonMetricValue();
    gui->objectiveFunctionValueEdit->setText(QString().sprintf("%.3f", value));
  } else {
    gui->objectiveFunctionValueEdit->setText(QString("-"));
  }

}


void
PSFEstimator
::on_optimizePSFParametersButton_clicked() {
  m_DataModel->Optimize();

  Sully();
  RefreshUI();

  // Load settings from the data model.
  m_BSFPropertyTableModel->InitializeSettingsCache();
  m_BSFPropertyTableModel->Refresh();

  on_applyButton_clicked();
}


void
PSFEstimator
::on_submitOptimizationJobToQueueButton_clicked() {
  // Create a working directory in the users home directory if it doesn't exist
  QDir dir;
  QString homeDirectory = dir.homePath();
  QString workingDirectory = homeDirectory;
  workingDirectory.append("/BatchPSFEstimator-Files");

  if (!dir.exists(workingDirectory))
    dir.mkdir(workingDirectory);

  QString dateTimeString = QDateTime::currentDateTime().toString(tr("MM-dd-yyyy-hh-mm-ss")); 
  QString sessionFile = workingDirectory + dir.separator() +
    tr("BatchPSFEstimator-") + dateTimeString + tr(".psfe");

  // Write the settings file to the working directory
  m_DataModel->SaveSessionFile(sessionFile.toStdString());

  // Wubmit the job to the queue
  QProcess qsub;
  QStringList arguments;
  arguments << "-N" << "BatchPSFEstimator";

  QString batchPSFEstimatorExecutable = QCoreApplication::applicationDirPath();
  batchPSFEstimatorExecutable.append("/BatchPSFEstimator");

  QStringList script;
  script <<
    "#$ -S /bin/bash\n" <<
    "#$ -o $HOME/BatchPSFEstimator-Files/$JOB_NAME.$JOB_ID\n" <<
    "#$ -j y\n" <<
    "#$ -pe smp 16\n" <<

    "# Run the job.\n" <<
    "date\n" <<
    batchPSFEstimatorExecutable << " " << sessionFile << "\n" <<
    "date\n"
    ;

  qsub.start("qsub", arguments);
  qsub.waitForStarted();

  QString stdInput = script.join(tr(""));
  std::cout << stdInput.toStdString() << std::endl;
  qsub.write(stdInput);
  qsub.closeWriteChannel();
  qsub.waitForFinished();

  QString result(qsub.readAllStandardOutput());

  if (qsub.exitCode() == 0) {
    gui->statusbar->showMessage(result);
    QMessageBox::information(this, tr("Optimization job submitted"), result);
  } else {
    QMessageBox::information(this, tr("Error"), tr("Could not submit job to queue."));
  }
}


void
PSFEstimator
::handle_imageInformationTableModel_dataChanged(const QModelIndex& topLeft,
                                                const QModelIndex& bottomRight) {

  if (topLeft != bottomRight) {
    return;
  }

  Sully();
}


void
PSFEstimator
::handle_BSFPropertyTableModel_dataChanged(const QModelIndex& topLeft,
                                           const QModelIndex& bottomRight) {

  if (topLeft != bottomRight) {
    return;
  }

  gui->applyButton->setEnabled(true);

  Sully();
}


void
PSFEstimator
::Sully() {
  m_Dirty = true;
}


void
PSFEstimator
::RefreshUI() {

  ///////////////// Update window title /////////////////
  QString version = QString().sprintf("%d.%d.%d",
                                      PSFEstimator_MAJOR_NUMBER,
                                      PSFEstimator_MINOR_NUMBER,
                                      PSFEstimator_REVISION_NUMBER);

  QString windowTitle("PSF Estimator ");
  windowTitle.append(version);

  QFileInfo fileInfo(m_DataModel->GetMeasuredImageFileName().c_str());
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
  int row = 0;
  QString dataMin = QString().sprintf(decimalFormat, GetDisplayedImageDataMinimum());
  QStandardItem* item;
  item = m_ImageInformationTableModel->item(row, 1);
  item->setText(dataMin);
  m_ImageInformationTableModel->setItem(row++, 1, item);

  QString dataMax = QString().sprintf(decimalFormat, GetDisplayedImageDataMaximum());
  item = m_ImageInformationTableModel->item(row, 1);
  item->setText(dataMax);
  m_ImageInformationTableModel->setItem(row++, 1, item);

  int dims[3];
  m_DataModel->GetMeasuredImageDimensions(dims);

  QString xDim = QString().sprintf(intFormat, dims[0]);
  item = m_ImageInformationTableModel->item(row, 1);
  item->setText(xDim);
  m_ImageInformationTableModel->setItem(row++, 1, item);

  QString yDim = QString().sprintf(intFormat, dims[1]);
  item = m_ImageInformationTableModel->item(row, 1);
  item->setText(yDim);
  m_ImageInformationTableModel->setItem(row++, 1, item);

  QString zDim = QString().sprintf(intFormat, dims[2]);
  item = m_ImageInformationTableModel->item(row, 1);
  item->setText(zDim);
  m_ImageInformationTableModel->setItem(row++, 1, item);

  gui->imageDataView->setModel(NULL);
  gui->imageDataView->setModel(m_ImageInformationTableModel);

  ///////////////// Other widgets //////////////////////////////
  gui->useCustomZSlicePositions->
    setCheckState(m_DataModel->GetUseCustomZCoordinates() ? Qt::Checked : Qt::Unchecked);

  ///////////////// Update visualization stuff /////////////////
  m_Renderer->RemoveAllViewProps();

  if (m_DataModel->GetMeasuredImageData()) {
    m_Visualization->AddToRenderer();
    SetMapsToBlackValueFromSliderPosition(gui->mapsToBlackSlider->sliderPosition());
    SetMapsToWhiteValueFromSliderPosition(gui->mapsToWhiteSlider->sliderPosition());
  }

  gui->qvtkWidget->GetRenderWindow()->Render();
}


double
PSFEstimator
::GetDisplayedImageDataMinimum() {
  double value = 0.0;
  if (m_DisplayedImage == MEASURED_PSF_IMAGE) {
    value = m_DataModel->GetMeasuredImageDataMinimum();
  } else if (m_DisplayedImage == CALCULATED_PSF_IMAGE) {
    value = m_DataModel->GetPSFImageDataMinimum();
  } else if (m_DisplayedImage == CALCULATED_BSF_IMAGE) {
    value = m_DataModel->GetBSFImageDataMinimum();
  } else if (m_DisplayedImage == BSF_DIFFERENCE_IMAGE) {
    value = m_DataModel->GetBSFDifferenceImageDataMinimum();
  }
  return value;
}


double
PSFEstimator
::GetDisplayedImageDataMaximum() {
  double value = 0.0;
  if (m_DisplayedImage == MEASURED_PSF_IMAGE) {
    value = m_DataModel->GetMeasuredImageDataMaximum();
  } else if (m_DisplayedImage == CALCULATED_PSF_IMAGE) {
    value = m_DataModel->GetPSFImageDataMaximum();
  } else if (m_DisplayedImage == CALCULATED_BSF_IMAGE) {
    value = m_DataModel->GetBSFImageDataMaximum();
  } else if (m_DisplayedImage == BSF_DIFFERENCE_IMAGE) {
    value = m_DataModel->GetBSFDifferenceImageDataMaximum();
  }
  return value;
}


void
PSFEstimator
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
PSFEstimator
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
PSFEstimator
::Exit() {
  // Ask if user really wants to quit.
  QMessageBox messageBox(this);
  messageBox.setText("Do you really want to exit?");
  messageBox.setInformativeText("If you exit now, all unsaved settings will be lost.");
  messageBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
  messageBox.setDefaultButton(QMessageBox::Cancel);
  int selected = messageBox.exec();

  if (selected == QMessageBox::Ok) {
    WriteProgramSettings();
    qApp->exit();
  }
}


void
PSFEstimator
::WriteProgramSettings() {
  QSettings settings;

  settings.beginGroup("MainWindow");
  settings.setValue("WindowSettings", saveState());
  settings.setValue("Geometry", saveGeometry());

  for (int i = 0; i < 4; i++) {
    QString colName;
    colName.sprintf("PSFSettingsTableViewColumnWidth%0d", i);
    settings.setValue(colName, gui->psfSettingsTableView->columnWidth(i));
  }

  settings.endGroup();
}


void
PSFEstimator
::ReadProgramSettings() {
  QSettings settings;

  settings.beginGroup("MainWindow");
  restoreState(settings.value("WindowSettings").toByteArray());
  restoreGeometry(settings.value("Geometry").toByteArray());

  for (int i = 0; i < 4; i++) {
    QString colName;
    colName.sprintf("PSFSettingsTableViewColumnWidth%0d", i);
    int width = settings.value(colName).toInt();
    if (width > 0)
      gui->psfSettingsTableView->setColumnWidth(i, width);
  }

  settings.endGroup();
}


void
PSFEstimator
::SaveFileChooserDirectory(const QString& path) {
  // Strip any file name from the end
  QString newPath = path;

  QFileInfo fileInfo(path);
  if (!fileInfo.isDir())
    newPath = fileInfo.dir().absolutePath();

  QSettings settings;
  settings.beginGroup("FileChooser");
  settings.setValue("path", newPath);
  settings.endGroup();
}


QString
PSFEstimator
::GetFileChooserDirectory() {
  QSettings settings;
  settings.beginGroup("FileChooser");
  QString path = settings.value("path", tr(".")).toString();
  settings.endGroup();

  return path;
}


void
PSFEstimator
::closeEvent(QCloseEvent* event) {
  Exit();

  // If we made it past the call above, the user clicked cancel.
  event->ignore();
}
