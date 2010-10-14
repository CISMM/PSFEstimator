#include <QPointSpreadFunctionPropertyTableModel.h>

#include <DataModel.h>


QPointSpreadFunctionPropertyTableModel
::QPointSpreadFunctionPropertyTableModel(QObject* parent)
  : QAbstractTableModel(parent) {

  m_DataModel = NULL;

  m_PropertyNameList.append(QVariant("X pixel size"));
  m_UnitsList.append(QVariant("nanometers"));

  m_PropertyNameList.append(QVariant("Y pixel size"));
  m_UnitsList.append(QVariant("nanometers"));

  m_PropertyNameList.append(QVariant("Z slice spacing"));
  m_UnitsList.append(QVariant("nanometers"));

  m_PropertyNameList.append(QVariant("Bead radius"));
  m_UnitsList.append(QVariant("nanometers"));

  m_PropertyNameList.append(QVariant("Bead center X"));
  m_UnitsList.append(QVariant("nanometers"));

  m_PropertyNameList.append(QVariant("Bead center Y"));
  m_UnitsList.append(QVariant("nanometers"));

  m_PropertyNameList.append(QVariant("Bead center Z"));
  m_UnitsList.append(QVariant("nanometers"));

  m_PropertyNameList.append(QVariant("Shear X"));
  m_UnitsList.append(QVariant("nanometers in X vs nanometers in Z"));

  m_PropertyNameList.append(QVariant("Shear Y"));
  m_UnitsList.append(QVariant("nanometers in Y vs nanometers in Z"));

  m_PropertyNameList.append(QVariant("Emission Wavelength"));
  m_UnitsList.append(QVariant("nanometers"));

  m_PropertyNameList.append(QVariant("Numerical Aperture"));
  m_UnitsList.append(QVariant("unitless"));

  m_PropertyNameList.append(QVariant("Magnification"));
  m_UnitsList.append(QVariant("unitless"));

  m_PropertyNameList.append(QVariant("Design Cover Slip Refractive Index"));
  m_UnitsList.append(QVariant("unitless"));

  m_PropertyNameList.append(QVariant("Actual Cover Slip Refractive Index"));
  m_UnitsList.append(QVariant("unitless"));

  m_PropertyNameList.append(QVariant("Design Cover Slip Thickness"));
  m_UnitsList.append(QVariant("micrometers"));

  m_PropertyNameList.append(QVariant("Actual Cover Slip Thickness"));
  m_UnitsList.append(QVariant("micrometers"));

  m_PropertyNameList.append(QVariant("Design Immersion Oil Refractive Index"));
  m_UnitsList.append(QVariant("unitless"));

  m_PropertyNameList.append(QVariant("Actual Immersion Oil Refractive Index"));
  m_UnitsList.append(QVariant("unitless"));

  m_PropertyNameList.append(QVariant("Design Immersion Oil Thickness"));
  m_UnitsList.append(QVariant("micrometers"));

  m_PropertyNameList.append(QVariant("Design Specimen Layer Refractive Index"));
  m_UnitsList.append(QVariant("unitless"));

  m_PropertyNameList.append(QVariant("Actual Specimen Layer Refractive Index"));
  m_UnitsList.append(QVariant("unitless"));

  m_PropertyNameList.append(QVariant("Actual Point Source Depth in Specimen Layer"));
  m_UnitsList.append(QVariant("micrometers"));

  m_PropertyNameList.append(QVariant("Design Distance from Back Focal Plane to Detector"));
  m_UnitsList.append(QVariant("millimeters"));

  m_PropertyNameList.append(QVariant("Actual Distance from Back Focal Plane to Detector"));
  m_UnitsList.append(QVariant("millimeters"));

  m_PropertyNameList.append(QVariant("Intensity Shift"));
  m_UnitsList.append(QVariant("-"));

  m_PropertyNameList.append(QVariant("Intensity Scale"));
  m_UnitsList.append(QVariant("-"));
}


QPointSpreadFunctionPropertyTableModel
::~QPointSpreadFunctionPropertyTableModel() {
}


void
QPointSpreadFunctionPropertyTableModel
::SetDataModel(DataModel* model) {
  m_DataModel = model;

  m_PropertyValues.clear();
  m_OptimizeValues.clear();
  for (int i = 0; i < m_DataModel->GetNumberOfProperties(); i++) {
    m_PropertyValues.append(0.0);
    m_OptimizeValues.append(false);
  }

  Refresh();
}


DataModel*
QPointSpreadFunctionPropertyTableModel
::GetDataModel() {
  return m_DataModel;
}


bool
QPointSpreadFunctionPropertyTableModel
::setData(const QModelIndex& index, const QVariant& value, int role) {
  int numPSFProperties = m_DataModel->GetNumberOfProperties();

  int row = index.row();
  int col = index.column();
  if (col == 1) {

    if (row < numPSFProperties) {
      m_PropertyValues[row] = value.toDouble();
    } else {
      if (m_DataModel->GetUseCustomZCoordinates())
        m_DataModel->SetZCoordinate(row - numPSFProperties, value.toDouble());
    }

  } else if (col == 3) {
    m_OptimizeValues[row] = value.toBool();
  }

  // Notify Qt items that the data has changed.
  emit dataChanged(index, index);

  return true;
}


QVariant
QPointSpreadFunctionPropertyTableModel
::data(const QModelIndex& index, int role) const {
  if (m_DataModel == NULL) {
    return QVariant();
  }

  int row = index.row();
  int col = index.column();
  if (row < 0 || row >= this->rowCount() || col < 0 || col >= this->columnCount())
    return QVariant();

  int numPSFProperties = m_DataModel->GetNumberOfProperties();

  if (role == Qt::DisplayRole || role == Qt::EditRole) {
    if (col == 0) { /** Property name **/

      if (row < numPSFProperties) {
        return m_PropertyNameList[row];
      } else {
        if (m_DataModel->GetUseCustomZCoordinates())
          return QString().sprintf("Z slice %d position", row - numPSFProperties);
      }

    } else if (col == 1) { /** Property value **/

      if (row < numPSFProperties) {
        return QVariant(m_PropertyValues[row]);
      } else {
        if (m_DataModel->GetUseCustomZCoordinates())
          return QVariant(m_DataModel->GetZCoordinate(row - numPSFProperties));
      }

    } else if (col == 2) { /** Property unit **/

      if (row < numPSFProperties) {
        return m_UnitsList[row];
      } else {
        if (m_DataModel->GetUseCustomZCoordinates())
          return QString("nanometers");
      }

    } else {
      return QVariant();
    }
  } else if (role == Qt::CheckStateRole) {

    if (col == 3 && row < numPSFProperties) { /** Optimize column **/
      return m_OptimizeValues[row] ? Qt::Checked : Qt::Unchecked;
    } else {
      return QVariant();
    }

  }

  return QVariant();
}


Qt::ItemFlags
QPointSpreadFunctionPropertyTableModel
::flags(const QModelIndex& index) const {
  Qt::ItemFlags flag = Qt::ItemIsSelectable | Qt::ItemIsEnabled;

  if (index.column() == 1) {
    flag = flag | Qt::ItemIsEditable;
    return flag;
  } else if (index.column() == 3) {
    return Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsUserCheckable; 
  }

  return Qt::ItemIsEnabled;
}


QVariant
QPointSpreadFunctionPropertyTableModel
::headerData(int section, Qt::Orientation orientation, int role) const {
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole) {
    switch (section) {
    case 0: return QVariant("Parameter"); break;
    case 1: return QVariant("Value"); break;
    case 2: return QVariant("Units"); break;
    case 3: return QVariant("Optimize?"); break;
    default: return QVariant(section); break;
    }
  } else if (orientation == Qt::Vertical && role == Qt::DisplayRole) {
    return QVariant(section);
  } else {
    return QVariant();
  }
}


int
QPointSpreadFunctionPropertyTableModel
::rowCount(const QModelIndex& parent) const {
  if (m_DataModel == NULL)
    return 0;

  int count, dims[3];
  count = m_DataModel->GetNumberOfProperties();
  if (m_DataModel->GetUseCustomZCoordinates()) {
    m_DataModel->GetMeasuredImageDimensions(dims);
    count += dims[2];
  }

  return count;
}


int
QPointSpreadFunctionPropertyTableModel
::columnCount(const QModelIndex& parent) const {
  return 4;
}


void
QPointSpreadFunctionPropertyTableModel
::Refresh() {
  reset();
}


void
QPointSpreadFunctionPropertyTableModel
::InitializeSettingsCache() {
  double triplet[3];

  m_PropertyValues.clear();
  m_OptimizeValues.clear();

  m_DataModel->GetMeasuredImageVoxelSpacing(triplet);
  for (int i = 0; i < 3; i++)
    m_PropertyValues.append(triplet[i]);

  m_PropertyValues.append(m_DataModel->GetBeadRadius());

  m_DataModel->GetBSFPointCenter(triplet);
  for (int i = 0; i < 3; i++)
    m_PropertyValues.append(triplet[i]);

  m_PropertyValues.append(m_DataModel->GetShearX());
  m_PropertyValues.append(m_DataModel->GetShearY());
  m_PropertyValues.append(m_DataModel->GetGLEmissionWavelength());
  m_PropertyValues.append(m_DataModel->GetGLNumericalAperture());
  m_PropertyValues.append(m_DataModel->GetGLMagnification());
  m_PropertyValues.append(m_DataModel->GetGLDesignCoverSlipRefractiveIndex());
  m_PropertyValues.append(m_DataModel->GetGLActualCoverSlipRefractiveIndex());
  m_PropertyValues.append(m_DataModel->GetGLDesignCoverSlipThickness());
  m_PropertyValues.append(m_DataModel->GetGLActualCoverSlipThickness());
  m_PropertyValues.append(m_DataModel->GetGLDesignImmersionOilRefractiveIndex());
  m_PropertyValues.append(m_DataModel->GetGLActualImmersionOilRefractiveIndex());
  m_PropertyValues.append(m_DataModel->GetGLDesignImmersionOilThickness());
  m_PropertyValues.append(m_DataModel->GetGLDesignSpecimenLayerRefractiveIndex());
  m_PropertyValues.append(m_DataModel->GetGLActualSpecimenLayerRefractiveIndex());
  m_PropertyValues.append(m_DataModel->GetGLActualPointSourceDepthInSpecimenLayer());
  m_PropertyValues.append(m_DataModel->GetGLDesignDistanceFromBackFocalPlaneToDetector());
  m_PropertyValues.append(m_DataModel->GetGLActualDistanceFromBackFocalPlaneToDetector());
  m_PropertyValues.append(m_DataModel->GetGLIntensityShift());
  m_PropertyValues.append(m_DataModel->GetGLIntensityScale());

  for (int i = 0; i < m_DataModel->GetNumberOfProperties(); i++) {
    m_OptimizeValues.append(m_DataModel->GetGLParameterEnabled(i));
  }
}



void
QPointSpreadFunctionPropertyTableModel
::SaveSettingsCache() {
  double triplet[3];
  int item = 0;

  for (int i = 0; i < 3; i++)
    triplet[i] = m_PropertyValues[item++];
  m_DataModel->SetMeasuredImageVoxelSpacing(triplet);
  m_DataModel->SetPSFImageVoxelSpacing(triplet);
  m_DataModel->SetBSFImageVoxelSpacing(triplet);

  m_DataModel->SetBeadRadius(m_PropertyValues[item++]);

  for (int i = 0; i < 3; i++)
    triplet[i] = m_PropertyValues[item++];
  m_DataModel->SetPSFPointCenter(triplet);
  m_DataModel->SetBSFPointCenter(triplet);

  m_DataModel->SetShearX(m_PropertyValues[item++]);
  m_DataModel->SetShearY(m_PropertyValues[item++]);
  m_DataModel->SetGLEmissionWavelength(m_PropertyValues[item++]);
  m_DataModel->SetGLNumericalAperture(m_PropertyValues[item++]);
  m_DataModel->SetGLMagnification(m_PropertyValues[item++]);
  m_DataModel->SetGLDesignCoverSlipRefractiveIndex(m_PropertyValues[item++]);
  m_DataModel->SetGLActualCoverSlipRefractiveIndex(m_PropertyValues[item++]);
  m_DataModel->SetGLDesignCoverSlipThickness(m_PropertyValues[item++]);
  m_DataModel->SetGLActualCoverSlipThickness(m_PropertyValues[item++]);
  m_DataModel->SetGLDesignImmersionOilRefractiveIndex(m_PropertyValues[item++]);
  m_DataModel->SetGLActualImmersionOilRefractiveIndex(m_PropertyValues[item++]);
  m_DataModel->SetGLDesignImmersionOilThickness(m_PropertyValues[item++]);
  m_DataModel->SetGLDesignSpecimenLayerRefractiveIndex(m_PropertyValues[item++]);
  m_DataModel->SetGLActualSpecimenLayerRefractiveIndex(m_PropertyValues[item++]);
  m_DataModel->SetGLActualPointSourceDepthInSpecimenLayer(m_PropertyValues[item++]);
  m_DataModel->SetGLDesignDistanceFromBackFocalPlaneToDetector(m_PropertyValues[item++]);
  m_DataModel->SetGLActualDistanceFromBackFocalPlaneToDetector(m_PropertyValues[item++]);
  m_DataModel->SetGLIntensityShift(m_PropertyValues[item++]);
  m_DataModel->SetGLIntensityScale(m_PropertyValues[item++]);

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


  for (int i = 0; i < m_DataModel->GetNumberOfProperties(); i++) {
    m_DataModel->SetGLParameterEnabled(i, m_OptimizeValues[i]);
  }
}
