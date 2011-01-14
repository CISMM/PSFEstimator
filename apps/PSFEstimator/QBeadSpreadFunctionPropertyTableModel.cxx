#include <QBeadSpreadFunctionPropertyTableModel.h>

#include <DataModel.h>


QBeadSpreadFunctionPropertyTableModel
::QBeadSpreadFunctionPropertyTableModel(QObject* parent)
  : QAbstractTableModel(parent) {

  m_DataModel = NULL;
}


QBeadSpreadFunctionPropertyTableModel
::~QBeadSpreadFunctionPropertyTableModel() {
}


void
QBeadSpreadFunctionPropertyTableModel
::SetDataModel(DataModel* model) {
  m_DataModel = model;

  InitializeSettingsCache();

  Refresh();
}


DataModel*
QBeadSpreadFunctionPropertyTableModel
::GetDataModel() {
  return m_DataModel;
}


bool
QBeadSpreadFunctionPropertyTableModel
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
QBeadSpreadFunctionPropertyTableModel
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
        return QVariant(m_DataModel->GetParameterName(row).c_str());
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
        return QVariant(m_DataModel->GetParameterUnit(row).c_str());
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
QBeadSpreadFunctionPropertyTableModel
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
QBeadSpreadFunctionPropertyTableModel
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
QBeadSpreadFunctionPropertyTableModel
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
QBeadSpreadFunctionPropertyTableModel
::columnCount(const QModelIndex& parent) const {
  return 4;
}


void
QBeadSpreadFunctionPropertyTableModel
::Refresh() {
  reset();
}


void
QBeadSpreadFunctionPropertyTableModel
::InitializeSettingsCache() {
  m_PropertyValues.clear();
  m_OptimizeValues.clear();

  for (int i = 0; i < m_DataModel->GetNumberOfProperties(); i++) {
    m_PropertyValues.append(m_DataModel->GetParameterValue(i));
    m_OptimizeValues.append(m_DataModel->GetParameterEnabled(i));
  }
}



void
QBeadSpreadFunctionPropertyTableModel
::SaveSettingsCache() {
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
    m_DataModel->SetParameterValue(i, m_PropertyValues[i]);
    m_DataModel->SetParameterEnabled(i, m_OptimizeValues[i]);
  }
}
