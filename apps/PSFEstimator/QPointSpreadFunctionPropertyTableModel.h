#ifndef _Q_POINT_SPREAD_FUNCTION_PROPERTY_TABLE_MODEL_H_
#define _Q_POINT_SPREAD_FUNCTION_PROPERTY_TABLE_MODEL_H_

#include <QAbstractTableModel>

// Forward declarations
class DataModel;


class QPointSpreadFunctionPropertyTableModel : public QAbstractTableModel {
  Q_OBJECT

 public:
  QPointSpreadFunctionPropertyTableModel(QObject* parent = 0);
  ~QPointSpreadFunctionPropertyTableModel();

  void SetDataModel(DataModel* model);
  DataModel* GetDataModel();

  bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole);
\
  QVariant data(const QModelIndex& index, int role) const;

  Qt::ItemFlags flags(const QModelIndex& index) const;

  QVariant headerData(int section, Qt::Orientation orientation,
                      int role = Qt::DisplayRole) const;

  int rowCount(const QModelIndex& parent = QModelIndex()) const;
  int columnCount(const QModelIndex& parent = QModelIndex()) const;

  void Refresh();

  void InitializeSettingsCache();
  void SaveSettingsCache();

 protected:
  DataModel* m_DataModel;

  QList<double>   m_PropertyValues;
  QList<bool>     m_OptimizeValues;

};

#endif // _Q_POINT_SPREAD_FUNCTION_PROPERTY_TABLE_MODEL_H_

