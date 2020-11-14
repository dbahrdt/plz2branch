#pragma once

#include <QAbstractTableModel>
#include "State.h"

namespace plz2branch {

class BranchTableModel: public QAbstractTableModel {
	Q_OBJECT
private:
	typedef enum { CN_ID=0, CN_LAT=1, CN_LON=2, CN_NAME=3, CN_PLZ=4, CN_SHOW=5, CN_COL_COUNT=CN_SHOW+1} ColNames;
public:
	BranchTableModel(QObject * parent, std::shared_ptr<State> state);
	virtual ~BranchTableModel();
	virtual int rowCount(const QModelIndex&) const;
	virtual int columnCount(const QModelIndex&) const;
	virtual QVariant data(const QModelIndex & index, int role) const;
	virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const;
	virtual Qt::ItemFlags flags(const QModelIndex& /*index*/) const {
		return (Qt::ItemIsEnabled);
	}
public Q_SLOTS:
	void doubleClicked(const QModelIndex&);
	void headerClicked(int);
	void resetData();
Q_SIGNALS:
	void branchCoordinateClicked(BranchId brId);
	void toggleShow(BranchId brId);
	void clearShownClicked();
private:
	std::shared_ptr<State> m_state;
};

}
