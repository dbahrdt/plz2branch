#pragma once

#include <QAbstractTableModel>
#include "State.h"

namespace plz2branch {

class BranchTableModel: public QAbstractTableModel {
	Q_OBJECT
private:
	typedef enum { CN_ID=0, CN_SHOW, CN_LAT, CN_LON, CN_NAME, CN_EMPLOYEES, CN_MAX_PLZ, CN_PLZ, CN_NUM_PLZ, CN_COST, CN_COL_COUNT} ColNames;
public:
	BranchTableModel(QObject * parent, std::shared_ptr<State> state);
	~BranchTableModel() override;
	int rowCount(const QModelIndex&) const override;
	int columnCount(const QModelIndex&) const override;
	QVariant data(const QModelIndex & index, int role) const override;
	QVariant headerData(int section, Qt::Orientation orientation, int role) const override;
	Qt::ItemFlags flags(const QModelIndex& /*index*/) const override;
	bool setData(const QModelIndex &index, const QVariant &value, int role = Qt::EditRole) override;
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
