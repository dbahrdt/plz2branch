#include "BranchTableModel.h"
#include <assert.h>

namespace plz2branch {

BranchTableModel::BranchTableModel(QObject* parent, std::shared_ptr<State> state) :
QAbstractTableModel(parent),
m_state(std::move(state))
{}

BranchTableModel::~BranchTableModel() {}

int BranchTableModel::columnCount(const QModelIndex&) const {
	return CN_COL_COUNT;
}

int BranchTableModel::rowCount(const QModelIndex&) const {
	std::shared_lock<std::shared_mutex> lck(m_state->dataMtx);
	return m_state->branches.size();
}

QVariant BranchTableModel::headerData(int section, Qt::Orientation orientation, int role) const {
	if (role == Qt::DisplayRole) {
		if (orientation == Qt::Horizontal) {
			switch (section) {
			case (CN_ID):
				return QVariant("Id");
			case (CN_LAT):
				return QVariant("Latitude");
			case (CN_LON):
				return QVariant("Longitude");
			case (CN_NAME):
				return QVariant("Name");
			case CN_PLZ:
				return QVariant("PLZ");
			case CN_NUM_PLZ:
				return QVariant("#PLZ");
			case (CN_SHOW):
				return QVariant("Show");
			default:
				return QVariant();
			}
		}
		else {
			return QVariant();
		}
	}
	return QVariant();
}

QVariant BranchTableModel::data(const QModelIndex& index, int role) const {
	if (role == Qt::DisplayRole) {

		int col = index.column();
		int row = index.row();
		
		assert(row < (int) m_state->graph.nodeCount());
		
		std::shared_lock<std::shared_mutex> lck(m_state->dataMtx);
		Branch const & info = m_state->branches.at(row);
		switch(col) {
		case CN_ID:
			return QVariant(QString::number(row));
		case CN_LAT:
			return QVariant(info.coord.lat());
		case CN_LON:
			return QVariant(info.coord.lon());
		case CN_NAME:
			return QVariant(info.name);
		case CN_PLZ:
		{
			std::stringstream ss;
			for(auto const & [rId, rDist] : info.assignedRegions) {
				ss << m_state->regionInfo.at(rId.value).plz << ": " << rDist.value  << ", ";
			}
			return QVariant( QString::fromStdString(ss.str()) );
		}
		case CN_NUM_PLZ:
		{
			return QVariant((uint32_t) info.assignedRegions.size());
		}
		case CN_SHOW:
		{
			std::shared_lock<std::shared_mutex> lck(m_state->dataMtx);
			if (m_state->enabledBranches.count(row)) {
				return QVariant("Yes");
			}
			else {
				return QVariant("No");
			}
		}
		default:
			return QVariant();
		}
	}
	return QVariant();
}

bool
BranchTableModel::setData(const QModelIndex &index, const QVariant &value, int role) {
	if (role == Qt::EditRole) {
		int col = index.column();
		int row = index.row();
		if (col == CN_NAME) {
			std::unique_lock<std::shared_mutex> lck(m_state->dataMtx);
			if (row < m_state->branches.size()) {
				m_state->branches.at(row).name = value.toString();
				return true;
			}
		}
	}
	return QAbstractTableModel::setData(index, value, role);
}

Qt::ItemFlags
BranchTableModel::flags(const QModelIndex &index) const
{
	int col = index.column();

	if (col == CN_NAME) {
		return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
	}
	return QAbstractTableModel::flags(index);
}

void BranchTableModel::resetData() {
	emit beginResetModel();
	emit endResetModel();
}

void BranchTableModel::doubleClicked(const QModelIndex & index) {
	int column = index.column();
	int row = index.row();
		
	assert(row >= 0 && (unsigned int) row < m_state->graph.nodeCount());

	switch (column) {
	case CN_LAT:
	case CN_LON:
	{
		emit branchCoordinateClicked(BranchId{.value = uint32_t(row)});
		break;
	}
	case CN_SHOW:
	{
		emit toggleShow(BranchId{.value = uint32_t(row)});
		break;
	}
	default:
		break;
	}
}

void BranchTableModel::headerClicked(int logicalindex) {
	switch(logicalindex) {
	case CN_SHOW:
		emit clearShownClicked();
		break;
	default:
		break;
	}
}

}//end namespace
