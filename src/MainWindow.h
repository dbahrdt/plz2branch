#pragma once

#include <QMainWindow>
#include <memory>
#include "State.h"

class QTableView;
class QComboBox;
class QLabel;
class QPushButton;

namespace plz2branch {
namespace detail {

class BackgroundComputer: public QObject {
	Q_OBJECT
public:
	BackgroundComputer(QObject * parent, std::shared_ptr<State> state);
	~BackgroundComputer() override;
public slots:
	void compute(DistanceWeightConfig const & dwc);
signals:
	void finished();
private:
	std::shared_ptr<State> m_state;
};

}//end namespace detail

class MarbleMap;
class BranchTableModel;

class MainWindow: public QMainWindow {
	Q_OBJECT
public:
	MainWindow(std::shared_ptr<State> state);
	~MainWindow() override;
public slots:
	void computeAssignments();
	void toggleBranch(BranchId brId);
	void clearShownBranches();
	void saveAssignments();
	void loadFile();
	void showTextInfo(QString const &);
signals:
	void dataChanged();
private slots:
	void computationFinished();
private://data stuff
	std::shared_ptr<State> m_state;
	detail::BackgroundComputer * m_bc;
private://gui stuff
	MarbleMap * m_map;
	BranchTableModel * m_branchTableModel;
	QTableView * m_branchesTableView;
	
	QPushButton * m_computeBtn;
	
	QComboBox * m_weightModel;
	QComboBox * m_nodeSelection;
	
	QLabel * m_infoLabel;
};

};
