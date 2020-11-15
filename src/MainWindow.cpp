#include "MainWindow.h"
#include <QMessageBox>
#include <QTableView>
#include <QComboBox>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QFileDialog>
#include <QPushButton>
#include <QMenuBar>
#include <iostream>
#include <future>

#include "MarbleMap.h"
#include "BranchTableModel.h"

namespace plz2branch {
namespace detail {

BackgroundComputer::BackgroundComputer(QObject * parent, std::shared_ptr<State> state) :
QObject(parent),
m_state(std::move(state))
{}

BackgroundComputer::~BackgroundComputer() {}

void BackgroundComputer::compute(DistanceWeightConfig const & dwc) {
	(void) std::async(std::launch::async, [this](DistanceWeightConfig const & dwc) {
		std::cout << "Computing branch assignments in thread " << std::this_thread::get_id() << std::endl;
		m_state->computeBranchAssignments(dwc);
		finished();
	}, dwc);
}

}//end namespace

MainWindow::MainWindow(std::shared_ptr<State> state):
QMainWindow(),
m_state(state),
m_bc(new detail::BackgroundComputer(this, state))
{
	m_map = new MarbleMap(this, m_state);
	m_map->setMapThemeId("earth/openstreetmap/openstreetmap.dgml");
	
	m_branchTableModel = new BranchTableModel(this, m_state);
	
	m_branchesTableView = new QTableView(this);
	m_branchesTableView->setModel(m_branchTableModel);


	m_branchesTableView->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
	
	QHBoxLayout * cfgLayout = new QHBoxLayout(this);
	
	using DWC = DistanceWeightConfig;
	m_weightModel = new QComboBox(this);
	m_weightModel->addItem("All Nodes, Min Distance", QVariant::fromValue(DWC(DWC::All, DWC::Min)));
	m_weightModel->addItem("All Nodes, Max Distance", QVariant::fromValue(DWC(DWC::All, DWC::Max)));
	m_weightModel->addItem("All Nodes, Mean Distance", QVariant::fromValue(DWC(DWC::All, DWC::Mean)));
	m_weightModel->addItem("All Nodes, Median Distance", QVariant::fromValue(DWC(DWC::All, DWC::Median)));
	m_weightModel->addItem("Only Crossroads, Min Distance", QVariant::fromValue(DWC(DWC::OnlyCrossRoads, DWC::Min)));
	m_weightModel->addItem("Only Crossroads, Max Distance", QVariant::fromValue(DWC(DWC::OnlyCrossRoads, DWC::Max)));
	m_weightModel->addItem("Only Crossroads, Mean Distance", QVariant::fromValue(DWC(DWC::OnlyCrossRoads, DWC::Mean)));
	m_weightModel->addItem("Only Crossroads, Median Distance", QVariant::fromValue(DWC(DWC::OnlyCrossRoads, DWC::Median)));
	
	m_computeBtn = new QPushButton("Energize!", this);
	
	cfgLayout->addWidget(m_weightModel);
	cfgLayout->addWidget(m_computeBtn);
	
	QWidget * cfgWidget = new QWidget(this);
	cfgWidget->setLayout(cfgLayout);
	
	QHBoxLayout * cfgStatsLayout = new QHBoxLayout(this);
	cfgStatsLayout->addWidget(cfgWidget);
	
	QWidget * cfgStatsWidget = new QWidget(this);
	cfgStatsWidget->setLayout(cfgStatsLayout);
	
	m_infoLabel = new QLabel();
	
	QVBoxLayout * leftLayout = new QVBoxLayout(this);
	leftLayout->addWidget(m_branchesTableView, 5);
	leftLayout->addWidget(cfgStatsWidget, 1);
	leftLayout->addWidget(m_infoLabel, 1);
	
	QWidget * leftLayoutWidget = new QWidget(this);
	leftLayoutWidget->setLayout(leftLayout);
	
	QHBoxLayout * mainLayout = new QHBoxLayout(this);
	mainLayout->addWidget(leftLayoutWidget, 2);
	mainLayout->addWidget(m_map, 3);
	
	connect(m_branchesTableView, SIGNAL(doubleClicked(QModelIndex)), m_branchTableModel, SLOT(doubleClicked(QModelIndex)));
	connect(m_branchesTableView->horizontalHeader(), SIGNAL(sectionClicked(int)), m_branchTableModel, SLOT(headerClicked(int)));

	connect(m_branchTableModel, SIGNAL(branchCoordinateClicked(BranchId)), m_map, SLOT(zoomToBranch(BranchId)));	
	connect(m_branchTableModel, SIGNAL(toggleShow(BranchId)), this, SLOT(toggleBranch(BranchId)));
	connect(m_branchTableModel, SIGNAL(clearShownClicked()), this, SLOT(clearShownBranches()));

	//map
	connect(m_map, SIGNAL(branchCreated(double, double)), m_state.get(), SLOT(createBranch(double, double)));
	
	//connect background router
	connect(m_computeBtn, SIGNAL(clicked()), this, SLOT(computeAssignments()));
	connect(m_bc, SIGNAL(finished()), this, SLOT(computationFinished()));
	
	//connect data changed
	connect(m_state.get(), SIGNAL(dataChanged()), this, SIGNAL(dataChanged()));
	connect(this, SIGNAL(dataChanged()), m_branchTableModel, SLOT(resetData()));
	connect(this, SIGNAL(dataChanged()), m_map, SLOT(dataChanged()));
	
	connect(m_state.get(), SIGNAL(textInfo(QString const&)), m_infoLabel, SLOT(setText(QString)));
	
	auto fileMenu = menuBar()->addMenu("&File");
	QAction * openAction = new QAction("&Open", this);
	fileMenu->addAction(openAction);
	QAction * saveAction = new QAction("&Save", this);
	fileMenu->addAction(saveAction);
	
	connect(openAction, SIGNAL(triggered(bool)), this, SLOT(loadFile()));
	connect(saveAction, SIGNAL(triggered(bool)), this, SLOT(saveAssignments()));
	
	QWidget * centralWidget = new QWidget(this);
	centralWidget->setLayout(mainLayout);
	setCentralWidget(centralWidget);
}

MainWindow::~MainWindow() {}

void MainWindow::computeAssignments() {
	m_bc->compute(m_weightModel->currentData().value<DistanceWeightConfig>());
}

void MainWindow::computationFinished() {
	emit dataChanged();
}

void MainWindow::clearShownBranches() {
	{
		std::shared_lock<std::shared_mutex> lck(m_state->dataMtx);
		m_state->enabledBranches.clear();
	}
	emit dataChanged();
}

void MainWindow::toggleBranch(BranchId brId) {
	{
		std::shared_lock<std::shared_mutex> lck(m_state->dataMtx);
		if (m_state->enabledBranches.count(brId.value)) {
			m_state->enabledBranches.erase(brId.value);
		}
		else {
			m_state->enabledBranches.insert(brId.value);
		}
	}
	emit dataChanged();
}

void MainWindow::saveAssignments() {
	std::cout << "Saving assignments" << std::endl;
	QString fileName = QFileDialog::getSaveFileName(this, tr("Open destiation"), "/", tr("Text Files (*.txt)"));
	std::ofstream out;
	out.open(fileName.toStdString());
	if (out.is_open()) {
		m_state->writeBranchAssignments(out);
		out.close();
	}
	else {
		QMessageBox msgBox;
		msgBox.setText("Could not open file: " + fileName);
		msgBox.exec();
	}
}

void MainWindow::loadFile() {
	std::cout << "Load branches" << std::endl;
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open destiation"), "/", tr("Text Files (*.txt)"));
	std::ifstream in;
	in.open(fileName.toStdString());
	if (in.is_open()) {
		m_state->createBranches(in);
		in.close();
	}
	else {
		QMessageBox msgBox;
		msgBox.setText("Could not open file: " + fileName);
		msgBox.exec();
	}
}

void
MainWindow::showTextInfo(QString const & str) {
	m_infoLabel->setText(str);
}

}//end namespace
