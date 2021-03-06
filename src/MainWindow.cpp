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
		try {
			m_state->computeBranchAssignments(dwc);
		}
		catch(std::exception const & e) {
			std::cerr << "Error occured: " << e.what() << std::endl;
		}
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
	
	QGridLayout * cfgLayout = new QGridLayout(this);
	
	using DWC = DistanceWeightConfig;
	
	
	m_computeBtn = new QPushButton("Energize!", this);
	
	m_nodeSelection = new QComboBox();
	m_nodeWeight = new QComboBox();
	m_regionWeight = new QComboBox();
	m_branchWeight = new QComboBox();
	m_algoSelection = new QComboBox();
	
	m_nodeSelection->setPlaceholderText("Straßenknotenauswahl");
	m_nodeSelection->addItem("Alle", QVariant::fromValue(DWC::NS::All));
	m_nodeSelection->addItem("Nur Kreuzungen", QVariant::fromValue(DWC::NS::OnlyCrossRoads));
	
	m_nodeWeight->setPlaceholderText("PLZ-Distanz-Gewichtung");
	m_nodeWeight->addItem("Min", QVariant::fromValue(DWC::NWM::Min));
	m_nodeWeight->addItem("Max", QVariant::fromValue(DWC::NWM::Max));
	m_nodeWeight->addItem("Mean", QVariant::fromValue(DWC::NWM::Mean));
	m_nodeWeight->addItem("Median", QVariant::fromValue(DWC::NWM::Median));
	
	m_regionWeight->setPlaceholderText("PLZ-Gewichtung");
	m_regionWeight->addItem("Ohne", QVariant::fromValue(DWC::RWM::Equal));
	m_regionWeight->addItem("Bevölkerungsdichte", QVariant::fromValue(DWC::RWM::Inhabitants));
	
	m_branchWeight->setPlaceholderText("Filial-Gewichtung");
	m_branchWeight->addItem("Ohne", QVariant::fromValue(DWC::BWM::Equal));
	m_branchWeight->addItem("Mitarbeiteranzahl", QVariant::fromValue(DWC::BWM::Employees));
	
	m_algoSelection->setPlaceholderText("Zuordnungsalgorithmus");
	m_algoSelection->addItem("Nach Distanz", QVariant::fromValue(DWC::BAM::ByDistance));
	m_algoSelection->addItem("Evolutionär Gesamtkostenreduzierung", QVariant::fromValue(DWC::BAM::Evolutionary));
	m_algoSelection->addItem("MinCostFlow Gesamtkostenreduzierung", QVariant::fromValue(DWC::BAM::Optimal));

	cfgLayout->addWidget(m_nodeSelection, 0, 0);
	cfgLayout->addWidget(m_nodeWeight, 1, 0);
	cfgLayout->addWidget(m_regionWeight, 2, 0);
	cfgLayout->addWidget(m_branchWeight, 0, 1);
	cfgLayout->addWidget(m_algoSelection, 1, 1);
	cfgLayout->addWidget(m_computeBtn, 2, 1);
	
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
	QAction * openAction = new QAction("&Load", this);
	fileMenu->addAction(openAction);
	QAction * saveAction = new QAction("&Save", this);
	fileMenu->addAction(saveAction);
	QAction * exportAction = new QAction("&Export", this);
	fileMenu->addAction(exportAction);
	QAction * closeAction = new QAction("&Close", this);
	fileMenu->addAction(closeAction);
	
	auto editMenu = menuBar()->addMenu("&Edit");
	QAction * clearAction = new QAction("&Clear", this);
	editMenu->addAction(clearAction);
	
	connect(openAction, SIGNAL(triggered(bool)), this, SLOT(load()));
	connect(saveAction, SIGNAL(triggered(bool)), this, SLOT(save()));
	connect(exportAction, SIGNAL(triggered(bool)), this, SLOT(exportResults()));
	connect(closeAction, SIGNAL(triggered(bool)), this, SIGNAL(quit()));
	connect(clearAction, SIGNAL(triggered(bool)), this, SLOT(clear()));
	
	QWidget * centralWidget = new QWidget(this);
	centralWidget->setLayout(mainLayout);
	setCentralWidget(centralWidget);
}

MainWindow::~MainWindow() {
	emit quit();
}

void MainWindow::computeAssignments() {
	DistanceWeightConfig dwc;
#define E(__D, __S) dwc.__D = __S->currentData().value<decltype(dwc.__D)>();
	E(nodeSelection, m_nodeSelection)
	E(nodeWeightModel, m_nodeWeight)
	E(regionWeightModel, m_regionWeight)
	E(branchWeightModel, m_branchWeight)
	E(branchAssignmentModel, m_algoSelection)
#undef E
	m_bc->compute(dwc);
}

void MainWindow::computationFinished() {
	m_infoLabel->setText(QString("Total cost: %1").arg(m_state->totalCost));
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

void MainWindow::save() {
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

void MainWindow::exportResults() {
	std::cout << "Exporting assignments" << std::endl;
	QString fileName = QFileDialog::getSaveFileName(this, tr("Open destiation"), "/", tr("Text Files (*.txt)"));
	std::ofstream out;
	out.open(fileName.toStdString());
	if (out.is_open()) {
		m_state->exportBranchAssignments(out);
		out.close();
	}
	else {
		QMessageBox msgBox;
		msgBox.setText("Could not open file: " + fileName);
		msgBox.exec();
	}
}

void MainWindow::load() {
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

void MainWindow::clear() {
	m_state->clear();
}

void
MainWindow::showTextInfo(QString const & str) {
	m_infoLabel->setText(str);
}

}//end namespace
