#include "MainWindow.h"
#include <iostream>
#include <QApplication>
#include <QFile>
#include <QMetaType>


void help() {
	std::cout << "plz2branch [options] file.osm.pbf\n";
	std::cout << std::endl;
}

int main(int argc, char ** argv) {
	plz2branch::registerTypesWithQt();
	
	QApplication app(argc, argv);
	QStringList cmdline_args = QCoreApplication::arguments();
	
	if (cmdline_args.size() < 2) {
		std::cerr << "Usage: plz2branch <path to osm pbf file> [<path to plz population>]" << std::endl;
		return -1;
	}
	if (!QFile::exists(cmdline_args.at(1))) {
		std::cerr << "Selected file " << cmdline_args.at(1).toStdString() << " does not exist or is not accessible" << std::endl;
		return -1;
	}
	
	plz2branch::Config cfg;
	cfg.osmFile = cmdline_args.at(1).toStdString();
	if (cmdline_args.size() > 2) {
		if (!QFile::exists(cmdline_args.at(2))) {
			std::cerr << "Selected file " << cmdline_args.at(2).toStdString() << " does not exist or is not accessible" << std::endl;
			return -1;
		}
		else {
			cfg.plzFile = cmdline_args.at(2).toStdString();
		}
	}
	
	
	auto state = std::make_shared<plz2branch::State>(cfg);

	state->importFiles();

	auto mainWindow = new plz2branch::MainWindow(state);
	mainWindow->setAttribute(Qt::WA_DeleteOnClose);
	
	QObject::connect(mainWindow, &plz2branch::MainWindow::quit, &app, &QApplication::quit, Qt::QueuedConnection);
	
	mainWindow->show();
	return app.exec();
}
