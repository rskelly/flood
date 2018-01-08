/*
 * flood_qgis.cpp
 *
 *  Created on: Dec 21, 2017
 *      Author: rob
 */

#include <QtCore/QString>
#include <QtCore/QObject>
#include <QtGui/QAction>

#include "flood_qgis.hpp"
#include "flood.hpp"

// Constants in flood_qgis.hpp needed in this include.
#include <qgsdockwidget.h>
#include <qgsproject.h>
#include <qgsrasterlayer.h>
#include <qgsmaplayerregistry.h>
#include <qgsmessagelog.h>

const QString FloodPlugin::s_name = QObject::tr("Flood Plugin");
const QString FloodPlugin::s_description = QObject::tr("Incremental flooding of raster terrains. Extracts raster and vector basins and locates spill points.");
const QString FloodPlugin::s_category = QObject::tr("Plugins");
const QString FloodPlugin::s_version = QObject::tr("Version 0.0.1");
const QString FloodPlugin::s_icon = "";
const QgisPlugin::PLUGINTYPE FloodPlugin::s_type = QgisPlugin::UI;

QGISEXTERN QgisPlugin* classFactory(QgisInterface* qgis_if) {
	return new FloodPlugin(qgis_if);
}

QGISEXTERN QString name() {
	return FloodPlugin::s_name;
}

QGISEXTERN QString description() {
	return FloodPlugin::s_description;
}

QGISEXTERN QString category() {
	return FloodPlugin::s_category;
}

QGISEXTERN QgisPlugin::PLUGINTYPE type() {
	return FloodPlugin::s_type;
}

QGISEXTERN QString version() {
	return FloodPlugin::s_version;
}

QGISEXTERN QString icon() {
	return FloodPlugin::s_icon;
}

QGISEXTERN void unload(QgisPlugin* plugin) {
	delete plugin;
}

FloodPlugin::~FloodPlugin() {
	removeInputRasterLayer();
}

void FloodPlugin::removeInputRasterLayer() {
	if(m_inputRasterLayer) {
		QgsMessageLog::instance()->logMessage("Removing input layer.");
		QgsMapLayerRegistry::instance()->removeMapLayer(m_inputRasterLayer); // Deletes m_inputRasterLayer.
		m_inputRasterLayer = nullptr;
	}
}

void FloodPlugin::addInputRasterLayer() {
	if(!m_inputRasterLayer) {
		QgsMessageLog::instance()->logMessage("Adding input layer.");
		m_inputRasterLayer = new QgsRasterLayer(m_inputRaster);
		QgsMapLayerRegistry::instance()->addMapLayer(m_inputRasterLayer, true);
	}
}

void FloodPlugin::unload() {
	removeInputRasterLayer();
}

void FloodPlugin::initGui() {

	QgsDockWidget* w = new QgsDockWidget(tr("Flood"), m_qgis->mainWindow());
	w->setObjectName("Flood");
	w->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	m_qgis->addDockWidget(Qt::LeftDockWidgetArea, w);

	m_dockWidget = new FloodDockWidget(this);
	m_dockWidget->setupUi(w);

	m_action = new QAction(QIcon(""), tr("Flood"), this);
	m_action->setWhatsThis(tr("Floods raster terrains to produce basins and spill points."));
	//connect(m_action, SIGNAL(triggered()), this, SLOT(StartOverlay()));
	m_qgis->addRasterToolBarIcon(m_action);
	m_qgis->addPluginToMenu(tr("&Flood"), m_action);
}

void FloodPlugin::StartOverlay() {

}

void FloodPlugin::setInputRaster(QString raster) {
	removeInputRasterLayer();
	m_inputRaster = raster;
	addInputRasterLayer();
}

void FloodPlugin::setInputSeeds(QString seeds) {
	m_inputSeeds = seeds;
}

void FloodPlugin::setOutputRaster(QString raster) {
	m_outputRaster = raster;
}

void FloodPlugin::setOutputVector(QString vector) {
	m_outputVector = vector;
}

void FloodPlugin::setOutputSpillPoints(QString spillPoints) {
	m_outputSpillPoints = spillPoints;
}

void FloodPlugin::setStartElevation(double elev) {
	m_startElev = elev;
}

void FloodPlugin::setEndElevation(double elev) {
	m_endElev = elev;
}

void FloodPlugin::setMinBasinArea(double area) {
	m_minBasinArea = area;
}
void FloodPlugin::setMaxSpillDistance(double dist) {
	m_maxSpillDist = dist;
}

void FloodPlugin::setStep(double step) {
	m_step = step;
}

void FloodPlugin::setThreads(int threads) {
	m_threads = threads;
}

void FloodPlugin::setOverwrite(bool overwrite) {
	m_overwrite = overwrite;
}

void FloodPlugin::start() {
	geo::flood::Flood config(m_inputRaster.toStdString(), m_inputSeeds.toStdString(),
			m_outputVector.toStdString(), m_outputRaster.toStdString(),
			m_outputSpillPoints.toStdString(),
			m_startElev, m_endElev, m_step, m_minBasinArea, m_maxSpillDist);
			config.flood(m_threads, false);
}

void FloodPlugin::cancel() {

}


void FloodDockWidget::setupUi(QgsDockWidget* widget) {
	Ui::FloodDockWidget::setupUi(widget);
	connect(fleInputRaster, SIGNAL(fileChanged(QString)), m_plugin, SLOT(setInputRaster(QString)));
	connect(fleInputSeeds, SIGNAL(fileChanged(QString)), m_plugin, SLOT(setInputSeeds(QString)));
	connect(fleOutputRaster, SIGNAL(fileChanged(QString)), m_plugin, SLOT(setOutputRaster(QString)));
	connect(fleOutputVector, SIGNAL(fileChanged(QString)), m_plugin, SLOT(setOutputVector(QString)));
	connect(fleOutputSpill, SIGNAL(fileChanged(QString)), m_plugin, SLOT(setOutputSpillPoints(QString)));
	connect(spnStartElev, SIGNAL(valueChanged(double)), m_plugin, SLOT(setStartElevation(double)));
	connect(spnEndElev, SIGNAL(valueChanged(double)), m_plugin, SLOT(setEndElevation(double)));
	connect(spnMinBasinArea, SIGNAL(valueChanged(int)), m_plugin, SLOT(setMinBasinArea(int)));
	connect(spnMaxSpillDist, SIGNAL(valueChanged(double)), m_plugin, SLOT(setMaxSpillDistance(double)));
	connect(spnStep, SIGNAL(valueChanged(double)), m_plugin, SLOT(setStep(double)));
	connect(spnThreads, SIGNAL(valueChanged(int)), m_plugin, SLOT(setThreads(int)));
	connect(chkOverwrite, SIGNAL(toggled(bool)), m_plugin, SLOT(setOverwrite(bool)));
	connect(btnStart, SIGNAL(clicked()), m_plugin, SLOT(start()));
	connect(btnCancel, SIGNAL(clicked()), m_plugin, SLOT(cancel()));
}


