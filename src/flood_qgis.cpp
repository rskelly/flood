/*
 * flood_qgis.cpp
 *
 *  Created on: Dec 21, 2017
 *      Author: rob
 */


#include "flood_qgis.hpp"

// Constants in flood_qgis.hpp needed in this include.
#include "qgsdockwidget.h"


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

void FloodPlugin::unload() {

}

void FloodPlugin::initGui() {

	QgsDockWidget* w = new QgsDockWidget(tr("Flood"), m_qgis->mainWindow());
	w->setObjectName("Flood");
	w->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	m_qgis->addDockWidget(Qt::LeftDockWidgetArea, w);

	m_dockWidget = new Ui::FloodDockWidget();
	m_dockWidget->setupUi(w);

	m_action = new QAction(QIcon(""), tr("Flood"), this);
	m_action->setWhatsThis(tr("Floods raster terrains to produce basins and spill points."));
	//connect(m_action, SIGNAL(triggered()), this, SLOT(StartOverlay()));
	m_qgis->addRasterToolBarIcon(m_action);
	m_qgis->addPluginToMenu(tr("&Flood"), m_action);
}

void FloodPlugin::StartOverlay() {

}




