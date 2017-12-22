/*
 * flood_qgis.hpp
 *
 *  Created on: Dec 21, 2017
 *      Author: rob
 */

#ifndef INCLUDE_FLOOD_QGIS_HPP_
#define INCLUDE_FLOOD_QGIS_HPP_

// These are ususally (un)defined in the build files. No idea what they do.
#define CORE_EXPORT
#define GUI_EXPORT

#include <QObject>
#include <QAction>

#include "qgisplugin.h"
#include "qgisinterface.h"
#include "qgsmapcanvas.h"
#include "qgsmessagelog.h"

#include "ui_flood_qgis.h"

class FloodPlugin : public QObject, public QgisPlugin {

	Q_OBJECT

private:
	QgisInterface* m_qgis;
	QAction* m_action;
	Ui::FloodDockWidget* m_dockWidget;

public:
	static const QString s_name, s_description, s_category, s_version, s_icon;
	static const QgisPlugin::PLUGINTYPE s_type;

	FloodPlugin(QgisInterface* qgis) :
		QgisPlugin(s_name, s_description, s_category, s_version, s_type),
		m_qgis(qgis),
		m_action(nullptr),
		m_dockWidget(nullptr) {
	}

	virtual void initGui();
	virtual void unload();

public slots:

	void StartOverlay();
	void DrawOverlay(QPainter* painter);

};


#endif /* INCLUDE_FLOOD_QGIS_HPP_ */
