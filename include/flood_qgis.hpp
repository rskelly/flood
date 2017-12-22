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

class FloodPlugin : public QObject, public QgisPlugin {

	Q_OBJECT

public:
	static const QString s_name, s_description, s_category, s_version, s_icon;
	static const QgisPlugin::PLUGINTYPE s_type;

	FloodPlugin(QgisInterface* qgis_if) :
		QgisPlugin(s_name, s_description, s_category, s_version, s_type),
		m_qgis_if(qgis_if) {
	}

	virtual void initGui();
	virtual void unload();

public slots:

	void StartOverlay();
	void DrawOverlay(QPainter* painter);

private:
	QgisInterface* m_qgis_if;
	QAction* m_action;

};


#endif /* INCLUDE_FLOOD_QGIS_HPP_ */
