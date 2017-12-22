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

#include <QtCore/QObject>
#include <QtGui/QAction>

#include "qgisplugin.h"
#include "qgisinterface.h"
#include "qgsmapcanvas.h"
#include "qgsmessagelog.h"

#include "ui_flood_qgis.h"

class FloodPlugin;

class FloodDockWidget : public QObject, public Ui::FloodDockWidget {
	Q_OBJECT
private:
	FloodPlugin* m_plugin;
public:
	FloodDockWidget(FloodPlugin* plugin) :
		Ui::FloodDockWidget(), m_plugin(plugin) {
	}

	void setupUi(QgsDockWidget* widget);

};

class FloodPlugin : public QObject, public QgisPlugin {
	Q_OBJECT
private:
	QgisInterface* m_qgis;
	QAction* m_action;
	FloodDockWidget* m_dockWidget;
	QgsRasterLayer* m_inputRasterLayer;

	QString m_inputRaster;
	QString m_inputSeeds;
	QString m_outputRaster;
	QString m_outputVector;
	QString m_outputSpillPoints;
	double m_startElev;
	double m_endElev;
	double m_minBasinArea;
	double m_maxSpillDist;
	int m_threads;
	bool m_overwrite;

public:
	static const QString s_name, s_description, s_category, s_version, s_icon;
	static const QgisPlugin::PLUGINTYPE s_type;

	FloodPlugin(QgisInterface* qgis) :
		QgisPlugin(s_name, s_description, s_category, s_version, s_type),
		m_qgis(qgis),
		m_action(nullptr),
		m_dockWidget(nullptr),
		m_inputRasterLayer(nullptr),
		m_startElev(0), m_endElev(0),
		m_minBasinArea(1), m_maxSpillDist(20),
		m_threads(1),
		m_overwrite(false) {

	}

	virtual ~FloodPlugin();
	virtual void addInputRasterLayer();
	virtual void removeInputRasterLayer();
	virtual void initGui();
	virtual void unload();

public slots:

	void StartOverlay();
	void DrawOverlay(QPainter* painter);

	void setInputRaster(QString);
	void setInputSeeds(QString);
	void setOutputRaster(QString);
	void setOutputVector(QString);
	void setOutputSpillPoints(QString);
	void setStartElevation(double);
	void setEndElevation(double);
	void setMinBasinArea(double);
	void setMaxSpillDistance(double);
	void setThreads(int);
	void setOverwrite(bool);
	void start();
	void cancel();
};




#endif /* INCLUDE_FLOOD_QGIS_HPP_ */
