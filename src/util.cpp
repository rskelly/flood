/*
 * util.cpp
 *
 *  Created on: Jun 12, 2020
 *      Author: rob
 */

#include "flood.hpp"

using namespace geo::flood::util;

heuristic::heuristic(Band<float>* dem) :
	dem(dem) {
}

float heuristic::operator()(int cc, int cr, int gc, int gr) {
	const GridProps& props = dem->props();

	// If the cell is out of range, or the start or end are null,
	// return maximum cost.
	if(!props.hasCell(cc, cr) || !props.hasCell(gc, gr) ||
			dem->get(cc, cr) == props.nodata() || dem->get(gc, gr) == props.nodata())
		return geo::maxvalue<float>();

	// Calculate and return the heuristic.
	float a = dem->get(cc, cr);
	float b = dem->get(gc, gr);
	float d = std::sqrt(geo::sq(cc - gc) + geo::sq(cr - gr));
	float m = std::sin(std::atan2(b - a, d));
	return d + m * d;
}


bool BasinOutput::valid() const {
	return !rdir.empty();
}

void BasinOutput::prepare() {
	if(!rdir.empty() && !isdir(rdir))
		makedir(rdir);
	if(!isdir(rdir))
		throw std::runtime_error("Raster output directory does not exist and could not be created.");
}

std::string BasinOutput::rasterFile(int elevation) const {
	std::stringstream ss;
	ss << rdir << "/" << elevation << ".tif";
	return ss.str();
}


bool BasinDBOutput::valid() const {
	return BasinOutput::valid() && !(rdir.empty() && conn.empty()
			&& layer.empty() && idField.empty() && elevationField.empty());
}

void BasinDBOutput::prepare() {
	if(!rdir.empty() && !isdir(rdir))
		makedir(rdir);
	if(!isdir(rdir))
		throw std::runtime_error("Raster output directory does not exist and could not be created.");
}

void BasinDBOutput::saveVectors(Band<int>& rast, float elev) {
	std::vector<PolygonValue> fields = {{elevationField, elev}};
	rast.polygonizeToTable(conn, layer, idField, fields);
}


bool SpillDBOutput::valid() const {
	return !(projection.empty() || conn.empty() || layer.empty() ||
			idField.empty() || elevationField.empty() || maxElevationField.empty());
}

void SpillDBOutput::prepare() {
	GDALAllRegister();
	sr.importFromWkt(projection.c_str());
	ds = (GDALDataset*) GDALOpenEx(conn.c_str(), GDAL_OF_VECTOR|GDAL_OF_UPDATE, nullptr, nullptr, nullptr);
	lyr = ds->GetLayerByName(layer.c_str());
}

void SpillDBOutput::saveSpillPoints(Band<float>& dem, const std::vector<SpillPoint>& spillPoints) {
	g_debug("Outputting spill points.");
	if(OGRERR_NONE != lyr->StartTransaction())
		g_runerr("Failed to start transaction.");
	for (const SpillPoint& sp : spillPoints) {
		const Cell& c1 = sp.cell1();
		const Cell& c2 = sp.cell2();
		std::string pathStr = "";
		float maxElev = geo::minvalue<float>();
		{
			std::vector<std::pair<float, float>> path;
			heuristic h(&dem);
			dem.searchAStar(c1.col(), c1.row(), c2.col(), c2.row(), h, std::back_inserter(path), 1);
			if(!path.empty()) {
				std::stringstream ss;
				ss << std::setprecision(12) << "LINESTRING(";
				float v;
				for(size_t i = 0; i < path.size(); ++i) {
					float x = path[i].first;
					float y = path[i].second;
					if(i > 0) ss << ',';
					ss << x << ' ' << y;
					// Get the path's max elevation.
					if((v = dem.get(x, y)) > maxElev)
						maxElev = v;
				}
				ss << ")ne";
				pathStr = ss.str();
			}
		}
		OGRFeature feat(lyr->GetLayerDefn());
		feat.SetField(idField.c_str(), (long long) ++id);
		feat.SetField(elevationField.c_str(), sp.elevation());
		feat.SetField(maxElevationField.c_str(), maxElev);
		feat.SetField(bid1Field.c_str(), (long long) c1.seedId());
		feat.SetField(bid2Field.c_str(), (long long) c2.seedId());
		OGRGeometry* geom;
		if(OGRERR_NONE != gf.createFromWkt(pathStr.c_str(), &sr, &geom))
			g_warn("Failed to make linestring.");
		if(OGRERR_NONE != feat.SetGeometry(geom))
			g_warn("Failed to set geometry.");
		if(OGRERR_NONE != lyr->CreateFeature(&feat))
			g_warn("Failed to add geometry.");
		gf.destroyGeometry(geom);
	}
	if(OGRERR_NONE != lyr->CommitTransaction())
		g_runerr("Failed to commit transaction.");
}

SpillDBOutput::~SpillDBOutput() {
	GDALClose(ds);
}


bool BasinFileOutput::valid() const {
	return BasinOutput::valid() && !vdir.empty();
}

void BasinFileOutput::prepare() {
	BasinOutput::prepare();
	if(!vdir.empty() && !isdir(vdir))
		makedir(vdir);
	if(!vdir.empty() && !isdir(vdir))
		throw std::runtime_error("Vector output directory does not exist and could not be created.");
}

void BasinFileOutput::saveVectors(Band<int>& rast, float elevation) {
	std::stringstream ss;
	ss << vdir << "/" << (int) (elevation * 10000) << ".sqlite";
	std::string vfile = ss.str();
	rast.polygonizeToFile(vfile, "basins", "bid", "SQLite");
}

bool DummyBasinOutput::valid() const {
	return true;
}

void DummyBasinOutput::saveVectors(Band<int>&, float) {
	// no-op
}


bool DummySpillOutput::valid() const {
	return true;
}

void DummySpillOutput::prepare() {
	// no-op
}

void DummySpillOutput::saveSpillPoints(Band<float>&, const std::vector<SpillPoint>&) {
	// no-op
}

bool SpillFileOutput::valid() const {
	return !spillFile.empty();
}

void SpillFileOutput::prepare() {
	out.open(spillFile, std::ios::out);
	out << "id,id1,x1,y1,id2,x2,y2,xmidpoint,ymidpoint,elevation,distance,maxElevation,path" << std::endl;
	out << std::setprecision(12);
}

void SpillFileOutput::saveSpillPoints(Band<float>& dem, const std::vector<SpillPoint>& spillPoints) {
	g_debug("Outputting spill points.");
	const GridProps& props = dem.props();
	for (const SpillPoint& sp : spillPoints) {
		const Cell& c1 = sp.cell1();
		const Cell& c2 = sp.cell2();
		float x1 = props.toX(c1.col());
		float y1 = props.toY(c1.row());
		float x2 = props.toX(c2.col());
		float y2 = props.toY(c2.row());
		float x3 = (x1 + x2) / 2.0;
		float y3 = (y1 + y2) / 2.0;
		float dist = std::sqrt(geo::sq(x1 - x2) + geo::sq(y1 - y2));

		std::string pathStr = "";
		float maxElev = geo::minvalue<float>();
		{
			std::vector<std::pair<float, float>> path;
			heuristic h(&dem);
			dem.searchAStar(c1.col(), c1.row(), c2.col(), c2.row(), h, std::back_inserter(path), 1);
			if(!path.empty()) {
				std::stringstream ss;
				ss << std::setprecision(12) << "\"LINESTRING(";
				float v;
				for(size_t i = 0; i < path.size(); ++i) {
					float x = path[i].first;
					float y = path[i].second;
					if(i > 0) ss << ',';
					ss << x << ' ' << y;
					// Get the path's max elevation.
					if((v = dem.get(x, y)) > maxElev)
						maxElev = v;
				}
				ss << ")\"";
				pathStr = ss.str();
			}
		}
		out << ++id << ", " << c1.seedId() << "," << x1 << "," << y1 << "," << c2.seedId() << ","
				<< x2 << "," << y2 << "," << x3 << "," << y3 << ","
				<< sp.elevation() << "," << dist << "," << maxElev << "," << pathStr << std::endl;
	}
}


BreakLine::BreakLine(float x0, float y0, float x1, float y1, float value) :
	x0(x0), y0(y0), x1(x1), y1(y1), value(value) {
}

