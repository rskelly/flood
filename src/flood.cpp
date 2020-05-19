/*
 * flood.cpp
 *
 *  Created on: Aug 2, 2018
 *      Author: rob
 */

#include <cmath>

#include "ds/mqtree.hpp"
#include "grid.hpp"
#include "util.hpp"

#include "flood.hpp"

using namespace geo::flood::util;
using namespace geo::grid;
using namespace geo::ds;
using namespace geo::util;

namespace {
	/**
	 * Maintains the IDs for cells in the raster.
	 */
	size_t cellId = 0;
	std::mutex cellMtx;

	size_t nextCellId() {
		std::lock_guard<std::mutex> lk(cellMtx);
		return ++cellId;
	}


	bool readline(std::istream &st, std::vector<std::string> &row) {
		row.clear();
		std::string line;
		std::getline(st, line);
		std::stringstream linest(line);
		std::string cell;
		bool data = false;
		while (std::getline(linest, cell, ',')) {
			row.push_back(cell);
			data = true;
		}
		return data;
	}

} // anon

geo::flood::Cell::Cell() :
	m_col(0), m_row(0),
	m_cellId(nextCellId()),
	m_value(0),
	m_seedId(0) {
}

geo::flood::Cell::Cell(size_t seedId, int col, int row) :
	m_col(col), m_row(row),
	m_cellId(nextCellId()),
	m_value(0),
	m_seedId(seedId) {
}

geo::flood::Cell::Cell(size_t id, size_t seedId, int col, int row) :
	m_col(col), m_row(row),
	m_cellId(id),
	m_value(0),
	m_seedId(seedId) {
}

geo::flood::Cell::Cell(size_t id, size_t seedId, int col, int row, double value) :
	m_col(col), m_row(row),
	m_cellId(id),
	m_value(value),
	m_seedId(seedId) {
}

geo::flood::Cell::Cell(size_t seedId, int col, int row, double value) :
	m_col(col), m_row(row),
	m_cellId(nextCellId()),
	m_value(value),
	m_seedId(seedId) {
}

double geo::flood::Cell::operator[](int idx) const {
	switch(idx % 2) {
	case 0: return m_col;
	case 1: return m_row;
	}
	return 0;
}

double geo::flood::Cell::x() const {
	return m_col;
}

double geo::flood::Cell::y() const {
	return m_row;
}

int geo::flood::Cell::row() const {
	return m_row;
}

int geo::flood::Cell::col() const {
	return m_col;
}

size_t geo::flood::Cell::cellId() const {
	return m_cellId;
}

size_t geo::flood::Cell::seedId() const {
	return m_seedId;
}

double geo::flood::Cell::value() const {
	return m_value;
}

double geo::flood::Cell::distance(const Cell& other, double resx, double resy) const {
	double x0 = col() * resx;
	double y0 = row() * resy;
	double x1 = other.col() * resx;
	double y1 = other.row() * resy;
	return std::sqrt(geo::sq(x0 - x1) + geo::sq(y0 - y1));
}



Basin::Basin(unsigned int id, int minc, int minr, int maxc, int maxr, int area) :
	m_id(id),
	m_minc(minc), m_minr(minr),
	m_maxc(maxc), m_maxr(maxr),
	m_area(area),
	m_band(1) {
}

int Basin::computeEdges(Band<int>& grd, mqtree<Cell>& edgeCells) {
	// Offsets for searching. {row, col}
	static int offset[4][2] = {{-1, 0}, {0, -1}, {0, 1}, {1, 0}};
	// Get the grid dimensions.
	const GridProps& props = grd.props();
	int cols = props.cols();
	int rows = props.rows();
	// Iterate over pixels.
	int count = 0;
	std::ofstream out("seed_" + std::to_string(m_id) + ".csv", std::ios::out);
	out << std::setprecision(12);
	for (int r = m_minr; r <= m_maxr; ++r) {
		for (int c = m_minc; c <= m_maxc; ++c) {
			// Only look at pixels that match the seed ID.
			if ((unsigned int) grd.get(c, r) == m_id) {
				for (int i = 0; i < 4; ++i) {
					int co = c + offset[i][1];
					int ro = r + offset[i][0];
					// Check if this is an edge pixel.
					if(co < 0 || ro < 0 || co >= cols || ro >= rows || (unsigned int) grd.get(co, ro) != m_id) {
						edgeCells.add(Cell(m_id, c, r));
						out << props.toX(c) << "," << props.toY(r) << "\n";
						++count;
						break;
					}
				}
			}
		}
	}
	return count;
}

unsigned int Basin::seedId() const {
	return m_id;
}


LEFillOperator::LEFillOperator(Band<float>* src, int srcBand, Band<int>* dst, int dstBand, double elevation, unsigned int id) :
	m_src(src),
	m_dst(dst),
	m_elevation(elevation),
	m_id(id),
	m_nofill(0),
	m_srcBand(srcBand),
	m_dstBand(dstBand) {

	m_nodata = src->props().nodata();
}

bool LEFillOperator::shouldFill(int col, int row) const {
	if(!m_dst->props().hasCell(col, row))
		return false;
	int fill = m_dst->get(col, row);
	if(fill == m_nofill) {
		double value = m_src->get(col, row);
		return value != m_nodata && value <= m_elevation;
	}
	return false;
}

void LEFillOperator::fill(int col, int row) const {
	if(m_dst->props().hasCell(col, row))
		m_dst->set(col, row, (int) m_id);
}

void LEFillOperator::setFill(unsigned int fill) {
	m_id = fill;
}

void LEFillOperator::setNoFill(int nofill) {
	m_nofill = nofill;
}

const GridProps& LEFillOperator::srcProps() const {
	return m_src->props();
}

const GridProps& LEFillOperator::dstProps() const {
	return m_dst->props();
}

LEFillOperator::~LEFillOperator(){}



SpillPoint::SpillPoint(const Cell& c1, const Cell& c2, double elevation) :
		m_c1(c1),
		m_c2(c2),
		m_elevation(elevation) {
}

const geo::flood::Cell& SpillPoint::cell1() const {
	return m_c1;
}

const geo::flood::Cell& SpillPoint::cell2() const {
	return m_c2;
}

double SpillPoint::elevation() const {
	return m_elevation;
}

void SpillPoint::centroid(int* col, int* row) const {
	*col = (int) (m_c1.col() + (m_c2.col() - m_c1.col()) / 2.0);
	*row = (int) (m_c1.row() + (m_c2.row() - m_c2.row()) / 2.0);
}

SpillPoint::~SpillPoint() {
}


using namespace geo::flood;

Flood::Flood(const std::string& input, int band, const Output& output,
		const std::string& spill, const std::string& seeds,
		double start, double end, double step,
		double minBasinArea, double maxSpillDist,
		const std::vector<BreakLine>& breakLines) :
			m_start(start), m_end(end), m_step(step),
			m_minBasinArea(minBasinArea),
			m_maxSpillDist(maxSpillDist),
			m_t(1),
			m_band(band),
			m_output(&output),
			m_input(input),
			m_spill(spill),
			m_fseeds(seeds) {
	m_breakLines = breakLines;
}

Flood::Flood(const Flood& conf) :
	Flood(conf.m_input, conf.m_band, *conf.m_output, conf.m_spill, conf.m_fseeds,
		conf.m_start, conf.m_end, conf.m_step, conf.m_minBasinArea, conf.m_maxSpillDist,
		conf.m_breakLines) {
	// Copy the thread count.
	m_t = conf.m_t;
}

void Flood::validateInputs() {

	g_debug("Checking...");

	if(m_start > m_end)
		throw std::runtime_error("End elevation must be larger than start.");

	if(m_step <= 0)
		throw std::runtime_error("Step must be greater than zero.");

	if(m_minBasinArea < 0)
		throw std::runtime_error("Minimum basin area must be greater than or equal to zero.");

	if(m_maxSpillDist <= 0)
		throw std::runtime_error("Maximum spill distance must be larger than zero.");

	if (!m_output || m_output->valid())
		g_warn("The output configuration is not valid.");

	if (m_spill.empty())
		g_warn("No spill file; not producing spill points.");

	if (std::isnan(m_start))
		g_warn("Start value not given; using raster minimum.");

	if (std::isnan(m_end))
		g_warn("End value not given; using raster maximum.");

	unsigned int threads = std::thread::hardware_concurrency();
	if(m_t < 1)
		throw std::runtime_error("Threads must be greater than zero.");

	if(threads > 0 && m_t > threads)
		throw std::runtime_error("Threads must be less than or equal to the number of available cores.");

	if(m_input.empty())
		throw std::runtime_error("Input raster not given.");

	if(!isfile(m_input))
		throw std::runtime_error("Input raster not found.");

	if(!m_fseeds.empty() && !isfile(m_fseeds))
		throw std::runtime_error("Seeds file does not exist.");

	m_output->prepare();
}

Flood::~Flood() {

}

const std::string& Flood::input() const {
	return m_input;
}

const Output& Flood::output() const {
	return *m_output;
}

const std::string& Flood::spill() const {
	return m_spill;
}

const std::string& Flood::seedsFile() const {
	return m_fseeds;
}

double Flood::start() const {
	return m_start;
}

double Flood::end() const {
	return m_end;
}

double Flood::step() const {
	return m_step;
}

double Flood::minBasinArea() const {
	return m_minBasinArea;
}

double Flood::maxSpillDist() const {
	return m_maxSpillDist;
}

void Flood::loadSeeds(bool header) {
	if (m_fseeds.empty())
		g_argerr("No seed file given.");
	std::ifstream csv(m_fseeds);
	std::vector<std::string> row;
	if (header)
		readline(csv, row);
	const GridProps &props = m_dem.props();
	while (readline(csv, row)) {
		unsigned int id = (unsigned int) std::strtol(row[0].c_str(), nullptr, 10);
		double x = std::strtod(row[1].c_str(), nullptr);
		double y = std::strtod(row[2].c_str(), nullptr);
		g_debug("Found seed: " << id << ", " << x << ", " << y << ", " << props.toCol(x) << ", " << props.toRow(y));
		m_seeds.emplace_back(id, props.toCol(x), props.toRow(y));
	}
}

const std::vector<geo::flood::Cell>& Flood::seeds() const {
	return m_seeds;
}

Band<float>& Flood::dem() {
	return m_dem;
}

void Flood::init(std::mutex& mtx) {

	validateInputs();

	g_debug("Initing...");
	{
		std::lock_guard<std::mutex> lk(mtx);
		m_dem.init(m_input, m_band - 1, false, true);
		if(!m_breakLines.empty()) {
			for(const BreakLine& l : m_breakLines)
				m_dem.drawLine(l.x0, l.y0, l.x1, l.y1, l.value);
		}
	}

	g_debug("Computing stats...");
	const GridStats stats = m_dem.stats();
	if (std::isnan(m_start)) {
		m_start = stats.min;
	} else {
		m_start = geo::max(stats.min, m_start);
	}
	if (std::isnan(m_end)) {
		m_end = stats.max;
	} else {
		m_end = geo::min(stats.max, m_end);
	}
	g_debug("Stats: " << stats.min << ", " << stats.max);
	if (m_end < m_start)
		g_argerr("The ending elevation must be larger than the starting elevation: " << m_start << ", " << m_end);
}

int Flood::fillBasins(Band<int>& basinRaster, double elevation) {
	g_debug("Filling basins: " << elevation);

	basinRaster.fill(0);
	m_basinList.clear();

	// Set up the fill operator.
	LEFillOperator op(&m_dem, 1, &basinRaster, 1, elevation, 0);
	op.setNoFill(0);

	// Iterate over the seeds.
	for (const Cell& seed : seeds()) {

		if(Flood::cancel)
			break;

		if (!basinRaster.props().hasCell(seed.col(), seed.row())) {
			g_warn("Found a seed out of bounds: " << seed.col() << ", " << seed.row());
			continue;
		}

		// Fill the basin based on the elevations in the DEM.
		int minc, minr, maxc, maxr, area;
		op.setFill(seed.cellId());
		Band<float>::floodFill(seed.col(), seed.row(), op, true, &minc, &minr, &maxc, &maxr, &area);

		double barea = area * std::abs(m_dem.props().resX() * m_dem.props().resY());

		if (barea >= minBasinArea()) {
			// If it's large enough, save the basin.
			g_debug("Good Seed: " << seed.col() << ", " << seed.row() << "; elevation: " << elevation << "; basin area: " << barea);
			m_basinList.emplace_back(seed.cellId(), minc, minr, maxc, maxr, area);
		} else if(area > 0){
			//g_debug("Bad Seed: " << seed.col() << ", " << seed.row() << "; elevation: " << elevation << "; basin area: " << barea);
			// If the basin is too small, fill it with nodata. Do not collect more spill points.
			TargetFillOperator<int, int> fop(&basinRaster, 0, &basinRaster, 0, seed.cellId(), basinRaster.props().nodata());
			Band<int>::floodFill(seed.col(), seed.row(), fop, true, &minc, &minr, &maxc, &maxr, &area);
		}

	}

	return m_basinList.size();
}

bool Flood::findSpillPoints(Band<int>& basinRaster, double elevation) {
	g_debug("Finding spill points.");

	m_spillPoints.clear();

	const geo::util::Bounds bounds(0, 0, m_dem.props().cols(), m_dem.props().rows());

	std::unordered_map<int, std::unique_ptr<mqtree<Cell>>> trees;
	std::unordered_set<int> seen;

	// Compare each basin to each other basin.
	for (size_t i = 0; i < m_basinList.size(); ++i) {
		for (size_t j = i + 1; j < m_basinList.size() && !Flood::cancel; ++j) {

			//g_debug("Basins " << i << ", " << j);
			Basin& b0 = m_basinList[i];
			Basin& b1 = m_basinList[j];
			unsigned int id0 = b0.seedId();
			unsigned int id1 = b1.seedId();

			if(trees.find(id0) == trees.end()) {
				std::unique_ptr<mqtree<Cell>> t(new mqtree<Cell>(bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy()));
				b0.computeEdges(basinRaster, *t);
				trees.emplace(id0, std::move(t));
			}
			if(trees.find(id1) == trees.end()) {
				std::unique_ptr<mqtree<Cell>> t(new mqtree<Cell>(bounds.minx(), bounds.miny(), bounds.maxx(), bounds.maxy()));
				b1.computeEdges(basinRaster, *t);
				trees.emplace(id1, std::move(t));
			}

			int count0 = trees.at(id0)->size();
			int count1 = trees.at(id1)->size();

			// Compare the distances; save the ones that are near enough.
			mqtree<Cell>& cells0 = count0 < count1 ? *trees.at(id0) : *trees.at(id1);
			mqtree<Cell>& cells1 = count0 < count1 ? *trees.at(id1) : *trees.at(id0);
			std::vector<Cell> result;
			double spillDist = m_maxSpillDist / std::abs(basinRaster.props().resX());

			Cell c0;
			while(cells0.next(c0)) {
				if(Flood::cancel)
					break;
				seen.insert(c0.cellId());
				cells1.search(c0, spillDist, std::back_inserter(result));
				for(const Cell& c : result) {
					seen.insert(c.cellId());
					// Create a spill point; copies the cells.
					m_spillPoints.emplace_back(c0, c, elevation);
				}
				result.resize(0);
			}
		}
	}

	//g_debug("Minimum distance: " << minDist);
	return m_spillPoints.size();
}

/**
 * The heuristic calculates the distance and modifies it using the slope. A horizontal
 * has the normal 2d length. The distance is multiplied by the sine of the angle
 * of the slope w/r/t 1 map unit.
 */
class heuristic {
public:
	Band<float>* dem;
	heuristic(Band<float>* dem) : dem(dem) {}
	double operator()(int cc, int cr, int gc, int gr) {
		if(!dem->props().hasCell(cc, cr) || !dem->props().hasCell(gc, gr) ||
				dem->get(cc, cr) == dem->props().nodata() || dem->get(gc, gr) == dem->props().nodata())
			return geo::maxvalue<double>();
		double a = dem->get(cc, cr);
		double b = dem->get(gc, gr);
		double d = std::sqrt(geo::sq(cc - gc) + geo::sq(cr - gr));
		double m = std::sin(std::atan2(b - a, d));
		return d + m * d;
	}
};

void Flood::saveSpillPoints(unsigned int* id, std::ostream &out) {
	g_debug("Outputting spill points.");
	out << std::setprecision(12);
	const GridProps& props = m_dem.props();
	for (const SpillPoint& sp : m_spillPoints) {
		const Cell& c1 = sp.cell1();
		const Cell& c2 = sp.cell2();
		double x1 = props.toX(c1.col());
		double y1 = props.toY(c1.row());
		double x2 = props.toX(c2.col());
		double y2 = props.toY(c2.row());
		double x3 = (x1 + x2) / 2.0;
		double y3 = (y1 + y2) / 2.0;
		double dist = std::sqrt(geo::sq(x1 - x2) + geo::sq(y1 - y2));

		std::string pathStr = "";
		double maxElev = geo::minvalue<double>();
		{
			std::vector<std::pair<double, double>> path;
			heuristic h(&m_dem);
			m_dem.searchAStar(c1.col(), c1.row(), c2.col(), c2.row(), h, std::back_inserter(path), 1);
			if(!path.empty()) {
				std::stringstream ss;
				ss << std::setprecision(12) << "\"LINESTRING(";
				double v;
				for(size_t i = 0; i < path.size(); ++i) {
					double x = path[i].first;
					double y = path[i].second;
					if(i > 0) ss << ',';
					ss << x << ' ' << y;
					// Get the path's max elevation.
					if((v = m_dem.get(x, y)) > maxElev)
						maxElev = v;
				}
				ss << ")\"";
				pathStr = ss.str();
			}
		}
		out << ++(*id) << ", " << c1.seedId() << "," << x1 << "," << y1 << "," << c2.seedId() << ","
				<< x2 << "," << y2 << "," << x3 << "," << y3 << ","
				<< sp.elevation() << "," << dist << "," << maxElev << "," << pathStr << std::endl;
	}
}

void Flood::findMinima() {
	g_debug("Finding minima.");

	m_seeds.clear();
	int cols = m_dem.props().cols();
	int rows = m_dem.props().rows();
	double nodata = m_dem.props().nodata();
	for (int r = 0; r < rows && !Flood::cancel; ++r) {
		for (int c = 0; c < cols; ++c) {
			bool skip = false;
			double v;
			if ((v = m_dem.get(c, r)) == nodata)
				continue;
			for (int rr = geo::max(0, r - 1); !skip && rr < geo::min(r + 2, rows); ++rr) {
				for (int cc = geo::max(0, c - 1); !skip && cc < geo::min(c + 2, cols) && !Flood::cancel; ++cc) {
					double v0;
					if ((cc == c && rr == r) || (v0 = m_dem.get(cc, rr)) == nodata)
						continue;
					if (m_dem.get(cc, rr) < m_dem.get(c, r))
						skip = true;
				}
			}
			if (!skip)
				m_seeds.push_back(Cell(0, c, r, m_dem.get(c, r)));
		}
	}
}

void Flood::worker(Flood* config, std::mutex* mtx, std::ofstream* ofs, std::queue<double>* elevations) {

	Flood conf(*config);
	conf.init(*mtx);

	g_debug("Building seed list...");
	if (!conf.seedsFile().empty()) {
		// Load the seeds if given.
		conf.loadSeeds(true);
	} else {
		// Otherwise find and use minima.
		conf.findMinima();
	}

	bool run = true;

	while(run && !Flood::cancel) {
		double elevation;
		{
			std::lock_guard<std::mutex> lk(*mtx);
			if(elevations->empty()) {
				run = false;
			} else {
				elevation = elevations->front();
				elevations->pop();
			}
		}
		if(!run)
			break;

		g_debug("Filling to " << elevation << " by step " << conf.step());

		// Spill point ID. Needed?
		unsigned int id = 0;

		// The basin filename.
		std::string rfile = config->output().rasterFile(elevation);

		// Set up the output raster.
		GridProps props(conf.dem().props());
		props.setBands(1);
		props.setNoData(0);
		props.setDataType(DataType::UInt32);
		props.setWritable(true);
		props.setCompress(true);

		Band<int> basinRaster(rfile, props);

		// Generate basins.
		int basins = conf.fillBasins(basinRaster, elevation);

		basinRaster.flush();

		// Find and output spill points.
		if (basins > 1 && !conf.spill().empty() && conf.findSpillPoints(basinRaster, elevation) > 0) {
			std::lock_guard<std::mutex> lk(*mtx);
			conf.saveSpillPoints(&id, *ofs);
		}

		g_debug("Writing basin raster " << rfile);

		basinRaster.flush();

		// If desired, generate vectors.
		if (basins > 0 && conf.output().doVectors())
			conf.output().saveVector(basinRaster, elevation);
	}

}

void Flood::flood(int numThreads) {

	using namespace geo::flood::util;

	g_debug("Flooding...");

	std::list<std::thread> threads;
	std::ofstream ofs;
	std::queue<double> elevations;

	// Build the config object.
	Flood config(*this);
	config.validateInputs();

	if(!config.spill().empty()) {
		rem(spill());
		ofs.open(spill());
		ofs << "id,id1,x1,y1,id2,x2,y2,xmidpoint,ymidpoint,elevation,distance,maxElevation,path" << std::endl;
	}

	for(double e = start(); e <= end(); e += step())
		elevations.push(e);

	std::mutex mtx;
	for(int i = 0; i < numThreads; ++i)
		threads.emplace_back(Flood::worker, &config, &mtx, &ofs, &elevations);

	for(std::thread& th : threads)
		th.join();

}


