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

using namespace geo::flood;
using namespace geo::flood::util;
using namespace geo::grid;
using namespace geo::ds;
using namespace geo::util;

namespace {

	int cellId = 0;			///<! Maintains the IDs for cells in the raster.
	std::mutex cellIdMtx;	///<! Protects the cell ID variable.

	/**
	 * \brief Return the next cell ID.
	 *
	 * \return The next cell ID.
	 */
	int nextCellId(int id) {
		std::lock_guard<std::mutex> lk(cellIdMtx);
		if(id == 0) {
			return ++cellId;
		} else {
			if(id >= cellId) cellId = id;
			return id;
		}
	}

	/**
	 * \brief Read a line of comma-delimited text from the stream into the row vector.
	 *
	 * \param st An input stream.
	 * \param[out] row A vector to contain the values of the current row.
	 * \return True if a line was read.
	 */
	bool readline(std::istream& st, std::vector<std::string> &row) {
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

	/**
	 * \brief Check if the column exists and is required.
	 *
	 * If the column is missing and required, an exception is raised.
	 * If the column is missing and not required returns false.
	 * Else, returns true.
	 *
	 * \param col The name of the column.
	 * \param colMap The mapping of columns to indices.
	 * \param colReq The mapping of columns to required.
	 * \return True if present. False if missing and not required.
	 */
	bool checkCol(const std::string& col,
			const std::vector<std::string>& row,
			const std::unordered_map<std::string, int>& colMap,
			const std::unordered_map<std::string, bool>& colReq) {
		if(colMap.find(col) == colMap.end() || colMap.at(col) >= (int) row.size()) {
			// Column is missing.
			if(colReq.find(col) != colReq.end() && colReq.at(col)) {
				// Column is required.
				g_runerr("Missing seed column " << col << " is required.");
			} else {
				return false;
			}
		}
		return true;
	}

	int getIntCol(const std::string& col, const std::vector<std::string>& row,
			const std::unordered_map<std::string, int>& colMap,
			const std::unordered_map<std::string, bool>& colReq,
			int def) {
		if(checkCol(col, row, colMap, colReq))
			return (int) std::strtol(row[colMap.at(col)].c_str(), nullptr, 10);
		return def;
	}

	int getFloatCol(const std::string& col, const std::vector<std::string>& row,
			const std::unordered_map<std::string, int>& colMap,
			const std::unordered_map<std::string, bool>& colReq,
			float def) {
		if(checkCol(col, row, colMap, colReq))
			return (float) std::strtod(row[colMap.at(col)].c_str(), nullptr);
		return def;
	}

	int getBoolCol(const std::string& col, const std::vector<std::string>& row,
			const std::unordered_map<std::string, int>& colMap,
			const std::unordered_map<std::string, bool>& colReq,
			bool def) {
		if(checkCol(col, row, colMap, colReq)) {
			std::string v = lowercase(row[colMap.at(col)]);
			return v == "true" || v == "t" || std::strtol(v.c_str(), nullptr, 10) != 0;
		}
		return def;
	}

	void worker(Flood* config, std::queue<int>* elevations) {

		while(!config->cancel) {
			int elev;
			{
				std::lock_guard<std::mutex> lk(config->qmtx());
				if(elevations->empty()) {
					break;
				} else {
					elev = elevations->front();
					elevations->pop();
				}
			}
			if(config->cancel)
				break;

			int p = std::pow(10, config->precision());
			float elevation = (float) elev / p;

			g_debug("Filling to " << elevation << " by step " << config->step());

			// The basin filename.
			std::string rfile = config->basinOutput().rasterFile(elev);

			// Do not overwite if not required.
			if(!config->overwrite() && geo::util::isfile(rfile))
				continue;

			// Set up the output raster.
			GridProps props(config->dem().props());
			props.setBands(1);
			props.setNoData(0);
			props.setDataType(DataType::UInt32);
			props.setWritable(true);
			props.setCompress(true);

			Band<int> basinRaster(rfile, props);

			// Generate basins.
			int basins = config->fillBasins(basinRaster, elevation);

			basinRaster.flush();

			// Find and output spill points.
			if (basins > 1 && config->findSpillPoints(basinRaster, elevation) > 0) {
				config->spillOutput().saveSpillPoints(config->dem(), config->spillPoints());
			}

			g_debug("Writing basin raster " << rfile);

			basinRaster.flush();

			// If desired, generate vectors.
			if (basins > 0)
				config->basinOutput().saveVectors(basinRaster, elevation);
		}

	}

} // anon

geo::flood::Cell::Cell() :
		Cell(0, 0, 0, 0, 0) {
}

geo::flood::Cell::Cell(int id, int col, int row, float value, int priority) :
	m_col(col), m_row(row),
	m_cellId(nextCellId(id)),
	m_value(value),
	m_priority(priority > 0 ? priority : geo::maxvalue<int>()) {
}

double geo::flood::Cell::operator[](int idx) const {
	// Return the dimension according to its index.
	switch(idx % 2) {
	case 0: return m_col;
	default: return m_row;
	}
}

bool geo::flood::Cell::operator<(const Cell& other) const {
	return m_priority < other.m_priority;
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

int geo::flood::Cell::cellId() const {
	return m_cellId;
}

float geo::flood::Cell::value() const {
	return m_value;
}

int geo::flood::Cell::priority() const {
	return m_priority;
}

void geo::flood::Cell::setPriority(int priority) {
	m_priority = priority;
}

float geo::flood::Cell::distance(const Cell& other, float resx, float resy) const {
	float x0 = col() * resx;
	float y0 = row() * resy;
	float x1 = other.col() * resx;
	float y1 = other.row() * resy;
	return std::sqrt(geo::sq(x0 - x1) + geo::sq(y0 - y1));
}



Basin::Basin(int id, int minc, int minr, int maxc, int maxr, int area) :
	m_id(id),
	m_minc(minc), m_minr(minr),
	m_maxc(maxc), m_maxr(maxr),
	m_area(area) {
}

int Basin::computeEdges(Band<int>& grd, mqtree<Cell>& edgeCells) {
	// Offsets for searching. {row, col}
	static int offset[4][2] = {{-1, 0}, {0, -1}, {0, 1}, {1, 0}};

	// Get the grid dimensions.
	const GridProps& props = grd.props();
	int cols = props.cols();
	int rows = props.rows();

	// Count the number of edge px found.
	int count = 0;

	// TODO: Maybe enable for a debug mode.
	// std::ofstream out("seed_" + std::to_string(m_id) + ".csv", std::ios::out);
	// out << std::setprecision(12);

	for (int r = m_minr; r <= m_maxr; ++r) {
		for (int c = m_minc; c <= m_maxc; ++c) {
			// Only look at pixels that match the seed ID.
			if (grd.get(c, r) == m_id) {
				for (int i = 0; i < 4; ++i) {
					int co = c + offset[i][1];
					int ro = r + offset[i][0];
					// Check if this is an edge pixel.
					if(co < 0 || ro < 0 || co >= cols || ro >= rows || grd.get(co, ro) != m_id) {
						edgeCells.add(Cell(m_id, c, r, 0, 0));
						++count;
						break;
					}
				}
			}
		}
	}
	return count;
}

int Basin::seedId() const {
	return m_id;
}


LEFillOperator::LEFillOperator(Band<float>* src, Band<int>* dst, float elevation, int id) :
	m_src(src),
	m_dst(dst),
	m_elevation(elevation),
	m_id(id),
	m_target(0) {

	m_nodata = src->props().nodata();
}

bool LEFillOperator::shouldFill(int col, int row) const {
	if(!m_dst->props().hasCell(col, row))
		return false;
	int fill = m_dst->get(col, row);
	if(fill == m_target) {
		float value = m_src->get(col, row);
		return value != m_nodata && value <= m_elevation;
	}
	return false;
}

void LEFillOperator::fill(int col, int row) const {
	if(m_dst->props().hasCell(col, row))
		m_dst->set(col, row, (int) m_id);
}

void LEFillOperator::setFill(int fill) {
	m_id = fill;
}

void LEFillOperator::setTarget(int target) {
	m_target = target;
}

const GridProps& LEFillOperator::srcProps() const {
	return m_src->props();
}

const GridProps& LEFillOperator::dstProps() const {
	return m_dst->props();
}

LEFillOperator::~LEFillOperator(){}



SpillPoint::SpillPoint(const Cell& c1, const Cell& c2, float elevation) :
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

float SpillPoint::elevation() const {
	return m_elevation;
}

void SpillPoint::centroid(int& col, int& row) const {
	col = (int) (m_c1.col() + (m_c2.col() - m_c1.col()) / 2.0);
	row = (int) (m_c1.row() + (m_c2.row() - m_c2.row()) / 2.0);
}

SpillPoint::~SpillPoint() {
}


Flood::Flood(const std::string& input, int band, bool overwrite,
		BasinOutput* basinOutput, SpillOutput* spillOutput,
		const std::string& seeds,
		int start, int end, int step, int precision,
		float minBasinArea, float maxSpillDist,
		const std::vector<BreakLine>& breakLines) :
			m_start(start), m_end(end), m_step(step),
			m_precision(precision),
			m_minBasinArea(minBasinArea),
			m_maxSpillDist(maxSpillDist),
			m_t(1),
			m_band(band),
			m_overwrite(overwrite),
			m_basinOutput(basinOutput),
			m_spillOutput(spillOutput),
			m_input(input),
			m_fseeds(seeds),
			cancel(false) {
	m_breakLines = breakLines;
}

void Flood::validateInputs() {

	g_debug("Checking...");

	if(m_start > m_end)
		g_runerr("End elevation must be larger than start.");

	if(m_step <= 0)
		g_runerr("Step must be greater than zero.");

	if(m_precision <= 0)
		g_runerr("Precision must be >0.");

	if(m_minBasinArea < 0)
		g_runerr("Minimum basin area must be greater than or equal to zero.");

	if(m_maxSpillDist <= 0)
		g_runerr("Maximum spill distance must be larger than zero.");

	if (!m_basinOutput || m_basinOutput->valid())
		g_warn("The basin output configuration is not valid.");

	if (!m_spillOutput || m_spillOutput->valid())
		g_warn("The spill output configuration is not valid.");

	if (std::isnan(m_start))
		g_warn("Start value not given; using raster minimum.");

	if (std::isnan(m_end))
		g_warn("End value not given; using raster maximum.");

	int threads = std::thread::hardware_concurrency();
	if(m_t < 1)
		g_runerr("Threads must be greater than zero.");

	if(threads > 0 && m_t > threads)
		g_runerr("Threads must be less than or equal to the number of available cores.");

	if(m_input.empty())
		g_runerr("Input raster not given.");

	if(!isfile(m_input))
		g_runerr("Input raster not found.");

	if(!m_fseeds.empty() && !isfile(m_fseeds))
		g_runerr("Seeds file does not exist.");

	m_basinOutput->prepare();
	m_spillOutput->prepare();
}

Flood::~Flood() {

}

std::mutex& Flood::qmtx() {
	return m_qmtx;
}

bool Flood::overwrite() const {
	return m_overwrite;
}

const std::string& Flood::input() const {
	return m_input;
}

BasinOutput& Flood::basinOutput() {
	return *m_basinOutput;
}

SpillOutput& Flood::spillOutput() {
	return *m_spillOutput;
}

const std::string& Flood::seedsFile() const {
	return m_fseeds;
}

int Flood::start() const {
	return m_start;
}

int Flood::precision() const {
	return m_precision;
}

int Flood::end() const {
	return m_end;
}

int Flood::step() const {
	return m_step;
}

float Flood::minBasinArea() const {
	return m_minBasinArea;
}

float Flood::maxSpillDist() const {
	return m_maxSpillDist;
}

const std::vector<SpillPoint>& Flood::spillPoints() const {
	return m_spillPoints;
}

void Flood::loadSeeds(bool header) {
	if (m_fseeds.empty())
		g_argerr("No seed file given.");

	std::ifstream csv(m_fseeds);
	std::vector<std::string> row;

	// Holds the row indices of each column.
	std::unordered_map<std::string, int> colMap;
	// Indicates which columns are required.
	std::unordered_map<std::string, bool> colReq = {{"gid", true}, {"x", true}, {"y", true}, {"priority", false}, {"keep", false}};

	if (header) {
		// Assign the indices to columns.
		readline(csv, row);
		for(size_t i = 0; i < row.size(); ++i) {
			std::string c = lowercase(row[i]);
			if(colReq.find(c) != colReq.end())
				colMap.emplace(c, i);
		}
	} else {
		// Columns in the seed file are expected in this order.
		std::vector<std::string> colOrder = {"gid", "x", "y", "priority", "keep"};
		for(size_t i = 0; i < colOrder.size(); ++i)
			colMap.emplace(colOrder[i], i);
	}

	const GridProps &props = m_dem.props();
	int ridx = 0;
	while (readline(csv, row)) {
		bool keep = getBoolCol("keep", row, colMap, colReq, true);
		if(!keep)
			continue;
		int id = getIntCol("gid", row, colMap, colReq, -1);
		if(id <= 0)
			g_runerr("Invalid or missing seed ID (" << id << ") on row " << ridx);
		float x = getFloatCol("x", row, colMap, colReq, std::nan(""));
		if(std::isnan(x))
			g_runerr("Invalid or missing seed X (" << x << ") on row " << ridx);
		float y = getFloatCol("y", row, colMap, colReq, std::nan(""));
		if(std::isnan(y))
			g_runerr("Invalid or missing seed Y (" << y << ") on row " << ridx);
		int priority = getIntCol("priority", row, colMap, colReq, 0);
		g_debug("Found seed: " << id << ", " << x << ", " << y << ", " << props.toCol(x) << ", " << props.toRow(y) << ", " << priority);
		m_seeds.emplace_back(id, props.toCol(x), props.toRow(y), 0, priority);
		++ridx;
	}

	// Sort the points without priority by sorting on their elevations and
	// assiging a lower priority to the ones with a higher elevation.
	g_debug("Sorting points...");
	int maxpriority = 0;
	float nd = (float) m_dem.props().nodata();
	std::vector<std::pair<Cell*, float>> ecells;
	for(Cell& c : m_seeds) {
		if(c.priority() == geo::maxvalue<int>()) {
			// Add the non-priority Cells and elevations to the list.
			float e = m_dem.get(c.col(), c.row());
			ecells.emplace_back(&c, e == nd ? geo::minvalue<float>() : e);
		}else if(c.priority() > maxpriority) {
			// Find the max set priority.
			maxpriority = c.priority();
		}
	}

	// Reverse-sort on elevation.
	std::sort(ecells.begin(), ecells.end(), [](const std::pair<Cell*, float>& a, const std::pair<Cell*, float>& b) {
		return a.second > b.second;
	});

	// Assign the priorities.
	for(size_t i = 0; i < ecells.size(); ++i)
		ecells[i].first->setPriority(i + maxpriority + 1);

	// Re-sort on priority.
	std::sort(m_seeds.begin(), m_seeds.end());
}

const std::vector<geo::flood::Cell>& Flood::seeds() const {
	return m_seeds;
}

Band<float>& Flood::dem() {
	return m_dem;
}

void Flood::init() {

	validateInputs();

	g_debug("Initing...");
	Band<float> dem(m_input, m_band - 1, false, true);
	m_dem.init("/tmp/smooth.tif", dem.props()); //m_input, m_band - 1, false, true);
	dem.smooth(m_dem);
	m_dem.flush();

	if(!m_breakLines.empty()) {
		for(const BreakLine& l : m_breakLines)
			m_dem.drawLine(l.x0, l.y0, l.x1, l.y1, l.value);
	}

	m_spillOutput->projection = m_dem.props().projection();

	g_debug("Computing stats...");
	const GridStats stats = m_dem.stats();
	int p = std::pow(10, m_precision);
	m_start = geo::max((int) (stats.min * p), m_start);
	m_end = geo::min((int) (stats.max * p), m_end);

	g_debug("Stats: " << stats.min << ", " << stats.max);
	if (m_end < m_start)
		g_argerr("The ending elevation must be larger than the starting elevation: " << m_start << ", " << m_end);
}

int Flood::fillBasins(Band<int>& basinRaster, float elevation) {
	g_debug("Filling basins: " << elevation);

	basinRaster.fill(0);
	m_basinList.clear();

	// Set up the fill operator.
	LEFillOperator op(&m_dem, &basinRaster, elevation, 0);
	op.setTarget(0);

	// Iterate over the seeds.
	for (const Cell& seed : seeds()) {

		if(cancel)
			break;

		if (!basinRaster.props().hasCell(seed.col(), seed.row())) {
			g_warn("Found a seed out of bounds: " << seed.col() << ", " << seed.row());
			continue;
		}

		// Fill the basin based on the elevations in the DEM.
		int minc, minr, maxc, maxr, area;
		op.setFill(seed.cellId());
		Band<float>::floodFill(seed.col(), seed.row(), op, true, &minc, &minr, &maxc, &maxr, &area);

		float barea = area * std::abs(m_dem.props().resX() * m_dem.props().resY());

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

bool Flood::findSpillPoints(Band<int>& basinRaster, float elevation) {
	g_debug("Finding spill points.");

	m_spillPoints.clear();

	const geo::util::Bounds bounds(0, 0, m_dem.props().cols(), m_dem.props().rows());

	std::unordered_map<int, std::unique_ptr<mqtree<Cell>>> trees;

	// Compare each basin to each other basin.
	for (size_t i = 0; i < m_basinList.size(); ++i) {
		for (size_t j = i + 1; j < m_basinList.size() && !cancel; ++j) {

			//g_debug("Basins " << i << ", " << j);
			Basin& b0 = m_basinList[i];
			Basin& b1 = m_basinList[j];
			int id0 = b0.seedId();
			int id1 = b1.seedId();

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
			float spillDist = m_maxSpillDist / std::abs(basinRaster.props().resX());

			Cell c0;
			cells0.reset();
			while(cells0.next(c0)) {
				if(cancel)
					break;
				cells1.search(c0, spillDist, std::back_inserter(result));
				for(const Cell& c : result) {
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

void Flood::findMinima() {
	g_debug("Finding minima.");

	m_seeds.clear();
	int cols = m_dem.props().cols();
	int rows = m_dem.props().rows();
	float nodata = m_dem.props().nodata();
	for (int r = 0; r < rows && !cancel; ++r) {
		for (int c = 0; c < cols; ++c) {
			bool skip = false;
			float v;
			if ((v = m_dem.get(c, r)) == nodata)
				continue;
			for (int rr = geo::max(0, r - 1); !skip && rr < geo::min(r + 2, rows); ++rr) {
				for (int cc = geo::max(0, c - 1); !skip && cc < geo::min(c + 2, cols) && !cancel; ++cc) {
					float v0;
					if ((cc == c && rr == r) || (v0 = m_dem.get(cc, rr)) == nodata)
						continue;
					if (m_dem.get(cc, rr) < m_dem.get(c, r))
						skip = true;
				}
			}
			if (!skip)
				m_seeds.emplace_back(0, c, r, m_dem.get(c, r), 0);
		}
	}
}

void Flood::flood(int numThreads) {

	using namespace geo::flood::util;

	g_debug("Flooding...");

	std::list<std::thread> threads;
	std::queue<int> elevations;

	init();

	g_debug("Building seed list...");
	if (!seedsFile().empty()) {
		// Load the seeds if given.
		loadSeeds(true);
	} else {
		// Otherwise find and use minima.
		findMinima();
	}

	for(int e = start(); e <= end(); e += step())
		elevations.push(e);

	cancel = false;

	for(int i = 0; i < numThreads; ++i)
		threads.emplace_back(worker, this, &elevations);

	for(std::thread& th : threads)
		th.join();

}


