/*
 * flood.cpp
 *
 *  Created on: Aug 2, 2018
 *      Author: rob
 */

#include "flood.hpp"

using namespace geo::flood::util;


/**
 * Maintains the IDs for cells in the raster.
 */
static size_t _cellId = 0;
static std::mutex _cellMtx;

static size_t _nextCellId() {
	std::lock_guard<std::mutex> lk(_cellMtx);
	return ++_cellId;
}


bool geo::flood::util::readline(std::istream &st, std::vector<std::string> &row) {
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

geo::flood::util::Cell::Cell() :
	m_col(0), m_row(0),
	m_cellId(_nextCellId()),
	m_value(0),
	m_seedId(0) {
}

geo::flood::util::Cell::Cell(size_t seedId, int col, int row) :
	m_col(col), m_row(row),
	m_cellId(++_cellId),
	m_value(0),
	m_seedId(seedId) {
}

geo::flood::util::Cell::Cell(size_t id, size_t seedId, int col, int row) :
	m_col(col), m_row(row),
	m_cellId(id),
	m_value(0),
	m_seedId(seedId) {
}

geo::flood::util::Cell::Cell(size_t id, size_t seedId, int col, int row, double value) :
	m_col(col), m_row(row),
	m_cellId(id),
	m_value(value),
	m_seedId(seedId) {
}

geo::flood::util::Cell::Cell(size_t seedId, int col, int row, double value) :
	m_col(col), m_row(row),
	m_cellId(++_cellId),
	m_value(value),
	m_seedId(seedId) {
}

double geo::flood::util::Cell::operator[](int idx) const {
	switch(idx % 2) {
	case 0: return m_col;
	case 1: return m_row;
	}
	return 0;
}

double geo::flood::util::Cell::x() const {
	return m_col;
}

double geo::flood::util::Cell::y() const {
	return m_row;
}

int geo::flood::util::Cell::row() const {
	return m_row;
}

int geo::flood::util::Cell::col() const {
	return m_col;
}

size_t geo::flood::util::Cell::cellId() const {
	return m_cellId;
}

size_t geo::flood::util::Cell::seedId() const {
	return m_seedId;
}

double geo::flood::util::Cell::value() const {
	return m_value;
}

double geo::flood::util::Cell::distance(const Cell& other, double resx, double resy) const {
	double x0 = col() * resx;
	double y0 = row() * resy;
	double x1 = other.col() * resx;
	double y1 = other.row() * resy;
	return std::sqrt(g_sq(x0 - x1) + g_sq(y0 - y1));
}



Basin::Basin(unsigned int id, int minc, int minr, int maxc, int maxr, int area) :
	m_id(id),
	m_minc(minc), m_minr(minr),
	m_maxc(maxc), m_maxr(maxr),
	m_area(area),
	m_band(1) {
}

int Basin::computeEdges(Grid& grd, Grid& dem, KDTree<Cell>& edgeCells) {
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
			if ((unsigned int) grd.getInt(c, r, m_band) == m_id) {
				for (int i = 0; i < 4; ++i) {
					int co = c + offset[i][1];
					int ro = r + offset[i][0];
					// Check if this is an edge pixel.
					if(co < 0 || ro < 0 || co >= cols || ro >= rows || (unsigned int) grd.getInt(co, ro, m_band) != m_id) {
						edgeCells.add(new Cell(m_id, c, r));
						out << props.toCentroidX(c) << "," << props.toCentroidY(r) << "\n";
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


LEFillOperator::LEFillOperator(Grid* src, int srcBand, Grid* dst, int dstBand, double elevation, unsigned int id) :
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
	int fill = m_dst->getInt(col, row, m_srcBand);
	if(fill == m_nofill) {
		double value = m_src->getFloat(col, row, m_srcBand);
		return value != m_nodata && value <= m_elevation;
	}
	return false;
}

void LEFillOperator::fill(int col, int row) const {
	m_dst->setInt(col, row, (int) m_id, m_dstBand);
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

const geo::flood::util::Cell& SpillPoint::cell1() const {
	return m_c1;
}

const geo::flood::util::Cell& SpillPoint::cell2() const {
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

Flood::Flood(const std::string& input, const std::string& vdir, const std::string& rdir,
		const std::string& spill, const std::string& seeds, const std::string& outfile,
		double start, double end, double step,
		double minBasinArea, double maxSpillDist) :
			m_start(start), m_end(end), m_step(step),
			m_minBasinArea(minBasinArea),
			m_maxSpillDist(maxSpillDist),
			m_t(1),
			m_band(1),
			m_dem(nullptr),
			m_input(input),
			m_vdir(vdir),
			m_rdir(rdir),
			m_spill(spill),
			m_fseeds(seeds),
			m_outfile(outfile) {
}

Flood::Flood(const Flood& conf) :
	Flood(conf.m_input, conf.m_vdir, conf.m_rdir, conf.m_spill, conf.m_fseeds, conf.m_outfile,
		conf.m_start, conf.m_end, conf.m_step, conf.m_minBasinArea, conf.m_maxSpillDist) {
	// Copy the thread count.
	m_t = conf.m_t;
}

void Flood::validateInputs() {

	g_debug("Checking...");

	if(m_start > m_end)
		throw std::runtime_error("End elevation must be larger than start.");

	if(m_step <= 0)
		throw std::runtime_error("Step must be greater than zero.");

	if(m_start != m_end && !m_outfile.empty())
		throw std::runtime_error("The output file is specified, but there will be more than one output. This is not allowed.");

	if(m_minBasinArea < 0)
		throw std::runtime_error("Minimum basin area must be greater than or equal to zero.");

	if(m_maxSpillDist <= 0)
		throw std::runtime_error("Maximum spill distance must be larger than zero.");

	if (m_vdir.empty())
		g_debug("WARNING: No vector directory; not producing flood vectors.");

	if (m_spill.empty())
		g_debug("WARNING: No spill file; not producing spill points.");

	if (std::isnan(m_start))
		g_debug("WARNING: Start value not given; using raster minimum.");

	if (std::isnan(m_end))
		g_debug("WARNING: End value not given; using raster maximum.");

	unsigned int threads = std::thread::hardware_concurrency();
	if(m_t < 1)
		throw std::runtime_error("Threads must be greater than zero.");

	if(threads > 0 && m_t > threads)
		throw std::runtime_error("Threads must be less than or equal to the number of available cores.");

	if(m_input.empty())
		throw std::runtime_error("Input raster not given.");

	if(!Util::exists(m_input))
		throw std::runtime_error("Input raster not found.");

	if(!m_vdir.empty() && !Util::exists(m_vdir))
		Util::mkdir(m_vdir);

	if(!m_vdir.empty() && !Util::exists(m_vdir))
		throw std::runtime_error("Vector output directory does not exist and could not be created.");

	if(m_outfile.empty() && m_rdir.empty())
		throw std::runtime_error("Raster output directory not given.");

	if(!m_outfile.empty()) {
		std::string dir = Util::parent(m_outfile);
		if(!Util::exists(dir))
			Util::mkdir(dir);
		if(!Util::exists(dir))
			throw std::runtime_error("Failed to create directory for output file.");
	}

	if(!m_rdir.empty()) {
		if(!Util::exists(m_rdir))
			Util::mkdir(m_rdir);
		if(!Util::exists(m_rdir))
			throw std::runtime_error("Raster output directory does not exist and could not be created.");
	}

	if(!Util::exists(m_fseeds))
		throw std::runtime_error("Seeds file does not exists.");
}

Flood::~Flood() {
	delete m_dem;
}

const std::string& Flood::input() const {
	return m_input;
}

const std::string& Flood::vdir() const {
	return m_vdir;
}

const std::string& Flood::rdir() const {
	return m_rdir;
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
		geo::flood::util::readline(csv, row);
	const GridProps &props = m_dem->props();
	while (geo::flood::util::readline(csv, row)) {
		unsigned int id = (unsigned int) std::strtol(row[0].c_str(), nullptr, 10);
		double x = std::strtod(row[1].c_str(), nullptr);
		double y = std::strtod(row[2].c_str(), nullptr);
		g_debug("Found seed: " << id << ", " << x << ", " << y << ", " << props.toCol(x) << ", " << props.toRow(y));
		m_seeds.emplace_back(id, props.toCol(x), props.toRow(y));
	}
}

const std::vector<geo::flood::util::Cell>& Flood::seeds() const {
	return m_seeds;
}

MemRaster* Flood::dem() const {
	return m_dem;
}

void Flood::init(std::mutex& mtx, bool mapped) {

	validateInputs();

	g_debug("Initing...");
	{
		std::lock_guard<std::mutex> lk(mtx);
		Raster tmp(m_input);
		GridProps dp(tmp.props());
		m_dem = new MemRaster(dp, mapped);
		for(int r = 0; r < dp.rows(); ++r)
			tmp.writeTo(*m_dem, dp.cols(), 1, 0, r, 0, r);
	}

	g_debug("Computing stats...");
	const GridStats stats = m_dem->stats(1);
	if (std::isnan(m_start)) {
		m_start = stats.min;
	} else {
		m_start = g_max(stats.min, m_start);
	}
	if (std::isnan(m_end)) {
		m_end = stats.max;
	} else {
		m_end = g_min(stats.max, m_end);
	}
	g_debug("Stats: " << stats.min << ", " << stats.max);
	if (m_end < m_start)
		g_argerr("The ending elevation must be larger than the starting elevation: " << m_start << ", " << m_end);
}

int Flood::fillBasins(std::unique_ptr<MemRaster>& basinRaster, double elevation, bool mapped) {
	g_debug("Filling basins: " << elevation);

	// Set up the output raster.
	GridProps props(m_dem->props());
	props.setBands(1);
	props.setNoData(0);
	props.setDataType(DataType::UInt32);
	props.setWritable(true);
	props.setCompress(true);
	basinRaster.reset(new MemRaster(props, mapped));
	basinRaster->fillInt(0, m_band);
	m_basinList.clear();

	// Set up the fill operator.
	LEFillOperator op(m_dem, 1, basinRaster.get(), 1, elevation, 0);
	op.setNoFill(0);

	// Iterate over the seeds.
	for (const Cell& seed : seeds()) {

		if(Flood::cancel)
			break;

		if (!props.hasCell(seed.col(), seed.row())) {
			g_warn("Found a seed out of bounds: " << seed.col() << ", " << seed.row());
			continue;
		}

		// Fill the basin based on the elevations in the DEM.
		int minc, minr, maxc, maxr, area;
		op.setFill(seed.cellId());
		Grid::floodFill(seed.col(), seed.row(), op, false, &minc, &minr, &maxc, &maxr, &area);

		if (area >= minBasinArea()) {
			// If it's large enough, save the basin.
			g_debug("Good Seed: " << seed.col() << ", " << seed.row() << "; elevation: " << elevation << "; basin area: " << area);
			m_basinList.emplace_back(seed.cellId(), minc, minr, maxc, maxr, area);
		} else if(area > 0){
			g_debug("Bad Seed: " << seed.col() << ", " << seed.row() << "; elevation: " << elevation << "; basin area: " << area);
			// If the basin is too small, fill it with nodata. Do not collect more spill points.
			TargetFillOperator<double, int> fop(m_dem, m_band, basinRaster.get(), m_band, seed.cellId(), props.nodata());
			Grid::floodFill(seed.col(), seed.row(), fop, false, &minc, &minr, &maxc, &maxr, &area);
		}

	}

	return m_basinList.size();
}

void Flood::saveBasinVector(const std::string& rfile, const std::string& vfile) {
	g_debug("Saving vector " << vfile);
	Raster rast(rfile);
	FloodCallbacks cb;
	geo::util::Status status(&cb, 0.0f, 1.0f);
	rast.polygonize(vfile, "basin", "SQLite", rast.props().hsrid(), 1, false, false, "", 0, 1, Flood::cancel, status);
}

bool Flood::findSpillPoints(std::unique_ptr<MemRaster>& basinRaster, double elevation, bool mapped) {
	g_debug("Finding spill points.");

	m_spillPoints.clear();

	geo::util::Bounds bounds(0, 0, m_dem->props().cols(), m_dem->props().rows());

	std::unordered_map<int, KDTree<Cell> > trees;
	std::unordered_set<int> seen;

	// Compare each basin to each other basin.
	for (size_t i = 0; i < m_basinList.size(); ++i) {
		for (size_t j = i + 1; j < m_basinList.size() && !Flood::cancel; ++j) {

			g_debug("Basins " << i << ", " << j);
			Basin& b0 = m_basinList[i];
			Basin& b1 = m_basinList[j];
			unsigned int id0 = b0.seedId();
			unsigned int id1 = b1.seedId();

			if(trees.find(id0) == trees.end()) {
				trees.emplace(id0, 2);
				KDTree<Cell>& q = trees[id0];
				b0.computeEdges(*basinRaster, *m_dem, q);
				q.build();
			}
			if(trees.find(id1) == trees.end()) {
				trees.emplace(id1, 2);
				KDTree<Cell>& q = trees[id1];
				b1.computeEdges(*basinRaster, *m_dem, q);
				q.build();
			}

			int count0 = trees[id0].size();
			int count1 = trees[id1].size();

			// Compare the distances; save the ones that are near enough.
			KDTree<Cell>& cells0 = count0 < count1 ? trees[id0] : trees[id1];
			KDTree<Cell>& cells1 = count0 < count1 ? trees[id1] : trees[id0];
			std::vector<Cell*> result;
			std::vector<double> dist;
			double spillDist = m_maxSpillDist / std::abs(basinRaster->props().resolutionX());

			for(const Cell* c0 : cells0.items()) {
				if(Flood::cancel)
					break;
				seen.insert(c0->cellId());
				cells1.knn(*c0, 1, std::back_inserter(result), std::back_inserter(dist));
				for(size_t i = 0; i < result.size(); ++i) {
					std::cerr << dist[i] << "\n";
					if(dist[i] <= spillDist) {
						seen.insert(result[i]->cellId());
						// Create a spill point; copies the cells.
						m_spillPoints.emplace_back(*c0, *result[i], elevation);
					}
				}
				result.resize(0);
				dist.resize(0);
			}
		}
	}

	//g_debug("Minimum distance: " << minDist);
	return m_spillPoints.size();
}

void Flood::saveSpillPoints(unsigned int* id, std::ostream &out) {
	g_debug("Outputting spill points.");
	out << std::setprecision(12);
	const GridProps& props = m_dem->props();
	for (const SpillPoint& sp : m_spillPoints) {
		const Cell& c1 = sp.cell1();
		const Cell& c2 = sp.cell2();
		double x1 = props.toCentroidX(c1.col());
		double y1 = props.toCentroidY(c1.row());
		double x2 = props.toCentroidX(c2.col());
		double y2 = props.toCentroidY(c2.row());
		double x3 = (x1 + x2) / 2.0;
		double y3 = (y1 + y2) / 2.0;
		double dist = std::sqrt(g_sq(x1 - x2) + g_sq(y1 - y2));
		out << ++(*id) << ", " << c1.seedId() << "," << x1 << "," << y1 << "," << c2.seedId() << ","
				<< x2 << "," << y2 << "," << x3 << "," << y3 << ","
				<< sp.elevation() << "," << dist << std::endl;
	}
}

void Flood::findMinima() {
	g_debug("Finding minima.");

	m_seeds.clear();
	int cols = m_dem->props().cols();
	int rows = m_dem->props().rows();
	double nodata = m_dem->props().nodata();
	for (int r = 0; r < rows && !Flood::cancel; ++r) {
		for (int c = 0; c < cols; ++c) {
			bool skip = false;
			double v;
			if ((v = m_dem->getFloat(c, r, m_band)) == nodata)
				continue;
			for (int rr = g_max(0, r - 1); !skip && rr < g_min(r + 2, rows); ++rr) {
				for (int cc = g_max(0, c - 1); !skip && cc < g_min(c + 2, cols) && !Flood::cancel; ++cc) {
					double v0;
					if ((cc == c && rr == r) || (v0 = m_dem->getFloat(cc, rr, m_band)) == nodata)
						continue;
					if (m_dem->getFloat(cc, rr, m_band) < m_dem->getFloat(c, r, m_band))
						skip = true;
				}
			}
			if (!skip)
				m_seeds.push_back(Cell(0, c, r, m_dem->getFloat(c, r, m_band)));
		}
	}
}

void Flood::worker(Flood* config, std::mutex* mtx, std::ofstream* ofs, std::queue<double>* elevations, bool mapped) {

	Flood conf(*config);
	conf.init(*mtx, mapped);

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
		// To get the number of places.
		double exp = std::pow(10, std::floor(-std::log10(conf.step())) + 1.0);

		// The basin filename.
		std::stringstream ss;
		ss << conf.rdir() << "/" << (int) std::round(elevation * exp) << ".tif";
		std::string rfile = ss.str();

		std::unique_ptr<MemRaster> basinRaster;

		// Generate basins.
		int basins = conf.fillBasins(basinRaster, elevation, mapped);

		// Find and output spill points.
		if (basins > 1 && !conf.spill().empty() && conf.findSpillPoints(basinRaster, elevation, mapped) > 0) {
			std::lock_guard<std::mutex> lk(*mtx);
			conf.saveSpillPoints(&id, *ofs);
		}

		g_debug("Writing basin raster " << rfile);

		const GridProps& props = basinRaster->props();
		Raster bout(rfile, props);
		for(int r = 0; r < props.rows() && !Flood::cancel; ++r)
			basinRaster->writeTo(bout, props.cols(), 1, 0, r, 0, r);

		// If desired, generate vectors.
		if (basins > 0 && !conf.vdir().empty()) {
			ss.str(std::string());
			ss << conf.vdir() << "/" << (int) std::round(elevation * exp) << ".sqlite";
			std::string vfile = ss.str();
			conf.saveBasinVector(rfile, vfile);
		}
	}

}

void Flood::flood(int numThreads, bool memMapped) {

	using namespace geo::flood::util;

	g_debug("Flooding...");

	std::list<std::thread> threads;
	std::ofstream ofs;
	std::queue<double> elevations;

	// Build the config object.
	Flood config(*this);
	config.validateInputs();

	if(!config.spill().empty()) {
		Util::rm(spill());
		ofs.open(spill());
		ofs << "id,id1,x1,y1,id2,x2,y2,xmidpoint,ymidpoint,elevation,distance" << std::endl;
	}

	for(double e = start(); e <= end(); e += step())
		elevations.push(e);

	if(m_start == m_end && !m_outfile.empty()) {

		g_debug("Filling to " << m_end);

		// Spill point ID. Needed?
		unsigned int id = 0;

		std::mutex mtx;
		config.init(mtx, memMapped);

		g_debug("Building seed list...");
		if (!config.seedsFile().empty()) {
			// Load the seeds if given.
			config.loadSeeds(true);
		} else {
			// Otherwise find and use minima.
			config.findMinima();
		}

		std::unique_ptr<MemRaster> basinRaster;

		// Generate basins.
		int basins = config.fillBasins(basinRaster, m_end, memMapped);

		// Find and output spill points.
		if (basins > 1 && !config.spill().empty() && config.findSpillPoints(basinRaster, m_end, memMapped) > 0) {
			config.saveSpillPoints(&id, ofs);
		}

		const GridProps& props = basinRaster->props();
		Raster bout(m_outfile, props);
		for(int r = 0; r < props.rows() && !Flood::cancel; ++r)
			basinRaster->writeTo(bout, props.cols(), 1, 0, r, 0, r);

		// If desired, generate vectors.
		if (basins > 0 && !config.vdir().empty()) {
			std::string vec = Util::pathJoin(Util::parent(m_outfile), Util::basename(m_outfile) + ".sqlite");
			config.saveBasinVector(m_outfile, vec);
		}

	} else {
		std::mutex mtx;
		for(int i = 0; i < numThreads; ++i)
			threads.emplace_back(Flood::worker, &config, &mtx, &ofs, &elevations, memMapped);

		for(std::thread& th : threads)
			th.join();
	}
}


