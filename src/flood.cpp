/*****************************************************************************************************
 * This program 'floods' a digital terrain model using the flood fill algorithm, generating raster
 * and vector maps of the basin extents and 'spill points' indicating locations where basins can 
 * be expected to connect. The user can provide seed points to the program, or the program
 * can generate random seeds. The user can provide start and stop elevations, otherwise the program
 * will start at the elevation of the lowest seed, and stop at the elevation above the highest
 * seed when all basins have converged to one.
 *
 * The output of this program can have many uses, such as informing the selection of outlets
 * for watershed mapping with r.water.outlet.
 *
 * Run the program with no arguments for instructions.
 *****************************************************************************************************/

#include <vector>
#include <ostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <sstream>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <thread>
#include <functional>
#include <mutex>

#include "ds/kdtree.hpp"
#include "raster.hpp"
#include "geo.hpp"

using namespace geo::raster;
using namespace geo::ds;

namespace geo {

	namespace flood {

		namespace util {

			/**
			 * Maintains the IDs for cells in the raster.
			 */
			static uint64_t _cellId = 0;

			/**
			 * Reads a line from a CSV stream into a vector of values.
			 * @param st The input stream.
			 * @param row The list to append rows to.
			 */
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

			/**
			 * Represents a grid cell.
			 * Contains the grid coordinates, an ID and a value, such as elevation.
			 */
			class Cell {
			private:
				int m_col;
				int m_row;
				uint64_t m_cellId;
				double m_value;

			public:

				/**
				 * Construct an empty cell with an automatic ID.
				 */
				Cell() :
						m_col(0), m_row(0),
						m_cellId(++_cellId),
						m_value(0) {
				}

				/**
				 * Construct a new cell with the grid coordinates and an automatic ID.
				 * @param col The column.
				 * @param row The row.
				 */
				Cell(int col, int row) :
						m_col(col), m_row(row),
						m_cellId(++_cellId),
						m_value(0) {
				}

				/**
				 * Construct a new cell with the grid coordinates and an explicit ID.
				 * @param id The ID.
				 * @param col The column.
				 * @param row The row.
				 */
				Cell(uint64_t id, int col, int row) :
						m_col(col), m_row(row),
						m_cellId(id),
						m_value(0) {
				}

				/**
				 * Construct a new cell with the grid coordinates and an explicit ID and value.
				 * @param id The ID.
				 * @param col The column.
				 * @param row The row.
				 * @param value The value.
				 */
				Cell(uint64_t id, int col, int row, double value) :
						m_col(col), m_row(row),
						m_cellId(id),
						m_value(value) {
				}

				/**
				 * Construct a new cell with the grid coordinates and an automatic ID and explicit value.
				 * @param col The column.
				 * @param row The row.
				 * @param value The value.
				 */
				Cell(int col, int row, double value) :
						m_col(col), m_row(row), 
						m_cellId(++_cellId),
						m_value(value) {
				}

				/**
				 * Return the coordinate by axis index.
				 * @param idx The index
				 */
				double operator[](int idx) const {
					switch(idx % 3) {
					case 0: return m_col;
					case 1: return m_row;
					default: return 0;
					}
				}

				double x() const {
					return m_col;
				}

				double y() const {
					return m_row;
				}

				int row() const {
					return m_row;
				}

				int col() const {
					return m_col;
				}

				uint64_t cellId() const {
					return m_cellId;
				}

				double value() const {
					return m_value;
				}

				double distance(const Cell& other, double resx, double resy) const {
					double x0 = col() * resx;
					double y0 = row() * resy;
					double x1 = other.col() * resx;
					double y1 = other.row() * resy;
					return std::sqrt(g_sq(x0 - x1) + g_sq(y0 - y1));
				}

			};

			// Represents a basin, with bounds and area and the ID of the seed that spawned it.
			class Basin {
			private:
				uint32_t m_id;
				int m_minc;
				int m_minr;
				int m_maxc;
				int m_maxr;
				int m_area;

			public:

				Basin(uint32_t id, int minc, int minr, int maxc, int maxr, int area) :
						m_id(id),
						m_minc(minc), m_minr(minr),
						m_maxc(maxc), m_maxr(maxr),
						m_area(area) {
				}

				// Make a list of all the cells that are on the edges of this basin.
				int computeEdges(Grid &grd, KDTree<Cell>& edgeCells) {
					int offset[4][2] = {{0, -1}, {-1, 0}, {1, 0}, {0, 1}};
					const GridProps& props = grd.props();
					int cols = props.cols();
					int rows = props.rows();
					int count = 0;
					for (int r = m_minr; r <= m_maxr; ++r) {
						for (int c = m_minc; c <= m_maxc; ++c) {
							if ((uint32_t) grd.getInt(c, r) == m_id) {
								for (int i = 0; i < 4; ++i) {
									int co = offset[i][0] + c;
									int ro = offset[i][1] + r;
									if(co < 0 || ro < 0 || co >= cols || ro >= rows || 
											(uint32_t) grd.getInt(co, ro) != m_id) {
										Cell cl(m_id, c, r);
										edgeCells.add(cl);
										++count;
										break;
									}
								}
							}
						}
					}
					return count;
				}

				uint32_t id() const {
					return m_id;
				}

			};

			// A flood fill operator that fills pixels whose values are lower than
			// the given elevation.
			class LEFillOperator: public FillOperator<double, int> {
			private:
				Grid* m_src;
				Grid* m_dst;
				double m_elevation;
				double m_nodata;
				uint32_t m_id;
				int m_nofill;

			public:

				LEFillOperator(Grid* src, Grid* dst, double elevation, uint32_t id) :
					m_src(src), 
					m_dst(dst),
					m_elevation(elevation), 
					m_id(id),
					m_nofill(0) {
					m_nodata = src->props().nodata();
				}

				bool shouldFill(int col, int row) const {
					int fill = m_dst->getInt(col, row);
					if(fill == m_nofill) {
						double value = m_src->getFloat(col, row);
						return value != m_nodata && value <= m_elevation;
					}
					return false;
				}

				void fill(int col, int row) const {
					m_dst->setInt(col, row, (int) m_id);
				}

				void setFill(uint32_t fill) {
					m_id = fill;
				}

				void setNoFill(int nofill) {
					m_nofill = nofill;
				}

				const GridProps& srcProps() const {
                	return m_src->props();
            	}
            	
            	const GridProps& dstProps() const {
                	return m_dst->props();
            	}

            	~LEFillOperator(){}

			};

			// Represents a spill point between two cells.

			class SpillPoint {
			private:
				Cell m_c1;
				Cell m_c2;
				double m_elevation;
			public:

				SpillPoint(const Cell& c1, const Cell& c2, double elevation) :
						m_c1(c1), 
						m_c2(c2), 
						m_elevation(elevation) {
				}

				const Cell& cell1() const {
					return m_c1;
				}

				const Cell& cell2() const {
					return m_c2;
				}

				double elevation() const {
					return m_elevation;
				}

				void centroid(int* col, int* row) const {
					*col = (int) (m_c1.col() + (m_c2.col() - m_c1.col()) / 2.0);
					*row = (int) (m_c1.row() + (m_c2.row() - m_c2.row()) / 2.0);
				}

				~SpillPoint() {
				}
			};

			class Config {
			private:
				double m_start;
				double m_end;
				double m_step;
				double m_minBasinArea;
				double m_maxSpillDist;
				int m_t; // number of threads
				MemRaster* m_dem;
				std::string m_input;
				std::string m_vdir;
				std::string m_rdir;
				std::string m_spill;
				std::string m_fseeds;
				std::vector<Cell> m_seeds;
				std::vector<Basin> m_basinList;
				std::vector<SpillPoint> m_spillPoints;

			public:
				Config(const std::string& input, const std::string& vdir, const std::string& rdir,
						const std::string& spill, const std::string& seeds,
						double start, double end, double step,
						double minBasinArea, double maxSpillDist) :
							m_start(start), m_end(end), m_step(step),
							m_minBasinArea(minBasinArea),
							m_maxSpillDist(maxSpillDist),
							m_t(1),
							m_dem(nullptr),
							m_input(input),
							m_vdir(vdir),
							m_rdir(rdir),
							m_spill(spill),
							m_fseeds(seeds) {

				}

				Config(const Config& conf) :
					Config(conf.m_input, conf.m_vdir, conf.m_rdir,
							conf.m_spill, conf.m_fseeds, conf.m_start, conf.m_end, conf.m_step,
							conf.m_minBasinArea, conf.m_maxSpillDist) {
					// Copy the thread count.
					m_t = conf.m_t;
				}

				~Config() {
					delete m_dem;
				}

				const std::string& input() const {
					return m_input;
				}

				const std::string& vdir() const {
					return m_vdir;
				}

				const std::string& rdir() const {
					return m_rdir;
				}

				const std::string& spill() const {
					return m_spill;
				}

				const std::string& seedsFile() const {
					return m_fseeds;
				}

				double start() const {
					return m_start;
				}

				double end() const {
					return m_end;
				}

				double step() const {
					return m_step;
				}

				double minBasinArea() const {
					return m_minBasinArea;
				}

				double maxSpillDist() const {
					return m_maxSpillDist;
				}

				void loadSeeds(bool header = false) {
					if (m_fseeds.empty())
						g_argerr("No seed file given.");
					std::ifstream csv(m_fseeds);
					std::vector<std::string> row;
					if (header)
						readline(csv, row);
					const GridProps &props = m_dem->props();
					while (readline(csv, row)) {
						uint32_t id = (uint32_t) atoi(row[0].c_str());
						double x = atof(row[1].c_str());
						double y = atof(row[2].c_str());
						g_debug("Found seed: " << id << ", " << x << ", " << y << ", " << props.toCol(x) << ", " << props.toRow(y));
						m_seeds.push_back(Cell(id, props.toCol(x), props.toRow(y)));
					}
				}

				const std::vector<Cell>& seeds() const {
					return m_seeds;
				}

				MemRaster* dem() const {
					return m_dem;
				}

				void init(std::mutex& mtx, bool mapped) {
					g_debug("Checking...");
					if (m_input.empty())
						g_argerr("Input DEM must be provided.");
					if (m_rdir.empty())
						g_argerr("No raster directory. It is required.");
					if (m_vdir.empty())
						g_debug("WARNING: No vector directory; not producing flood vectors.");
					if (m_spill.empty())
						g_debug("WARNING: No spill file; not producing spill points.");
					if (std::isnan(m_start))
						g_debug("WARNING: Start value not given; using raster minimum.");
					if (std::isnan(m_end))
						g_debug("WARNING: End value not given; using raster maximum.");
					if (m_step <= 0.0)
						g_argerr("The step elevation must be greater than zero.");
					if (m_t < 1)
						g_debug("WARNING: Invalid number of threads. Using 1.");
					if (m_minBasinArea <= 0.0)
						g_argerr("Min basin area must be greater than zero.");
					if (m_maxSpillDist <= 0.0)
						g_argerr("Max spill distance must be greater than zero.");

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
					if (m_end <= m_start)
						g_argerr("The ending elevation must be larger than the starting elevation: " << m_start << ", " << m_end);
				}

				/**
				 * Perform flood filling and identify basins.
				 * @param filename The output file for the basin raster.
				 * @param elevation The elevation to fill to.
				 * @param mapped Set to true to used mapped memory. Slow.
				 * @return An int representing the number of basins created.
				 */
				int fillBasins(std::string &filename, double elevation, bool mapped) {
					g_debug("Filling basins: " << filename << "; " << elevation);

					GridProps props(m_dem->props());
					props.setBands(1);
					props.setNoData(0);
					props.setDataType(DataType::UInt32);
					props.setWritable(true);
					MemRaster basins(props, mapped);
					basins.fillInt(0);
					m_basinList.clear();

					LEFillOperator op(m_dem, &basins, elevation, 0);
					op.setNoFill(0);

					for (const Cell& seed : seeds()) {

						if (!props.hasCell(seed.col(), seed.row())) {
							g_warn("Found a seed out of bounds: " << seed.col() << ", " << seed.row());
							continue;
						}

						// Fill the basin based on the elevations in dem.
						int minc, minr, maxc, maxr, area;
						op.setFill(seed.cellId());
						Grid::floodFill(seed.col(), seed.row(), op, false, &minc, &minr, &maxc, &maxr, &area);

						if (area >= minBasinArea()) {
							// If it's large enough, save the basin.
							g_debug("Good Seed: " << seed.col() << ", " << seed.row() << "; elevation: " << elevation << "; basin area: " << area);
							m_basinList.push_back(Basin(seed.cellId(), minc, minr, maxc, maxr, area));
						} else if(area > 0){
							g_debug("Bad Seed: " << seed.col() << ", " << seed.row() << "; elevation: " << elevation << "; basin area: " << area);
							// If the basin is too small, fill it with nodata. Do not collect more spill points.
							TargetFillOperator<double, int> fop(m_dem, &basins, seed.cellId(), props.nodata());
							Grid::floodFill(seed.col(), seed.row(), fop, false, &minc, &minr, &maxc, &maxr, &area);
						}

					}

					g_debug("Writing basin raster " << filename);
					Raster bout(filename, props);
					for(int r = 0; r < props.rows(); ++r)
						basins.writeTo(bout, props.cols(), 1, 0, r, 0, r);

					return m_basinList.size();
				}

				/**
				 * Vectorize the given raster.
				 * @param rfile the raster file.
				 * @param vfile the vector file.
				 */
				void saveBasinVector(const std::string& rfile, const std::string& vfile) {
					g_debug("Saving vector " << vfile);
					Raster rast(rfile);
					rast.polygonize(vfile, "basin", "SQLite");
				}

				bool findSpillPoints(const std::string& rfile, double elevation) {
					g_debug("Finding spill points.");

					m_spillPoints.clear();

					Raster basins(rfile);
					MemRaster buf(basins.props(), false);
					for(int r = 0; r < basins.props().rows(); ++r)
						basins.writeTo(buf, basins.props().cols(), 1, 0, r, 0, r);

					geo::util::Bounds bounds(0, 0, m_dem->props().cols(), m_dem->props().rows());

					std::unordered_map<int, std::unique_ptr<KDTree<Cell> > > trees;
					std::unordered_set<int> seen;

					// Compare each basin to each other basin.
					for (size_t i = 0; i < m_basinList.size(); ++i) {
						for (size_t j = i + 1; j < m_basinList.size(); ++j) {

							g_debug("Basins " << i << ", " << j);
							Basin& b0 = m_basinList[i];
							Basin& b1 = m_basinList[j];
							uint32_t id0 = b0.id();
							uint32_t id1 = b1.id();

							int count0 = 0;
							int count1 = 0;

							if(trees.find(id0) == trees.end()) {
								std::unique_ptr<KDTree<Cell> > q(new KDTree<Cell>(2));
								b0.computeEdges(buf, *q);
								q->build();
								count0 = q->size();
								trees[id0] = std::move(q);
							}
							if(trees.find(id1) == trees.end()) {
								std::unique_ptr<KDTree<Cell> > q(new KDTree<Cell>(2));
								b1.computeEdges(buf, *q);
								q->build();
								count1 = q->size();
								trees[id1] = std::move(q);
							}

							// Compare the distances; save the ones that are near enough.
							std::list<Cell> result;
							std::unique_ptr<KDTree<Cell> >& cells0 = count0 < count1 ? trees[id0] : trees[id1];
							std::unique_ptr<KDTree<Cell> >& cells1 = count0 < count1 ? trees[id1] : trees[id0];

							for(const Cell& c0 : cells0->items()) {
							//while(cells0->next(c0)) {
								//if(seen.find(c0->cellId()) != seen.end())
								//	continue;
								seen.insert(c0.cellId());
								std::vector<Cell> result;
								std::vector<double> dist;
								cells1->knn(c0, 1, std::back_inserter(result), std::back_inserter(dist));
								for(size_t i = 0; i < result.size(); ++i) {
									if(dist[i] > m_maxSpillDist)
										continue;
									//if(seen.find(c1->cellId()) != seen.end())
									//	continue;
									seen.insert(result[i].cellId());
									m_spillPoints.push_back(SpillPoint(c0, result[i], elevation));
								}
							}
						}
					}

					//g_debug("Minimum distance: " << minDist);
					return m_spillPoints.size();
				}

				// Output the spill points to a stream, with comma delimiters.
				// The fields are: ID1, x1, y1, ID2, x2, y2, midpoint x, midpoint y, distance

				void saveSpillPoints(uint32_t* id, std::ostream &out) {
					g_debug("Outputting spill points.");
					out << std::setprecision(12);
					double resX = m_dem->props().resolutionX();
					double resY = m_dem->props().resolutionY();
					const Bounds& bounds = m_dem->props().bounds();
					for (const SpillPoint& sp : m_spillPoints) {
						const Cell& c1 = sp.cell1();
						const Cell& c2 = sp.cell2();
						double x1 = c1.col() * resX + bounds.minx();
						double y1 = c1.row() * resY + bounds.maxy();
						double x2 = c2.col() * resX + bounds.minx();
						double y2 = c2.row() * resY + bounds.maxy();
						double x3 = (x1 + x2) / 2.0;
						double y3 = (y1 + y2) / 2.0;
						double dist = std::sqrt(g_sq(x1 - x2) + g_sq(y1 - y2));
						out << ++(*id) << ", " << c1.cellId() << "," << x1 << "," << y1 << "," << c2.cellId() << ","
								<< x2 << "," << y2 << "," << x3 << "," << y3 << ","
								<< sp.elevation() << "," << dist << "\n";
					}
				}

				/**
				 * Find the cells at the bottoms of depressions.
				 */
				void findMinima() {
					g_debug("Finding minima.");

					m_seeds.clear();
					int cols = m_dem->props().cols();
					int rows = m_dem->props().rows();
					double nodata = m_dem->props().nodata();
					for (int r = 0; r < rows; ++r) {
						for (int c = 0; c < cols; ++c) {
							bool skip = false;
							double v;
							if ((v = m_dem->getFloat(c, r)) == nodata)
								continue;
							for (int rr = g_max(0, r - 1);
									!skip && rr < g_min(r + 2, rows); ++rr) {
								for (int cc = g_max(0, c - 1);
										!skip && cc < g_min(c + 2, cols); ++cc) {
									double v0;
									if ((cc == c && rr == r) || (v0 = m_dem->getFloat(cc, rr)) == nodata)
										continue;
									if (m_dem->getFloat(cc, rr) < m_dem->getFloat(c, r))
										skip = true;
								}
							}
							if (!skip)
								m_seeds.push_back(Cell(c, r, m_dem->getFloat(c, r)));
						}
					}
				}
			};

			void worker(Config* config, std::mutex* mtx, std::ofstream* ofs, std::queue<double>* elevations, bool mapped) {

				Config conf(*config);
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

				while(run) {
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
					uint32_t id = 0;
					// To get the number of places.
					double exp = std::pow(10, std::floor(-std::log10(conf.step())) + 1.0);

					// The basin filename.
					std::stringstream ss;
					ss << conf.rdir() << "/" << (int) std::round(elevation * exp) << ".tif";
					std::string rfile = ss.str();

					// Generate basins.
					int basins = conf.fillBasins(rfile, elevation, mapped);

					// If desired, generate vectors.
					if (basins > 0 && !conf.vdir().empty()) {
						ss.str(std::string());
						ss << conf.vdir() << "/" << (int) std::round(elevation * exp) << ".sqlite";
						std::string vfile = ss.str();
						conf.saveBasinVector(rfile, vfile);
					}

					// Find and output spill points.
					if (basins > 1 && !conf.spill().empty() && conf.findSpillPoints(rfile, elevation) > 0) {
						std::lock_guard<std::mutex> lk(*mtx);
						conf.saveSpillPoints(&id, *ofs);
					}

				}
			}

		} // util

		/**
		 */
		void flood(std::string& input, std::string& seeds, std::string& vdir,
				std::string& rdir, std::string& spill, double start, double end,
				double step, double minBasinArea, double maxSpillDist, int t, bool mapped) {

			using namespace geo::flood::util;

			g_debug("Flooding...");

			std::list<std::thread> threads;
			std::ofstream ofs;
			std::queue<double> elevations;

			// Build the config object.
			Config config(input, vdir, rdir, spill, seeds, start, end, step, minBasinArea, maxSpillDist);

			if(!config.spill().empty()) {
				Util::rm(spill);
				ofs.open(spill);
				ofs << "id,id1,x1,y1,id2,x2,y2,xmidpoint,ymidpoint,elevation,distance\n";
			}

			for(double e = start; e <= end; e += step)
				elevations.push(e);

			std::mutex mtx;
			for(int i = 0; i < t; ++i)
				threads.push_back(std::thread(worker, &config, &mtx, &ofs, &elevations, mapped));

			for(std::thread& th : threads)
				th.join();
		}

	} // flood

} // geo

void usage() {
	std::cerr << "Usage: flood <options>\n"
			<< " -i     Input elevation raster.\n"
			<< " -s     CSV containing seed points. Row format, ID, X, Y. If not given, minima are used.\n"
			<< " -v     Directory for basin vectors. If not given, they are not produced.\n"
			<< " -r     Directory for basin rasters. If not given, they are produced.\n"
			<< " -p     Spill point file. If not given, they are not produced. A Shapefile or CSV.\n"
			<< " -start Starting elevation, in same units as input raster.\n"
			<< " -end   Ending elevation.\n" << " -step  Step elevation.\n"
			<< " -t     Number of threads to use. Default 1.\n"
			<< " -b     Minimum basin area.\n"
			<< " -d     Maximum spill distance.\n";
}

int main(int argc, char **argv) {

	g_loglevel (G_LOG_DEBUG);

	std::string input;
	std::string seeds;
	std::string vdir;
	std::string rdir;
	std::string spill;
	double start = 0.0;
	double end = 0.0;
	double step = 0.0;
	double maxSpillDist = 100.0;
	double minBasinArea = 100.0;
	int t = 1;

	for (int i = 1; i < argc; ++i) {
		std::string a(argv[i]);
		if (a == "-i") {
			input.assign(argv[++i]);
		} else if (a == "-s") {
			seeds.assign(argv[++i]);
		} else if (a == "-v") {
			vdir.assign(argv[++i]);
		} else if (a == "-r") {
			rdir.assign(argv[++i]);
		} else if (a == "-p") {
			spill.assign(argv[++i]);
		} else if (a == "-start") {
			start = atof(argv[++i]);
		} else if (a == "-end") {
			end = atof(argv[++i]);
		} else if (a == "-step") {
			step = atof(argv[++i]);
		} else if (a == "-t") {
			t = atoi(argv[++i]);
		} else if (a == "-b") {
			minBasinArea = atof(argv[++i]);
		} else if (a == "-d") {
			maxSpillDist = atof(argv[++i]);
		}
	}

	try {

		geo::flood::flood(input, seeds, vdir, rdir, spill, start, end,
				step, minBasinArea, maxSpillDist, t, false);

	} catch (const std::exception &e) {
		std::cerr << e.what() << "\n";
		usage();
		return 1;
	}

	return 0;
}


