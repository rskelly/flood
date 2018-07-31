/*
 * flood.hpp
 *
 *  Created on: Dec 24, 2017
 *      Author: rob
 */

#ifndef INCLUDE_FLOOD_HPP_
#define INCLUDE_FLOOD_HPP_


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
			static size_t _cellId = 0;

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
				size_t m_cellId;
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
				Cell(size_t id, int col, int row) :
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
				Cell(size_t id, int col, int row, double value) :
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

				size_t cellId() const {
					return m_cellId;
				}

				double value() const {
					return m_value;
				}

				/**
				 * Return the distance between cells in map units, given the resolution of the original raster.
				 */
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
				int m_band;

			public:

				Basin(uint32_t id, int minc, int minr, int maxc, int maxr, int area) :
						m_id(id),
						m_minc(minc), m_minr(minr),
						m_maxc(maxc), m_maxr(maxr),
						m_area(area),
						m_band(1) {
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
							if ((uint32_t) grd.getInt(c, r, m_band) == m_id) {
								for (int i = 0; i < 4; ++i) {
									int co = offset[i][0] + c;
									int ro = offset[i][1] + r;
									if(co < 0 || ro < 0 || co >= cols || ro >= rows ||
											(uint32_t) grd.getInt(co, ro, m_band) != m_id) {
										edgeCells.add(new Cell(m_id, c, r));
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
				int m_srcBand;
				int m_dstBand;

			public:

				LEFillOperator(Grid* src, int srcBand, Grid* dst, int dstBand, double elevation, uint32_t id) :
					m_src(src),
					m_dst(dst),
					m_elevation(elevation),
					m_id(id),
					m_nofill(0),
					m_srcBand(srcBand),
					m_dstBand(dstBand) {

					m_nodata = src->props().nodata();
				}

				bool shouldFill(int col, int row) const {
					int fill = m_dst->getInt(col, row, m_srcBand);
					if(fill == m_nofill) {
						double value = m_src->getFloat(col, row, m_srcBand);
						return value != m_nodata && value <= m_elevation;
					}
					return false;
				}

				void fill(int col, int row) const {
					m_dst->setInt(col, row, (int) m_id, m_dstBand);
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

				/**
				 * Create a spill point. Copies the given cells.
				 */
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

				/**
				 * Return the centroid of the line connecting the two cells.
				 */
				void centroid(int* col, int* row) const {
					*col = (int) (m_c1.col() + (m_c2.col() - m_c1.col()) / 2.0);
					*row = (int) (m_c1.row() + (m_c2.row() - m_c2.row()) / 2.0);
				}

				~SpillPoint() {
				}
			};

		} // util

		using namespace flood::util;

		class Flood {
		private:
			double m_start;
			double m_end;
			double m_step;
			double m_minBasinArea;
			double m_maxSpillDist;
			int m_t; // number of threads
			int m_band;
			MemRaster* m_dem;
			std::string m_input;
			std::string m_vdir;
			std::string m_rdir;
			std::string m_spill;
			std::string m_fseeds;
			std::string m_outfile;
			std::vector<Cell> m_seeds;
			std::vector<Basin> m_basinList;
			std::vector<SpillPoint> m_spillPoints;

		public:
			static bool cancel;

			/**
			 * Create a Flood object.
			 *
			 * @param input The DEM file.
			 * @param vdir The vector output directory. If not given, no vectors are produced.
			 * @param rdir The raster output directory.
			 * @param spill The spill points file. CSV or SHP.
			 * @param seeds The seeds file. CSV.
			 * @param outfile If the height flag was specified, the output file can be used to direct the single output file.
			 * @param start The starting elevation.
			 * @param end The ending elevation.
			 * @param step The step elevation.
			 * @param minBasinArea The minimum area of a region to be considered a basing for connectivity purposes.
			 * @param maxSpillDist The maximum distance between two basins before they are suspected of connecting.
			 */
			Flood(const std::string& input, const std::string& vdir, const std::string& rdir,
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

			/**
			 * Copy the given Flood object.
			 */
			Flood(const Flood& conf) :
				Flood(conf.m_input, conf.m_vdir, conf.m_rdir, conf.m_spill, conf.m_fseeds, conf.m_outfile,
					conf.m_start, conf.m_end, conf.m_step, conf.m_minBasinArea, conf.m_maxSpillDist) {
				// Copy the thread count.
				m_t = conf.m_t;
			}

			void validateInputs() {

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

			~Flood() {
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

			/**
			 * Load the seeds from the seeds file. If header is given,
			 * the first row is read and discarded.
			 */
			void loadSeeds(bool header = false) {
				if (m_fseeds.empty())
					g_argerr("No seed file given.");
				std::ifstream csv(m_fseeds);
				std::vector<std::string> row;
				if (header)
					readline(csv, row);
				const GridProps &props = m_dem->props();
				while (readline(csv, row)) {
					uint32_t id = (uint32_t) std::strtol(row[0].c_str(), nullptr, 10);
					double x = std::strtod(row[1].c_str(), nullptr);
					double y = std::strtod(row[2].c_str(), nullptr);
					g_debug("Found seed: " << id << ", " << x << ", " << y << ", " << props.toCol(x) << ", " << props.toRow(y));
					m_seeds.emplace_back(id, props.toCol(x), props.toRow(y));
				}
			}

			const std::vector<Cell>& seeds() const {
				return m_seeds;
			}

			MemRaster* dem() const {
				return m_dem;
			}

			/**
			 * Initialize the flood processor.
			 *
			 * @param mtx A mutex to protect resources.
			 * @param mapped Set to true to use file-mapped memory.
			 */
			void init(std::mutex& mtx, bool mapped) {

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

			/**
			 * Perform flood filling and identify basins.
			 * @param filename The output file for the basin raster.
			 * @param elevation The elevation to fill to.
			 * @param mapped Set to true to used mapped memory. Slow.
			 * @return An int representing the number of basins created.
			 */
			int fillBasins(std::string &filename, double elevation, bool mapped) {
				g_debug("Filling basins: " << filename << "; " << elevation);

				// Set up the output raster.
				GridProps props(m_dem->props());
				props.setBands(1);
				props.setNoData(0);
				props.setDataType(DataType::UInt32);
				props.setWritable(true);
				props.setCompress(true);
				MemRaster basins(props, mapped);
				basins.fillInt(0, m_band);
				m_basinList.clear();

				// Set up the fill operator.
				LEFillOperator op(m_dem, 1, &basins, 1, elevation, 0);
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
						TargetFillOperator<double, int> fop(m_dem, m_band, &basins, m_band, seed.cellId(), props.nodata());
						Grid::floodFill(seed.col(), seed.row(), fop, false, &minc, &minr, &maxc, &maxr, &area);
					}

				}

				g_debug("Writing basin raster " << filename);
				Raster bout(filename, props);
				for(int r = 0; r < props.rows() && !Flood::cancel; ++r)
					basins.writeTo(bout, props.cols(), 1, 0, r, 0, r);

				return m_basinList.size();
			}

			class FloodCallbacks : public geo::util::Callbacks {
			public:
				void stepCallback(float status) const {}
				void overallCallback(float status) const {}
				void statusCallback(const std::string& msg) const {}
		   };

			/**
			 * Vectorize the given raster.
			 * @param rfile the raster file.
			 * @param vfile the vector file.
			 */
			void saveBasinVector(const std::string& rfile, const std::string& vfile) {
				g_debug("Saving vector " << vfile);
				Raster rast(rfile);
				FloodCallbacks cb;
				geo::util::Status status(&cb, 0.0f, 1.0f);
				rast.polygonize(vfile, "basin", "SQLite", rast.props().hsrid(), 1, false, false, "", 0, 1, Flood::cancel, status);
			}

			bool findSpillPoints(const std::string& rfile, double elevation) {
				g_debug("Finding spill points.");

				m_spillPoints.clear();

				Raster basins(rfile);
				MemRaster buf(basins.props(), true); // TODO: Configure mapped.
				for(int r = 0; r < basins.props().rows(); ++r)
					basins.writeTo(buf, basins.props().cols(), 1, 0, r, 0, r);

				geo::util::Bounds bounds(0, 0, m_dem->props().cols(), m_dem->props().rows());

				std::unordered_map<int, KDTree<Cell> > trees;
				std::unordered_set<int> seen;

				// Compare each basin to each other basin.
				for (size_t i = 0; i < m_basinList.size(); ++i) {
					for (size_t j = i + 1; j < m_basinList.size() && !Flood::cancel; ++j) {

						g_debug("Basins " << i << ", " << j);
						Basin& b0 = m_basinList[i];
						Basin& b1 = m_basinList[j];
						uint32_t id0 = b0.id();
						uint32_t id1 = b1.id();

						if(trees.find(id0) == trees.end()) {
							trees.emplace(id0, 2);
							KDTree<Cell>& q = trees[id0];
							b0.computeEdges(buf, q);
							q.build();
						}
						if(trees.find(id1) == trees.end()) {
							trees.emplace(id1, 2);
							KDTree<Cell>& q = trees[id1];
							b1.computeEdges(buf, q);
							q.build();
						}

						int count0 = trees[id0].size();
						int count1 = trees[id1].size();

						// Compare the distances; save the ones that are near enough.
						std::list<Cell> result;
						KDTree<Cell>& cells0 = count0 < count1 ? trees[id0] : trees[id1];
						KDTree<Cell>& cells1 = count0 < count1 ? trees[id1] : trees[id0];

						for(const Cell* c0 : cells0.items()) {
							if(Flood::cancel)
								break;
							seen.insert(c0->cellId());
							std::vector<Cell*> result;
							std::vector<double> dist;
							cells1.knn(*c0, 1, std::back_inserter(result), std::back_inserter(dist));
							for(size_t i = 0; i < result.size(); ++i) {
								if(dist[i] > m_maxSpillDist)
									continue;
								seen.insert(result[i]->cellId());
								// Create a spill point; copies the cells.
								m_spillPoints.push_back(SpillPoint(*c0, *result[i], elevation));
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
							m_seeds.push_back(Cell(c, r, m_dem->getFloat(c, r, m_band)));
					}
				}
			}

			static void worker(Flood* config, std::mutex* mtx, std::ofstream* ofs, std::queue<double>* elevations, bool mapped) {

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

			/**
			 */
			void flood(int numThreads, bool memMapped) {

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
					ofs << "id,id1,x1,y1,id2,x2,y2,xmidpoint,ymidpoint,elevation,distance\n";
				}

				for(double e = start(); e <= end(); e += step())
					elevations.push(e);

				if(m_start == m_end && !m_outfile.empty()) {

					g_debug("Filling to " << m_end);

					// Spill point ID. Needed?
					uint32_t id = 0;

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

					// Generate basins.
					int basins = config.fillBasins(m_outfile, m_end, memMapped);

					// If desired, generate vectors.
					if (basins > 0 && !config.vdir().empty()) {
						std::string vec = Util::pathJoin(Util::parent(m_outfile), Util::basename(m_outfile) + ".sqlite");
						config.saveBasinVector(m_outfile, vec);
					}

					// Find and output spill points.
					if (basins > 1 && !config.spill().empty() && config.findSpillPoints(m_outfile, m_end) > 0) {
						config.saveSpillPoints(&id, ofs);
					}

				} else {
					std::mutex mtx;
					for(int i = 0; i < numThreads; ++i)
						threads.emplace_back(Flood::worker, &config, &mtx, &ofs, &elevations, memMapped);

					for(std::thread& th : threads)
						th.join();
				}
			}

		};

	} // flood

} // geo

#endif /* INCLUDE_FLOOD_HPP_ */
