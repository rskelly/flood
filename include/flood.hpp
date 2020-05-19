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

#include "ds/mqtree.hpp"
#include "grid.hpp"
#include "geo.hpp"

using namespace geo::grid;
using namespace geo::ds;

namespace geo {

	namespace flood {

		namespace util {

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
				size_t m_seedId;

			public:

				/**
				 * Construct an empty cell with an automatic ID.
				 */
				Cell();

				/**
				 * Construct a new cell with the grid coordinates and an automatic ID.
				 * \param col The column.
				 * \param row The row.
				 */
				Cell(size_t seedId, int col, int row);

				/**
				 * Construct a new cell with the grid coordinates and an explicit ID.
				 * \param id The ID.
				 * \param col The column.
				 * \param row The row.
				 */
				Cell(size_t id, size_t seedId, int col, int row);

				/**
				 * Construct a new cell with the grid coordinates and an explicit ID and value.
				 * \param id The ID.
				 * \param col The column.
				 * \param row The row.
				 * \param value The value.
				 */
				Cell(size_t id, size_t seedId, int col, int row, double value);

				/**
				 * Construct a new cell with the grid coordinates and an automatic ID and explicit value.
				 * \param col The column.
				 * \param row The row.
				 * \param value The value.
				 */
				Cell(size_t seedId, int col, int row, double value);

				/**
				 * Return the coordinate by axis index.
				 * \param idx The index
				 */
				double operator[](int idx) const;

				double x() const;

				double y() const;

				int row() const;

				int col() const;

				size_t cellId() const;

				size_t seedId() const;

				double value() const;

				/**
				 * Return the distance between cells in map units, given the resolution of the original raster.
				 */
				double distance(const Cell& other, double resx, double resy) const;

			};

			// Represents a basin, with bounds and area and the ID of the seed that spawned it.
			class Basin {
			private:
				unsigned int m_id;
				int m_minc;
				int m_minr;
				int m_maxc;
				int m_maxr;
				int m_area;
				int m_band;

			public:

				/**
				 * Construct a basin with the given seed ID, extent and area.
				 *
				 * \param seedId The seed ID.
				 * \param minc The minimum column index.
				 * \param minr The minimum row index.
				 * \param maxc The maximum column index.
				 * \param maxr The maximum row index.
				 * \param area The basin area.
				 */
				Basin(unsigned int seedId, int minc, int minr, int maxc, int maxr, int area);

				/**
				 * Make a list of all the cells that are on the edges of this basin and are
				 * also a local minimum.
				 *
				 * \param grd The flood fill grid.
				 * \param edgeCells A ,mqtree to insert the cells into.
				 * @return The number of cells found.
				 */
				int computeEdges(Band<int>& grd, mqtree<Cell>& edgeCells);

				/**
				 * Return the basins seed ID.
				 *
				 * @return The seed ID.
				 */
				unsigned int seedId() const;

			};

			// A flood fill operator that fills pixels whose values are lower than
			// the given elevation.
			class LEFillOperator: public FillOperator<float, int> {
			private:
				Band<float>* m_src;
				Band<int>* m_dst;
				double m_elevation;
				double m_nodata;
				unsigned int m_id;
				int m_nofill;
				int m_srcBand;
				int m_dstBand;

			public:

				LEFillOperator(Band<float>* src, int srcBand, Band<int>* dst, int dstBand, double elevation, unsigned int id);

				bool shouldFill(int col, int row) const;

				void fill(int col, int row) const;

				void setFill(unsigned int fill);

				void setNoFill(int nofill);

				const GridProps& srcProps() const;

            	const GridProps& dstProps() const;

            	~LEFillOperator();

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
				SpillPoint(const Cell& c1, const Cell& c2, double elevation);

				const Cell& cell1() const;

				const Cell& cell2() const;

				double elevation() const;

				/**
				 * Return the centroid of the line connecting the two cells.
				 */
				void centroid(int* col, int* row) const;

				~SpillPoint();

			};

		} // util

		using namespace flood::util;

		class Output {
		public:
			std::string rdir;
			virtual bool valid() const { return false; }
			virtual void prepare() const {};
			virtual bool doVectors() const { return false; }
			virtual void saveVector(Band<int>& basins, float elevation) const {}
			std::string rasterFile(float elevation) const {
				std::stringstream ss;
				ss << rdir << "/" << (int) (elevation * 10000) << ".tif";
				return ss.str();
			}
			virtual ~Output() {}
		};

		class DBOutput : public Output {
		public:
			std::string conn;
			std::string layer;
			std::string id;
			std::string elevation;
			bool valid() const {
				return !rdir.empty() && !conn.empty() && !layer.empty() && !id.empty() && !elevation.empty();
			}
			void prepare() const {
				if(!rdir.empty() && !isdir(rdir))
					makedir(rdir);
				if(!isdir(rdir))
					throw std::runtime_error("Raster output directory does not exist and could not be created.");
			}
			bool doVectors() const { return true; }
			void saveVector(Band<int>& rast, float elev) const {
				std::vector<PolygonValue> fields = {{elevation, elev}};
				rast.polygonizeToTable(conn, layer, id, fields);
			}
		};

		class FileOutput : public Output {
		public:
			std::string vdir;
			bool valid() const {
				return !rdir.empty();
			}
			void prepare() const {
				if(!vdir.empty() && !isdir(vdir))
					makedir(vdir);
				if(!vdir.empty() && !isdir(vdir))
					throw std::runtime_error("Vector output directory does not exist and could not be created.");
				if(!rdir.empty() && !isdir(rdir))
					makedir(rdir);
				if(!isdir(rdir))
					throw std::runtime_error("Raster output directory does not exist and could not be created.");
			}
			bool doVectors() const { return !vdir.empty(); }
			void saveVector(Band<int>& rast, float elevation) const {
				std::stringstream ss;
				ss << vdir << "/" << (int) (elevation * 10000) << ".sqlite";
				std::string vfile = ss.str();
				rast.polygonizeToFile(vfile, "basins", "bid", "SQLite");
			}
		};

		class BreakLine {
		public:
			double x0;
			double y0;
			double x1;
			double y1;
			float value;
			BreakLine(double x0, double y0, double x1, double y1, float value) :
				x0(x0), y0(y0), x1(x1), y1(y1), value(value) {}
		};

		class Flood {
		private:
			double m_start;
			double m_end;
			double m_step;
			double m_minBasinArea;
			double m_maxSpillDist;
			unsigned int m_t; // number of threads
			unsigned int m_band;
			Band<float> m_dem;
			const Output* m_output;
			std::string m_input;
			std::string m_spill;
			std::string m_fseeds;
			std::vector<Cell> m_seeds;
			std::vector<Basin> m_basinList;
			std::vector<SpillPoint> m_spillPoints;
			std::vector<BreakLine> m_breakLines;

		public:
			static bool cancel;

			/**
			 * Create a Flood object.
			 *
			 * \param input The DEM file.
			 * \param band Input raster band.
			 * \param output Configuration for file/DB/etc. output.
			 * \param spill The spill points file. CSV or SHP.
			 * \param seeds The seeds file. CSV.
			 * \param start The starting elevation.
			 * \param end The ending elevation.
			 * \param step The step elevation.
			 * \param minBasinArea The minimum area of a region to be considered a basing for connectivity purposes.
			 * \param maxSpillDist The maximum distance between two basins before they are suspected of connecting.
			 * \param breakLines A list of breaklines to apply to the input raster.
			 */
			Flood(const std::string& input, int band, const Output& output,
					const std::string& spill, const std::string& seeds,
					double start, double end, double step,
					double minBasinArea, double maxSpillDist,
					const std::vector<BreakLine>& breakLines);

			/**
			 * Copy the given Flood object.
			 */
			Flood(const Flood& conf);

			/**
			 * Check the inputs, throw an exception if any are invalid.
			 */
			void validateInputs();

			~Flood();

			const std::string& input() const;

			const Output& output() const;

			const std::string& rdir() const;

			const std::string& spill() const;

			const std::string& seedsFile() const;

			double start() const;

			double end() const;

			double step() const;

			double minBasinArea() const;

			double maxSpillDist() const;

			/**
			 * Load the seeds from the seeds file. If header is given,
			 * the first row is read and discarded.
			 */
			void loadSeeds(bool header = false);

			const std::vector<Cell>& seeds() const;

			Band<float>& dem();

			/**
			 * Initialize the flood processor.
			 *
			 * \param mtx A mutex to protect resources.
			 */
			void init(std::mutex& mtx);

			/**
			 * Perform flood filling and identify basins.
			 * \param basinRaster An empty unique_ptr to contain a Band<float>.
			 * \param elevation The elevation to fill to.
			 * \param mapped Set to true to used mapped memory. Slow.
			 * @return An int representing the number of basins created.
			 */
			int fillBasins(Band<int>& basinRaster, double elevation);

			bool findSpillPoints(Band<int>& basinRaster, double elevation);

			/**
			 * Output the spill points to a stream, with comma delimiters.
			 * The fields are: ID1, x1, y1, ID2, x2, y2, midpoint x, midpoint y, distance
			 */
			void saveSpillPoints(unsigned int* id, std::ostream &out);

			/**
			 * Find the cells at the bottoms of depressions.
			 */
			void findMinima();

			/**
			 * Called by flood, used by threads.
			 */
			static void worker(Flood* config, std::mutex* mtx, std::ofstream* ofs, std::queue<double>* elevations);

			/**
			 * Start the flood process.
			 *
			 * \param numThreads The number of threads to run.
			 */
			void flood(int numThreads);

		};

	} // flood

} // geo

#endif /* INCLUDE_FLOOD_HPP_ */
