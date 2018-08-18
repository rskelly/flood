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
#include "db.hpp"
#include "geo.hpp"

using namespace geo::raster;
using namespace geo::ds;
using namespace geo::db;

namespace geo {

	namespace flood {

		namespace util {

			/**
			 * Reads a line from a CSV stream into a vector of values.
			 * @param st The input stream.
			 * @param row The list to append rows to.
			 */
			bool readline(std::istream &st, std::vector<std::string> &row);

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
				 * @param col The column.
				 * @param row The row.
				 */
				Cell(size_t seedId, int col, int row);

				/**
				 * Construct a new cell with the grid coordinates and an explicit ID.
				 * @param id The ID.
				 * @param col The column.
				 * @param row The row.
				 */
				Cell(size_t id, size_t seedId, int col, int row);

				/**
				 * Construct a new cell with the grid coordinates and an explicit ID and value.
				 * @param id The ID.
				 * @param col The column.
				 * @param row The row.
				 * @param value The value.
				 */
				Cell(size_t id, size_t seedId, int col, int row, double value);

				/**
				 * Construct a new cell with the grid coordinates and an automatic ID and explicit value.
				 * @param col The column.
				 * @param row The row.
				 * @param value The value.
				 */
				Cell(size_t seedId, int col, int row, double value);

				/**
				 * Return the coordinate by axis index.
				 * @param idx The index
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
				 * @param seedId The seed ID.
				 * @param minc The minimum column index.
				 * @param minr The minimum row index.
				 * @param maxc The maximum column index.
				 * @param maxr The maximum row index.
				 * @param area The basin area.
				 */
				Basin(unsigned int seedId, int minc, int minr, int maxc, int maxr, int area);

				/**
				 * Make a list of all the cells that are on the edges of this basin and are
				 * also a local minimum.
				 *
				 * @param grd The flood fill grid.
				 * @param dem The original elevation raster.
				 * @param edgeCells A KDTree to insert the cells into.
				 * @return The number of cells found.
				 */
				int computeEdges(Grid& grd, Grid& dem, KDTree<Cell>& edgeCells);

				/**
				 * Return the basins seed ID.
				 *
				 * @return The seed ID.
				 */
				unsigned int seedId() const;

			};

			// A flood fill operator that fills pixels whose values are lower than
			// the given elevation.
			class LEFillOperator: public FillOperator<double, int> {
			private:
				Grid* m_src;
				Grid* m_dst;
				double m_elevation;
				double m_nodata;
				unsigned int m_id;
				int m_nofill;
				int m_srcBand;
				int m_dstBand;

			public:

				LEFillOperator(Grid* src, int srcBand, Grid* dst, int dstBand, double elevation, unsigned int id);

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
				size_t m_id;
				Cell m_c1;
				Cell m_c2;
				double m_elevation;
				double m_x;
				double m_y;
				double m_z;

			public:
				std::vector<double> path; ///<! The optimal spill path as triples of x, y, z.

				size_t id() const;

				/**
				 * Create a spill point. Copies the given cells.
				 */
				SpillPoint(const Cell& c1, const Cell& c2, double elevation, double x, double y, double z);

				const Cell& cell1() const;

				const Cell& cell2() const;

				double x() const;

				double y() const;

				double z() const;

				double elevation() const;

				/**
				 * Return the centroid of the line connecting the two cells.
				 */
				void centroid(int* col, int* row) const;

				~SpillPoint();

			};

			class SPDB : public geo::db::DB {
			public:
				 SPDB(const std::string &file, const std::string &layer, const std::string &driver,
					GeomType type, int srid = 0, bool replace = false);
				void addSpillPoint(const SpillPoint& sp, const GridProps& props);
			};

		} // util

		using namespace flood::util;

		class FloodCallbacks : public geo::util::Callbacks {
		public:
			void stepCallback(float status) const {}
			void overallCallback(float status) const {}
			void statusCallback(const std::string& msg) const {}
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
			MemRaster* m_dem;
			std::string m_input;
			std::string m_vdir;
			std::string m_rdir;
			std::string m_spill;
			std::string m_fseeds;
			std::string m_outfile;
			std::vector<Cell> m_seeds;
			std::vector<Basin> m_basinList;

			std::mutex m_mtxRast;
			std::mutex m_mtxDb;
			std::mutex m_mtxQueue;
			std::mutex m_mtxTree;

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
					double minBasinArea, double maxSpillDist);

			/**
			 * Copy the given Flood object.
			 */
			Flood(const Flood& conf);

			/**
			 * Check the inputs, throw an exception if any are invalid.
			 */
			void validateInputs();

			~Flood();

			std::mutex& mtxRast();

			std::mutex& mtxDb();

			std::mutex& mtxQueue();

			const std::string& input() const;

			const std::string& vdir() const;

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

			MemRaster* dem() const;

			/**
			 * Initialize the flood processor.
			 *
			 * @param mapped Set to true to use file-mapped memory.
			 */
			void init(bool mapped);

			/**
			 * Perform flood filling and identify basins.
			 * @param basinRaster An empty unique_ptr to contain a MemRaster.
			 * @param elevation The elevation to fill to.
			 * @param mapped Set to true to used mapped memory. Slow.
			 * @return An int representing the number of basins created.
			 */
			int fillBasins(std::unique_ptr<MemRaster>& basinRaster, double elevation, bool mapped);

			/**
			 * Vectorize the given raster.
			 * @param rfile the raster file.
			 * @param vfile the vector file.
			 */
			void saveBasinVector(const std::string& rfile, const std::string& vfile);

			/**
			 * Find the spill points between basins in the given raster.
			 *
			 * @param basinRaster The raster containing basins. Basins are contiguous regions with the same integer ID.
			 * @param elevation The flood elevation used to create the basins.
			 * @return True on success.
			 */
			bool findSpillPoints(MemRaster& basinRaster, std::vector<SpillPoint>& spillPoints, double elevation);

			/**
			 * Output the spill points to a stream, with comma delimiters.
			 * The fields are: ID1, x1, y1, ID2, x2, y2, midpoint x, midpoint y, distance
			 */
			void saveSpillPoints(unsigned int* id, std::vector<SpillPoint>& spillPoints, SPDB& db);

			/**
			 * Find the cells at the bottoms of depressions.
			 */
			void findMinima();

			/**
			 * Called by flood, used by threads.
			 */
			static void worker(Flood* config, SPDB* ofs, std::queue<double>* elevations, bool mapped);

			/**
			 * Start the flood process.
			 *
			 * @param numThreads The number of threads to run.
			 * @param memMapped Set to true to use a file-backed temporary raster.
			 */
			void flood(int numThreads, bool memMapped);

		};

	} // flood

} // geo

#endif /* INCLUDE_FLOOD_HPP_ */
