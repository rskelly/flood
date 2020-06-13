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
			 * \brief Represents a grid cell.
			 *
			 * Contains the grid coordinates, an ID and a value, such as elevation.
			 *
			 * The cell's priority dictates the order in which basins are merged. A basin
			 * with a lower priority value will impart its ID to a higher-valued basins.
			 */
			class Cell {
			private:
				int m_col;			///<! The cell column.
				int m_row;			///<! The cell row.
				int m_cellId;		///<! The cell ID.
				float m_value;		///<! The cell value.
				int m_seedId;		///<! The ID of the seed that originated this cell.
				int m_priority;		///<! The cell's priority.

			public:

				/**
				 * \brief Construct an empty cell with an automatic ID.
				 */
				Cell();

				/**
				 * \brief Construct a new cell with the grid coordinates and an explicit ID and value.
				 *
				 * \param id The ID.
				 * \param col The column.
				 * \param row The row.
				 * \param value The value.
				 * \param priority The basin's priority.
				 */
				Cell(int id, int seedId, int col, int row, float value, int priority);

				/**
				 * \brief Construct a new cell with the grid coordinates and an automatic ID and explicit value.
				 *
				 * \param col The column.
				 * \param row The row.
				 * \param value The value.
				 */
				Cell(int seedId, int col, int row, float value, int proiority = geo::maxvalue<int>());

				/**
				 * \brief Return the coordinate by axis index.
				 *
				 * The value is a double to accommodate the tree class template.
				 *
				 * \param idx The index
				 * \return The coordinate value.
				 */
				double operator[](int idx) const;

				/**
				 * \brief Return true if this Cell is strictly less than the other Cell for ordering.
				 *
				 * \param other Another Cell.
				 * \return True if this Cell is strictly less than the other Cell for ordering.
				 */
				bool operator<(const Cell& other) const;

				/**
				 * \brief Return the x coordinate as a double.
				 *
				 * This is really the column value.
				 *
				 * \return The x coordinate as a double.
				 */
				double x() const;

				/**
				 * \brief Return the y coordinate as a double.
				 *
				 * This is really the row value.
				 *
				 * \return The y coordinate as a double.
				 */
				double y() const;

				/**
				 * \brief Return the cell's row index in the raster.
				 *
				 * \return The cell's row.
				 */
				int row() const;

				/**
				 * \brief Return the cell's column index in the raster.
				 *
				 * \return The cell's column.
				 */
				int col() const;

				/**
				 * \brief Return the cell's ID.
				 *
				 * \return The cell's ID.
				 */
				int cellId() const;

				/**
				 * \brief Return the ID of the seed that originated this cell.
				 *
				 * \return The ID of the seed that originated this cell.
				 */
				int seedId() const;

				/**
				 * \brief Return the cell's value.
				 *
				 * \return The cell's value.
				 */
				float value() const;

				/**
				 * \brief Return the Cartesian distance between cells in map units, given the resolution of the original raster.
				 *
				 * \param other The other Cell.
				 * \param resx The x-resolution of the raster.
				 * \param resy The y-resolution of the raster.
				 * \return The Cartesian distance between this Cell and another.
				 */
				float distance(const Cell& other, float resx, float resy) const;

			};

			/**
			 * \brief Represents a basin, with bounds and area and the ID of the seed that spawned it.
			 */
			class Basin {
			private:
				int m_id;		///<! The ID of this basin, derived from the ID of the seed that originated it.
				int m_minc;		///<! The minimum raster column of the basin's bounding box.
				int m_minr;		///<! The minimum raster row of the basin's bounding box.
				int m_maxc;		///<! The maximum raster column of the basin's bounding box.
				int m_maxr;		///<! The maximum raster row of the basin's bounding box.
				int m_area;		///<! The the pixel area of the basin (not of its bounding box).

			public:

				/**
				 * \brief Construct a basin with the given seed ID, extent and area.
				 *
				 * \param seedId The seed ID.
				 * \param minc The minimum column index.
				 * \param minr The minimum row index.
				 * \param maxc The maximum column index.
				 * \param maxr The maximum row index.
				 * \param area The basin area.
				 */
				Basin(int seedId, int minc, int minr, int maxc, int maxr, int area);

				/**
				 * \brief Make a list of all the cells that are on the edges of this basin and are also a local minimum.
				 *
				 * \param grd The flood fill grid.
				 * \param edgeCells An mqtree to insert the cells into.
				 * \return The number of cells found.
				 */
				int computeEdges(Band<int>& grd, mqtree<Cell>& edgeCells);

				/**
				 * \brief Return the basin's seed ID.
				 *
				 * \return The seed ID.
				 */
				int seedId() const;

			};

			/**
			 * \brief A flood fill operator that fills pixels whose values are lower than the given elevation.
			 */
			class LEFillOperator: public FillOperator<float, int> {
			private:
				Band<float>* m_src;		///<! The source raster (DEM).
				Band<int>* m_dst;		///<! The destination raster (mask).
				float m_elevation;		///<! The fill elevation.
				float m_nodata;			///<! The source nodata value.
				int m_id;				///<! The ID to fill a region with.
				int m_target;			///<! The target value to fill with the ID (zero).

			public:

				/**
				 * \brief Construct the operator with the given rasters, target elevation and fill ID.
				 *
				 * \param src The source raster (DEM).
				 * \param dst The destination raster (mask).
				 * \param elevation The target elevation.
				 * \param id The ID to fill a region with.
				 */
				LEFillOperator(Band<float>* src, Band<int>* dst, float elevation, int id);

				/**
				 * \brief Return true if the cell should be filled.
				 *
				 * \param col The column to check.
				 * \param row The row to check.
				 * \return True if the cell should be filled.
				 */
				bool shouldFill(int col, int row) const;

				/**
				 * \brief Fill the region starting at the given column and row.
				 *
				 * \param col The column to fill from.
				 * \param row The row to fill from.
				 */
				void fill(int col, int row) const;

				/**
				 * \brief Set the fill value.
				 *
				 * \param fill The fill value.
				 */
				void setFill(int fill);

				/**
				 * \brief Set the target value in the destination.
				 *
				 * Pixels without this value are ignored.
				 *
				 * \param target The target value in the destination.
				 */
				void setTarget(int target);

				/**
				 * \brief Return the properties of the source raster.
				 *
				 * \return The properties of the source raster.
				 */
				const GridProps& srcProps() const;

				/**
				 * \brief Return the properties of the destination raster.
				 *
				 * \return The properties of the destination raster.
				 */
            	const GridProps& dstProps() const;

            	~LEFillOperator();

			};

			/**
			 * \brief Represents a spill point between two cells.
			 */
			class SpillPoint {
			private:
				Cell m_c1;			///<! A Cell representing one side of the spill.
				Cell m_c2;			///<! A Cell representing one side of the spill.
				float m_elevation;	///<! The spill elevation.
			public:

				/**
				 * \brief Create a spill point. Copies the given cells.
				 *
				 * \param c1 A Cell.
				 * \param c2 A Cell.
				 * \param elevation The spill elevation.
				 */
				SpillPoint(const Cell& c1, const Cell& c2, float elevation);

				/**
				 * \brief Return on of the Cells.
				 *
				 * \return A Cell.
				 */
				const Cell& cell1() const;

				/**
				 * \brief Return on of the Cells.
				 *
				 * \return A Cell.
				 */
				const Cell& cell2() const;

				/**
				 * \brief Return the spill elevation.
				 *
				 * \return The spill elevation.
				 */
				float elevation() const;

				/**
				 * \brief Return the centroid of the line connecting the two cells.
				 *
				 * \param[out] col The centroid column.
				 * \param[out] row The centroid row.
				 */
				void centroid(int& col, int& row) const;

				~SpillPoint();

			};

			/**
			 * \brief A heuristic for use in the A* search algorithm for spill paths.
			 *
			 * The heuristic calculates the distance and modifies it using the slope. A horizontal
			 * has the normal 2d length. The distance is multiplied by the sine of the angle
			 * of the slope w/r/t 1 map unit.
			 */
			class heuristic {
			public:
				Band<float>* dem;	///<! The raster to search on.

				/**
				 * \brief Construct the heuristic on the given raster.
				 *
				 * \param dem A raster to search.
				 */
				heuristic(Band<float>* dem);

				/**
				 * \brief The search operator.
				 *
				 * \param cc The start column.
				 * \param cr The start row.
				 * \param gc The goal column.
				 * \param gr The goal row.
				 * \return The cost associated with the path from cc/cr to gc/gr.
				 */
				float operator()(int cc, int cr, int gc, int gr);

			};

			/**
			 * \brief An abstract class for saving basin rasters and geometries.
			 *
			 * Subclasses will implement ways to save vectors. The base class
			 * provides the ability to save rasters.
			 */
			class BasinOutput {
			public:
				int id;				///<! The current record ID.
				std::string rdir;	///<! The raster output directory.

				/**
				 * \brief Return true if the object is configured correctly.
				 *
				 * \return True if the object is configured correctly.
				 */
				virtual bool valid() const;

				/**
				 * \brief Save vector versions of the basin rasters.
				 *
				 * \param The grid containing basins.
				 * \param elevation The elevation from which the basins are derived.
				 */
				virtual void saveVectors(Band<int>& basins, float elevation) = 0;

				/**
				 * \brief Prepare the object for writing output.
				 */
				void prepare();

				/**
				 * \brief Return the path of the file for saving basins corresponding to the given elevation.
				 *
				 * \param elevation The flood elevation.
				 * \return The path of the file for saving basins corresponding to the given elevation.
				 */
				std::string rasterFile(float elevation) const;

				virtual ~BasinOutput() {}
			};

			/**
			 * \brief A dummy implementation of the BasinOutput class.
			 *
			 * Saving vectors is a no-op. Rasters only.
			 */
			class DummyBasinOutput : public BasinOutput {
			public:

				bool valid() const;

				void saveVectors(Band<int>&, float);

			};

			/**
			 * \brief An abstract class for saving spill paths between basins.
			 *
			 * Subclasses will implement ways to save vectors.
			 */
			class SpillOutput {
			public:
				int id;					///<! The current record ID.
				std::string projection;	///<! The projection of the dataset.

				/**
				 * \brief Return true if the object is configured correctly.
				 *
				 * \return True if the object is configured correctly.
				 */
				virtual bool valid() const = 0;

				/**
				 * \brief Prepare the object for writing output.
				 */
				virtual void prepare() = 0;

				/**
				 * \brief Save the spill points/paths.
				 *
				 * \param dem The raster from which basins are derived.
				 * \param spillPoints The list of spill points.
				 */
				virtual void saveSpillPoints(Band<float>& dem, const std::vector<SpillPoint>& spillPoints) = 0;

				virtual ~SpillOutput() {}
			};

			/**
			 * \brief A dummy class for spill output.
			 *
			 * Saving spill points is a no-op.
			 */
			class DummySpillOutput : public SpillOutput {
			public:

				bool valid() const;

				void prepare();

				void saveSpillPoints(Band<float>& dem, const std::vector<SpillPoint>& spillPoints);

			};

			/**
			 * \brief Vectorizes and saves basin outlines to a relational database.
			 */
			class BasinDBOutput : public BasinOutput {
			public:
				std::string conn;			///<! The DB connection string. (e.g., "PG:dbname=... etc.")
				std::string layer;			///<! The DB layer (table) name.
				std::string idField;		///<! The primary key (integer) field.
				std::string elevationField;	///<! The elevation field.

				bool valid() const;

				void prepare();

				void saveVectors(Band<int>& rast, float elev);

			};

			/**
			 * \brief Saves spill points to a relational database.
			 *
			 * Spill paths are saved as linestring geometries.
			 */
			class SpillDBOutput : public SpillOutput {
			private:
				GDALDataset* ds;				///<! The dataset for writing.
				OGRLayer* lyr;					///<! The layer for writing.
				OGRGeometryFactory gf;			///<! Factory for geometries.
				OGRSpatialReference sr;			///<! Spatial reference.
			public:
				std::string conn;				///<! The database connection string. (e.g., "PG:dbname=... etc.")
				std::string layer;				///<! The database layer (table).
				std::string idField;			///<! The primary key (integer) field.
				std::string bid1Field;			///<! The ID of basin1.
				std::string bid2Field;			///<! The ID of basin2.
				std::string elevationField;		///<! The elevation field.
				std::string maxElevationField;	///<! The field for maximum elevation crossed by the path.

				bool valid() const;

				void prepare();

				void saveSpillPoints(Band<float>& dem, const std::vector<SpillPoint>& spillPoints);

				~SpillDBOutput();

			};

			/**
			 * \brief Vectorizes and saves basin rasters to files in a folder.
			 */
			class BasinFileOutput : public BasinOutput {
			public:
				std::string vdir;	///<! The vector output folder.

				bool valid() const;

				void prepare();

				void saveVectors(Band<int>& rast, float elevation);

			};

			/**
			 * \brief Saves spill points to a CSV file.
			 *
			 * Paths are serialized as WKT linestrings.
			 */
			class SpillFileOutput : public SpillOutput {
			public:
				std::string spillFile;		///!< The output file.
				std::ofstream out;			///!< The output stream.

				bool valid() const;

				void prepare();

				void saveSpillPoints(Band<float>& dem, const std::vector<SpillPoint>& spillPoints);
			};

			/**
			 * \brief Represents a break line.
			 *
			 * Use to create, e.g., an artificial blockage, such as
			 * an ice jam.
			 */
			class BreakLine {
			public:
				float x0;		///<! Start x-coordinate.
				float y0;		///<! Start y-coordinate.
				float x1;		///<! End x-coordinate.
				float y1;		///<! End y-coordinate.
				float value;	///<! The "height" of the breakline.

				/**
				 * \brief Construct a breakline.
				 *
				 * \param Start x-coordinate.
				 * \param Start y-coordinate.
				 * \param End x-coordinate.
				 * \param End y-coordinate.
				 * \param The "height" of the breakline.
				 */
				BreakLine(float x0, float y0, float x1, float y1, float value);

			};

		} // util


		using namespace flood::util;


		class Flood {
		private:
			float m_start;								///<! The start elevation.
			float m_end;								///<! The end elevation (inclusive).
			float m_step;								///<! The elevation step.
			float m_minBasinArea;						///<! Ignore basins smaller than this.
			float m_maxSpillDist;						///<! Ignore spill paths longer than this (as the crow flies).
			int m_t; 									///<! Number of threads
			int m_band;									///<! The source band in the DEM.
			Band<float> m_dem;							///<! The DEM.
			BasinOutput* m_basinOutput;					///<! The basin output object.
			SpillOutput* m_spillOutput;					///<! The spill point output object.
			std::string m_input;						///<! The filename of the DEM.
			std::string m_fseeds;						///<! The input handle for the seed file.
			std::vector<Cell> m_seeds;					///<! List of seed cells.
			std::vector<Basin> m_basinList;				///<! List of basins.
			std::vector<SpillPoint> m_spillPoints;		///<! List of spill points.
			std::vector<BreakLine> m_breakLines;		///<! List of breaklines.
			std::mutex m_qmtx; 							///<! Mutex to protect the job queue.

		public:

			bool cancel;	///<! If true, current operations are cancelled.

			/**
			 * \brief Create a Flood object.
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
			Flood(const std::string& input, int band,
					BasinOutput* basinOutput, SpillOutput* spillOutput,
					const std::string& seeds,
					float start, float end, float step,
					float minBasinArea, float maxSpillDist,
					const std::vector<BreakLine>& breakLines);

			/**
			 * \brief Return the queue mutex.
			 *
			 * \return The queue mutex.
			 */
			std::mutex& qmtx();

			/**
			 * \brief Check the inputs, throw an exception if any are invalid.
			 */
			void validateInputs();

			~Flood();

			/**
			 * \brief Return the input DEM file.
			 */
			const std::string& input() const;

			/**
			 * \brief Return a reference to the BasinOutput object.
			 *
			 * \return A reference to the BasinOutput object.
			 */
			BasinOutput& basinOutput();

			/**
			 * \brief Return a reference to the SpillOutput object.
			 *
			 * \return A reference to the SpillOutput object.
			 */
			SpillOutput& spillOutput();

			/**
			 * \brief Return a reference to the file containing seeds.
			 *
			 * \return A reference to the file containing seeds.
			 */
			const std::string& seedsFile() const;

			/**
			 * \brief Return a reference to the spill points list.
			 *
			 * \return A reference to the spill points list.
			 */
			const std::vector<SpillPoint>& spillPoints() const;

			/**
			 * \brief Return the start elevation.
			 *
			 * \return The start elevation.
			 */
			float start() const;

			/**
			 * \brief Return the end elevation (inclusive).
			 *
			 * \return The end elevation.
			 */
			float end() const;

			/**
			 * \brief Return the elevation step.
			 *
			 * \return The elevation step.
			 */
			float step() const;

			/**
			 * \brief Return the minimum basin area.
			 *
			 * \return The minimum basin area.
			 */
			float minBasinArea() const;

			/**
			 * \brief Return the maximum spill distance.
			 *
			 * \return The maximum spill distnance.
			 */
			float maxSpillDist() const;

			/**
			 * \brief Load the seeds from the seeds file.
			 *
			 * If header is given, the first row is read and discarded.
			 *
			 * \param header Set to true if the first line of the file contains column labels.
			 */
			void loadSeeds(bool header = false);

			/**
			 * \brief Return a reference to the list of seeds.
			 *
			 * \return A reference to the list of seeds.
			 */
			const std::vector<Cell>& seeds() const;

			/**
			 * \brief Return a reference to the input DEM.
			 *
			 * \return A reference to the input DEM.
			 */
			Band<float>& dem();

			/**
			 * \brief Initialize the flood processor.
			 */
			void init();

			/**
			 * \brief Perform flood filling and identify basins.
			 *
			 * \param basinRaster An empty unique_ptr to contain a Band<float>.
			 * \param elevation The elevation to fill to.
			 * \param mapped Set to true to used mapped memory. Slow.
			 * \return An int representing the number of basins created.
			 */
			int fillBasins(Band<int>& basinRaster, float elevation);

			/**
			 * \brief Find the spill points between basins represented in the raster.
			 *
			 * \param basinRaster The raster containing delineated basins.
			 * \param elevation The elevation for which basins are derived.
			 * \return True if spill points are found.
			 */
			bool findSpillPoints(Band<int>& basinRaster, float elevation);

			/**
			 * \brief Find the cells at the bottoms of depressions.
			 *
			 * This is used when no seed file is provided.
			 */
			void findMinima();

			/**
			 * \brief Start the flood process.
			 *
			 * \param numThreads The number of threads to run.
			 */
			void flood(int numThreads);

		};

	} // flood

} // geo

#endif /* INCLUDE_FLOOD_HPP_ */
