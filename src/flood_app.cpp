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

#include "geo.hpp"
#include "flood.hpp"

bool geo::flood::Flood::cancel = false;

void usage() {
	std::cerr << "Usage: flood <options>\n"
			<< " -i     Input elevation raster.\n"
			<< " -s     CSV containing seed points. Row format, ID, X, Y. If not given, minima are used.\n"
			<< " -dbc   The database connection string for output vectors.\n"
			<< " -dbl   The database layer name.\n"
			<< " -dbi   The database ID field.\n"
			<< " -dbe   The database elevation field.\n"
			<< " -v     Directory for basin vectors. If not given, they are not produced.\n"
			<< " -r     Directory for basin rasters. If not given, they are produced.\n"
			<< " -p     Spill point file. If not given, they are not produced. A Shapefile or CSV.\n"
			<< " -start Starting elevation, in same units as input raster.\n"
			<< " -end   Ending elevation.\n"
			<< " -step  Step elevation.\n"
			<< " -t     Number of threads to use. Default 1.\n"
			<< " -b     Minimum basin area.\n"
			<< " -d     Maximum spill distance.\n"
			<< " -m     Use file-backed memory for intermediate rasters.\n"
			<< " -bl    Add a break line of the specified height. In the form <x0,y0,x1,y1,height>.";
}

int main(int argc, char **argv) {

	geo::loglevel(G_LOG_DEBUG);

	std::string input;
	std::string seeds;
	std::string rdir;
	std::string spill;
	double start = 0.0;
	double end = 0.0;
	double step = 0.0;
	double maxSpillDist = 100.0;
	double minBasinArea = 100.0;
	int t = 1;
	std::vector<geo::flood::BreakLine> breakLines;
	geo::flood::DBOutput db;
	geo::flood::FileOutput file;

	for (int i = 1; i < argc; ++i) {
		std::string a(argv[i]);
		if (a == "-i") {
			input.assign(argv[++i]);
		} else if (a == "-s") {
			seeds.assign(argv[++i]);
		} else if (a == "-v") {
			file.vdir = argv[++i];
		} else if (a == "-r") {
			rdir = argv[++i];
			db.rdir = rdir;
			file.rdir = rdir;
		} else if (a == "-p") {
			spill = argv[++i];
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
		} else if (a == "-bl") {
			std::vector<std::string> list;
			geo::util::split(std::back_inserter(list), argv[++i], ",");
			if(list.size() < 5)
				g_runerr("Need five elements in the breakline.");
			breakLines.emplace_back(std::stod(list[0]), std::stod(list[1]), std::stod(list[2]), std::stod(list[3]), std::stof(list[4]));
		} else if(a == "-dbc") {
			db.conn = argv[++i];
		} else if(a == "-dbl") {
			db.layer = argv[++i];
		} else if(a == "-dbi") {
			db.id = argv[++i];
		} else if(a == "-dbe") {
			db.elevation = argv[++i];
		}
	}

	try {

		if(db.valid()) {
			geo::flood::Flood config(input, 1, db, spill, seeds, start, end,
				step, minBasinArea, maxSpillDist, breakLines);
			config.flood(t);
		} else if(file.valid()) {
			geo::flood::Flood config(input, 1, file, spill, seeds, start, end,
				step, minBasinArea, maxSpillDist, breakLines);
			config.flood(t);
		} else {
			g_runerr("No valid output configuration.");
		}

	} catch (const std::exception &e) {
		std::cerr << e.what() << "\n";
		usage();
		return 1;
	}

	return 0;
}


