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

void usage() {
	std::cerr << "Usage: flood <options>\n"
			<< " -i     Input elevation raster.\n"
			<< " -s     CSV containing seed points. Row format, ID, X, Y. If not given, minima are used.\n"
			<< " -dbc   The database connection string for output vectors.\n"
			<< " -dbl   The basin database layer name (MultiPolygon).\n"
			<< " -dbi   The basin database ID field.\n"
			<< " -dbe   The basin database elevation field.\n"
			<< " -dbg   The spill database geometry (MultiPolygon) field.\n"
			<< " -dsl   The spill database layer name (LineString).\n"
			<< " -dsi   The spill database ID field.\n"
			<< " -dse   The spill database elevation field.\n"
			<< " -dsm   The spill database maximum elevation field.\n"
			<< " -dsb   The basin ID fields. Comma-separated list of two column names.\n"
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
			<< " -bl    Add a break line of the specified height. In the form <x0,y0,x1,y1,height>.\n"
			<< " -o     If given, each elevation will be computed even if raster file exists.\n"
			<< "        Otherwise, the elevation will be skipped.\n"
			<< " -a     The source raster band. Default 1.\n";
}

using namespace geo::flood;

int main(int argc, char **argv) {

	geo::loglevel(G_LOG_DEBUG);

	std::string input;
	std::string seeds;
	std::string rdir;
	std::string vdir;
	std::string spill;
	double start = 0.0;
	double end = 0.0;
	double step = 0.0;
	double maxSpillDist = 100.0;
	double minBasinArea = 100.0;
	int t = 1;
	int band = 1;
	bool overwrite = false;
	std::vector<geo::flood::BreakLine> breakLines;
	SpillOutput* spillOutput;
	BasinOutput* basinOutput;

	std::string bconn;
	std::string blayer;
	std::string bid;
	std::string belev;

	std::string sconn;
	std::string slayer;
	std::string sid;
	std::string selev;
	std::string smaxElev;
	std::vector<std::string> sbid;

	for (int i = 1; i < argc; ++i) {
		std::string a(argv[i]);
		if (a == "-i") {
			input.assign(argv[++i]);
		} else if (a == "-s") {
			seeds.assign(argv[++i]);
		} else if (a == "-v") {
			vdir = argv[++i];
		} else if (a == "-r") {
			rdir = argv[++i];
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
		} else if (a == "-o") {
			overwrite = true;
		} else if (a == "-a") {
			band = atoi(argv[++i]);
		} else if (a == "-bl") {
			std::vector<std::string> list;
			geo::util::split(std::back_inserter(list), argv[++i], ",");
			if(list.size() < 5)
				g_runerr("Need five elements in the breakline.");
			breakLines.emplace_back(std::stod(list[0]), std::stod(list[1]), std::stod(list[2]), std::stod(list[3]), std::stof(list[4]));
		} else if(a == "-dbc") {
			bconn = argv[++i];
		} else if(a == "-dbl") {
			blayer = argv[++i];
		} else if(a == "-dbi") {
			bid = argv[++i];
		} else if(a == "-dbe") {
			belev = argv[++i];
		} else if(a == "-dsc") {
			sconn = argv[++i];
		} else if(a == "-dsl") {
			slayer = argv[++i];
		} else if(a == "-dsi") {
			sid = argv[++i];
		} else if(a == "-dse") {
			selev = argv[++i];
		} else if(a == "-dsm") {
			smaxElev = argv[++i];
		} else if(a == "-dsb") {
			split(std::back_inserter(sbid), argv[++i], ",");
		}
	}

	if(!bconn.empty()) {
		BasinDBOutput* db = dynamic_cast<BasinDBOutput*>(basinOutput = new BasinDBOutput());
		db->conn = bconn;
		db->layer = blayer;
		db->idField = bid;
		db->elevationField = belev;
	} else if(!vdir.empty()){
		BasinFileOutput* db = dynamic_cast<BasinFileOutput*>(basinOutput = new BasinFileOutput());
		db->vdir = vdir;
	} else {
		basinOutput = new DummyBasinOutput();
	}
	basinOutput->rdir = rdir;

	if(!sconn.empty()) {
		SpillDBOutput* db = dynamic_cast<SpillDBOutput*>(spillOutput = new SpillDBOutput());
		db->conn = sconn;
		db->layer = slayer;
		db->idField = sid;
		db->elevationField = selev;
		db->maxElevationField = smaxElev;
		db->bid1Field = sbid[0];
		db->bid2Field = sbid[1];
	} else if(!spill.empty()){
		SpillFileOutput* db = dynamic_cast<SpillFileOutput*>(spillOutput = new SpillFileOutput());
		db->spillFile = spill;
	} else {
		spillOutput = new DummySpillOutput();
	}

	try {
		Flood config(input, band, overwrite, basinOutput, spillOutput, seeds, start, end,
			step, minBasinArea, maxSpillDist, breakLines);
		config.flood(t);
	} catch (const std::exception &e) {
		std::cerr << e.what() << "\n";
		usage();
		return 1;
	}

	return 0;
}


