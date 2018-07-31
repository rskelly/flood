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
#include "flood.hpp"

bool geo::flood::Flood::cancel = false;

void usage() {
	std::cerr << "Usage: flood <options>\n"
			<< " -i     Input elevation raster.\n"
			<< " -s     CSV containing seed points. Row format, ID, X, Y. If not given, minima are used.\n"
			<< " -v     Directory for basin vectors. If not given, they are not produced.\n"
			<< " -r     Directory for basin rasters. If not given, they are produced.\n"
			<< " -p     Spill point file. If not given, they are not produced. A Shapefile or CSV.\n"
			<< " -o     Output file. This is only available if the -h (height) flag is used.\n"
			<< " -start Starting elevation, in same units as input raster.\n"
			<< " -end   Ending elevation.\n"
			<< " -step  Step elevation.\n"
			<< " -h     Height. Performs only one flood operation at the specified height. If the -o flag \n"
			<< "        is used, the output is saved to that file."
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
	std::string outfile;
	double start = 0.0;
	double end = 0.0;
	double step = 0.0;
	double height = 0.0;
	bool hasHeight = false; // If true, start, end and step are disallowed.
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
			spill = argv[++i];
		} else if (a == "-o") {
			outfile = argv[++i];
		} else if (a == "-h") {
			height = atof(argv[++i]);
			hasHeight = true;
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

		if(hasHeight) {
			start = end = height;
			step = 1;
		} else {
			if(!outfile.empty())
				throw std::runtime_error("Output file is only for when the -h flag is specified.");
		}

		geo::flood::Flood config(input, vdir, rdir, spill, seeds, outfile, start, end,
				step, minBasinArea, maxSpillDist);
		config.flood(t, true);

	} catch (const std::exception &e) {
		std::cerr << e.what() << "\n";
		usage();
		return 1;
	}

	return 0;
}


