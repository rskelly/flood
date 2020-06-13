# Flood

`Flood` is an application that delinates basins on a digital terrain 
model (DTM) and attempts to locate the 'spill points' between them.

Typical watershed and channel delination techniques us a least-cost 
search algorithm to locate the path that a hypothetical water droplet 
will follow as it moves down slope. This algorithm "floods" basins, 
from user-defined seeds or surface minima, and locates the 
pixels on their edges which are within a specified distance of
each other. When such pairs of pixels are found, the least-cost path
(A*) between them is found and recorded as a *potential* spill path,
along with the flood elevation and the maximum height attained by
the path.

This algorithm was inspired by the unique hydrological conditions
that prevail in the Peace-Athabasca Delta of Northern Alberta, Canada:
the Delta is littered with "perched" basins that recharge neither
from precipitation nor groundwater, but from overland spillage due to ice jam
flooding from the Peace and Athabasca Rivers. 

The relief of the delta is extraordinarily flat. Excluding Shield outcrops, 
over the 5000k expanse of the Delta, the elevation change might be 1-2m. For 
each basin, there is a definite least-cost path, according to the algorithm 
(there can't not be!) but this strategy isn't necessarily reflective of 
real-world processes. The `Flood` algorithm is designed to suggest the locations
and elevations of *likely* spill points.

This work was begun as part of an internship with Dr. Daniel Peters, with 
Environment and Climate Change Canada's Water and Climate Impacts Research 
Centre (W-CIRC) at the University of Victoria.

## Operation

`Flood` takes as input a DTM and an optional list of seed points representing locations
within basins. If seeds are not provided, the minima in the raster are used. A minimum basin
area filters out basins that are too small to be of interest (as they grow, basins are rapidly
subsumed by neighbouring basins so that the many initial seeds coalesce into a more manageable
number of basins). 

A start and end elevation, and an elevation step are also given. The program
will "flood" each elevation, step-by-step, from the start to the end (inclusive). 

The default output is a series of integer rasters, one for each elevation, containing the 
delinated basins. The pixel values are the IDs of the original seeds. Optionally, basins 
can be vectorized and output to files or a spatial database, though this step can be time-
consuming.

Spill points and path geometries can be output to a CSV file (with line geometries represented
as WKT linestrings) or to a spatial database.
 
## Installation

`Flood` is designed to run on Linux systems, using the usual `cmake` process:

1) `$ git clone https://github.com/rskelly/flood`
2) `$ cd flood && mkdir build && cd build`
3) `$ cmake ..`
4) `$ make && sudo make install`

Run `flood` to see a usage message.
Flooding and connectivity app.

## Eye Candy

*Soon*