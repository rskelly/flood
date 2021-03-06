# Flood

`Flood` is an application that delinates basins on a digital terrain 
model (DTM) and attempts to locate the 'spill points' between them.

## Background

Typical watershed and channel delination techniques us a least-cost 
search algorithm to locate the path that a hypothetical water droplet 
will follow as it moves down slope. This algorithm "floods" basins, 
from user-defined seeds or surface minima, and locates the 
pixels on their edges which are within a specified distance of
each other. When such pairs of pixels are found, the least-cost path
(A*) between them is found and recorded as a *potential* spill path,
along with the flood elevation and the maximum height attained by
the path.

<img src="assets/pad_map.png" alt="Peace-Athabasca Delta" align="right" /> This algorithm was inspired by the unique hydrological conditions ([1](https://onlinelibrary.wiley.com/doi/abs/10.1002/hyp.6420), [2](https://link.springer.com/article/10.1007/s11273-005-1114-1))
that prevail in the [Peace-Athabasca Delta](https://en.wikipedia.org/wiki/Peace%E2%80%93Athabasca_Delta) of Northern Alberta, Canada:
the Delta is littered with "perched" basins that recharge neither
from precipitation nor groundwater, but from overland spillage due to ice jam
flooding from the Peace and Athabasca Rivers. 

The relief of the delta is extraordinarily flat. Excluding Shield outcrops, 
over the 5000km<sup>2</sup> expanse of the Delta, the elevation change might be 1-2m. For 
each basin, there is a definite least-cost path, according to the algorithm 
(there can't not be!) but this strategy isn't necessarily reflective of 
real-world processes. The `Flood` algorithm is designed to suggest the locations
and elevations of *likely* spill points.

This work was begun as part of an internship with [Dr. Daniel Peters](https://profils-profiles.science.gc.ca/en/profile/daniel-l-peters-phd-pgeo), with 
Environment and Climate Change Canada's [Water and Climate Impacts Research 
Centre](https://www.uvic.ca/research/centres/wcirc/) (W-CIRC) at the University of Victoria.

## Operation

`Flood` takes as input a DTM and an optional list of seed points representing locations
within basins. If seeds are not provided, the minima in the raster are used. Starting from the `start`
elevation and ending at the `end`, the program iteratively
fills pixels below the current elevation using the [flood fill](https://en.wikipedia.org/wiki/Flood_fill) algorithm. 
At each iteration, the current elevation is increased by the value of the `step`.
Regions of filled pixels comprise a "basin" and are given the ID of the seed point. A minimum basin
area filters out basins that are too small to be of interest (as they grow, basins are rapidly
subsumed by neighbouring basins so that the many initial seeds coalesce into a more manageable
number of basins). 

At each iteration, the edge pixels of each basin are placed in a [quadtree](https://en.wikipedia.org/wiki/Quadtree), 
one for each basin. Pairs of trees are searched for pixels within a configured range. For each pair,
the [A* least-cost search](https://en.wikipedia.org/wiki/A*_search_algorithm) is used to find the lowest-elevation path from one basin to the other. The end points,
the path geometry, the current elevation, and the maximum elevation traversed by the path are recorded.

The program also outputs a series of integer rasters, one for each elevation, containing the 
delinated basins. The pixel values are the IDs of the original seeds. Optionally, basins 
can be vectorized and output to files or a spatial database (PostGIS), though this step can be time-
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

## Eye Candy

[![Youtube Video](https://img.youtube.com/vi/UHseercT3Zg/0.jpg)](https://youtu.be/UHseercT3Zg)

Click image to view on Youtube.
