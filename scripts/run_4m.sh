#!/bin/bash

# This script runs the flood program on a DEM.

dem=/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif
seeds=v_poi.csv
rdir=4m/r
spill=4m/spill.csv
estart=209
eend=213
estep=0.01
minarea=10000
maxdist=40

flood \
	-i $dem \
	-s $seeds \
	-r $rdir \
	-p $spill \
	-dsc "PG:dbname=ec user=rob password=river" \
	-dsl spill_4m \
	-dsi sid \
	-dse elevation \
	-dsm max_elevation \
    -dsb bid1,bid2 \
	-start $estart \
	-end $eend \
	-step $estep \
	-b $minarea \
	-d $maxdist

# Basin output.
#	-dbc "PG:dbname=ec user=rob password=river" \
#	-dbl basins_4m \
#	-dbi bid \
#	-dbe elevation \

# Breaklines (Ice Jams).
#	-bl 466415.3,6529654.5,467098,6529732,211.5 \
#	-bl 477151.9,6507498.0,477072.6,6507575.3,212.25 \
