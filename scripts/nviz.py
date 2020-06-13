#!/usr/bin/env python3

'''
*** This program must be run in a GRASS GIS environment. ***

Produces 3D renderings of flood basins for the production of videos.

The flood basin rasters are produced by the flood program (https://github.com/rskelly/flood).

A flood basin raster contains regions of contiguous integral values delineating basins. The value is the ID of the basin.

The basins are rendered over a LiDAR-derived terrain model with an RGB image draped over it.

Viewpoints, targets, angles of view, etc. are configured for each video sequence.
'''

import os
import sys
import gdal
import math
import numpy as np
import threading as th
import traceback as tb
import psycopg2 as pg

# Set to true to overwrite existing outputs. Otherwise, 
# set to false and merely delete the images you want redone.
overwrite = False

# Database connection parameters for vectors.
dbname = 'ec'
user = 'rob'
password = 'river'
host = 'localhost'
layer = 'spill_4m'

# Rendermode. Coarse is useful for first-run checking.
nviz_mode = 'fine' # 'coarse'

# Sequence of render jobs. The params are:
# active: Set to true to run the job; false to ignore it.
# name: A name for the job.
# out_dir: The output directory for rendered imaged (*.ppm).
# out_pattern: The pattern for the output file names. Formatted with e -> elevation and n -> frame number.
# basin_dir: The location of basin rasters produced by the flood tool.
# basin_pattern: The pattern for the file names of the basin rasters. Same format parameters as out_pattern.
# dem_file: The elevation model raster.
# dem_name: The GRASS name for the DEM.
# rgb_file: An RGB orthophoto to drape over the DEM.
# rgb_name: A GRASS name for the orthophoto (colour bands will be composited).
# elev_range: The start and end elevations. Must correspond to the basin file names.
# elev_step: The elevation step.
# elev_mult: A multiplier for the elevation to produce the value used in basin file names.
# size: The size of the output images.
# position: The position of the viewer above the model in model coordinates at start and end. Will be graduated.
# height: The height of the viewer above the model.
# zexag: Z-exaggeration for elevations.
# focus: The position of focus within the scene at the start and end. Will be graduated.
# steps_per_elev: How many frames produced for each elevation step.
# perspective: Angle of view.
# color_table: Colour table for basins.
# basin_ids: The basin IDs of interest. Determines the loading of spill paths.
sequences = [
	{
		'active': False,
		'name': 'Horseshoe Slough',
		'out_dir': '/home/rob/Desktop/ec/videos/horseshoe_slough',
		'out_pattern': '{e:d}_{n}',
		'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
		'basin_pattern': '{e:d}.tif',
		'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
		'dem_name': 'pad_dem_4m',
		'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
		'rgb_name': 'pad_ortho',
		'elev_range': (209.5, 211.),
		'elev_step': 0.01,
		'elev_mult': 10000., 
		'size': (1920, 1080),
		#'position': ((465087.273, 6525263.298), (465290.378, 6525676.061)), 
		'position': ((.41, .455), (.42, .465)),
		'height': 500,
		'zexag': 2,
		'focus': ((466488.145, 6525389.335, 250), (466488.145, 6525889.335, 200)),
		'steps_per_elev': 3,
		'perspective': 20.,
		'color_table': '/home/rob/Desktop/ec/videos/color.txt',
		'basin_ids': (16, 23)
	},
	{
		'active': True,
		'name': 'PAD 58/Lake 540',
		'out_dir': '/home/rob/Desktop/ec/videos/pad58_lake540',
		'out_pattern': '{e:d}_{n}',
		'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
		'basin_pattern': '{e:d}.tif',
		'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
		'dem_name': 'pad_dem_4m',
		'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
		'rgb_name': 'pad_ortho',
		'elev_range': (209.5, 211.5),
		'elev_step': 0.01,
		'elev_mult': 10000., 
		'size': (1920, 1080),
		'position': ((.433, .5), (.437, .49)),
		'height': 500,
		'zexag': 2,
		'focus': ((468002.823, 6521240.493, 250), (468261.619, 6521286.356, 200)),
		'steps_per_elev': 3,
		'perspective': 30.,
		'color_table': '/home/rob/Desktop/ec/videos/color.txt',
		'basin_ids': (17, 18, 23)
	},
	{
		'active': False,
		'name': 'Lake 50',
		'out_dir': '/home/rob/Desktop/ec/videos/lake50',
		'out_pattern': '{e:d}_{n}',
		'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
		'basin_pattern': '{e:d}.tif',
		'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
		'dem_name': 'pad_dem_4m',
		'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
		'rgb_name': 'pad_ortho',
		'elev_range': (209.5, 211.5),
		'elev_step': 0.01,
		'elev_mult': 10000., 
		'size': (1920, 1080),
		'position': ((.417, .475), (.420, .48)),
		'height': 500,
		'zexag': 2,
		'focus': ((482157.983, 6521967.743, 350), (482655.920, 6522439.472, 300)),
		'steps_per_elev': 3,
		'perspective': 45.,
		'color_table': '/home/rob/Desktop/ec/videos/color.txt',
	},
	{
		'active': False,
		'name': 'PAD 37 - Rocher Pond',
		'out_dir': '/home/rob/Desktop/ec/videos/pad37',
		'out_pattern': '{e:d}_{n}',
		'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
		'basin_pattern': '{e:d}.tif',
		'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
		'dem_name': 'pad_dem_4m',
		'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
		'rgb_name': 'pad_ortho',
		'elev_range': (209.5, 211.5),
		'elev_step': 0.01,
		'elev_mult': 10000., 
		'size': (1920, 1080),
		'position': ((.417, .475), (.420, .48)),
		'height': 500,
		'zexag': 2,
		'focus': ((482957.303, 6522072.572, 350), (482393.848, 6521443.599, 300)),
		'steps_per_elev': 3,
		'perspective': 45.,
		'color_table': '/home/rob/Desktop/ec/videos/color.txt',
	},
	{
		'active': False,
		'name': 'PAD 5b - Mud Lake',
		'out_dir': '/home/rob/Desktop/ec/videos/pad5b',
		'out_pattern': '{e:d}_{n}',
		'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
		'basin_pattern': '{e:d}.tif',
		'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
		'dem_name': 'pad_dem_4m',
		'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
		'rgb_name': 'pad_ortho',
		'elev_range': (208.5, 210.2),
		'elev_step': 0.01,
		'elev_mult': 10000., 
		'size': (1920, 1080),
		'position': ((.417, .475), (.420, .48)),
		'height': 500,
		'zexag': 2,
		'focus': ((481044.178, 6517702.522, 300), (480185.892, 6517204.586, 300)),
		'steps_per_elev': 3,
		'perspective': 45.,
		'color_table': '/home/rob/Desktop/ec/videos/color.txt',
	},
	{
		'active': False,
		'name': 'Lake 582',
		'out_dir': '/home/rob/Desktop/ec/videos/lake582',
		'out_pattern': '{e:d}_{n}',
		'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
		'basin_pattern': '{e:d}.tif',
		'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
		'dem_name': 'pad_dem_4m',
		'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
		'rgb_name': 'pad_ortho',
		'elev_range': (208.5, 210.0),
		'elev_step': 0.01,
		'elev_mult': 10000., 
		'size': (1920, 1080),
		'position': ((.417, .475), (.420, .48)),
		'height': 500,
		'zexag': 2,
		'focus': ((479063.897, 6520596.779, 300), (478959.068, 6520039.876, 300)),
		'steps_per_elev': 3,
		'perspective': 45.,
		'color_table': '/home/rob/Desktop/ec/videos/color.txt',
	}		
]

class Job(th.Thread):
	'''
	A processing job in the grass environment.
	'''

	def __init__(self, tid, queue):
		'''
		Create a job with the given thread id, from configs in the queue.'
		'''
		super(Job, self).__init__()
		self.queue = queue;
		self.tid = tid

	def run(self):
		'''
		Run the job.
		'''
		while len(self.queue):
			# Extract configs from the queue.
			elev, frame, frames, basin, basin_elev, out_tpl, steps, seq = self.queue.pop(0)
			if overwrite or not has_file(elev, frame, **seq):
				try:
					# Load the basin rasters.
					load_basin(self.tid, elev, basin, basin_elev, **seq)
					# Render the scene.
					render_scene(self.tid, elev, frame, frames, basin, basin_elev, out_tpl, steps, **seq)
				except Exception as e: 
					print(tb.format_exc())
					break

def run():
	'''
	Entry point. Run all the jobs.
	'''

	# Iterate over active sequences.
	for seq in sequences:
		if not seq['active']:
			continue

		print('Processing {n}'.format(n = seq['name']))

		# Prepare output dir and load the basemaps.
		create_out_dir(**seq)
		load_basemap(**seq)

		# Prepare start, end and counters.
		elev, elev_end = list(map(float, seq['elev_range']))
		step = float(seq['elev_step'])
		steps = int(math.ceil((elev_end - elev) / step))
		spe = int(seq['steps_per_elev'])

		# Prepare the output file template.
		out_tpl = os.path.join(seq['out_dir'], seq['out_pattern'])

		# Prepare frame counter.
		frame = 1
		frames = steps * spe

		# Create the list of job configurations.
		jobs = []
		while elev <= elev_end:
			jobs.append((elev, frame, frames, 'basin_{id}', 'basin_elev_{id}', out_tpl, steps, seq))
			frame += spe
			elev += step

		# Run the jobs and wait for completion.
		# Using threads here because the work happens native code and the GIL 
		# helps us here by managing the queue's access.
		numthreads = 2
		threads = []
		for i in range(numthreads):
			threads.append(Job(i, jobs))
			threads[i].start()
		for i in range(numthreads):
			threads[i].join()


def has_file(elev, frame, elev_mult, steps_per_elev, out_dir, **seq):
	'''
	See if the sequence of files associated with the given elevation
	already exists. If any of the files is missing, all will be redone.
	'''
	e = int(round(elev * 100.)) * int(round(elev_mult / 100.))
	for i in range(frame, frame + steps_per_elev):
		f = os.path.join(out_dir, '{0}_{1}.ppm'.format(e, i))
		print(f)
		if not os.path.exists(f):
			return False
	return True

def get_raster_value(ds, band, x, y):
	'''
	Get the raster value at the given position.
	'''
	trans = ds.GetGeoTransform()
	b = ds.GetRasterBand(band)
	c = int((x - trans[0]) / trans[1])
	r = int((y - trans[3]) / trans[5])
	d = b.ReadAsArray(c, r, 1, 1)
	return d[0,0]

def get_position(x, y, zexag, height, raster):
	'''
	Calculate the model coordinate for nviz, which is 0-1 in x/y, from the 
	southwest corner of the image.

	TODO: This doesn't work. Not clear how it should.
	'''
	ds = gdal.Open(raster)
	trans = ds.GetGeoTransform()
	w = abs(ds.RasterXSize * trans[1])
	h = abs(ds.RasterYSize * trans[5])
	wx = abs(ds.RasterXSize * trans[1])
	wy = abs(ds.RasterYSize * trans[5])
	r = min(wx, wy)
	range = 5000
	range_offset = 2000
	# Center position.
	cx = trans[0] + (ds.RasterXSize / 2.) * trans[1]
	cy = trans[3] + (ds.RasterYSize / 2.) * trans[5]
	# Distance from center.
	dx = (x - cx) * (1 if trans[1] > 0 else -1)
	dy = (y - cy) * (1 if trans[5] > 0 else -1)
	z = get_raster_value(ds, 1, x, y)
	print(dx, dy, cx, cy, z)
	px = math.atan2((r - z - height) * zexag, dx) / math.pi
	py = math.atan2((r - z - height) * zexag, dy) / math.pi
	print(px, py)
	px = (px + range_offset) / range
	py = (py + range_offset) / range
	print(px, py)
	return (px, py)

def get_focus(x, y, z, raster):
	'''
	Calculate the focus for nviz, which is metres in x/y/z, from the 
	southwest corner of the image.
	'''
	ds = gdal.Open(raster)
	trans = ds.GetGeoTransform()
	bly = trans[3] + ds.RasterYSize * trans[5] if trans[5] < 0. else trans[3]
	blx = trans[0] + ds.RasterXSize * trans[1] if trans[1] < 0. else trans[0]
	#print('foc', x, y, blx, bly)
	return (x - blx, y - bly, z)


def create_out_dir(out_dir, **kwargs):
	'''
	Make the output directory. Make it and clear it out first if necessary.
	'''
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	for f in os.listdir(out_dir):
		try:
			os.path.unlink(os.path.join(out_dir, f))
		except: pass

def load_basemap(dem_file, dem_name, rgb_file, rgb_name, **kwargs):
	'''
	Load the dem and ortho, composite the ortho.
	'''

	# Load dem if not exists.
	try:
		os.system('r.in.gdal input={i} output={o}'.format(i = dem_file, o = dem_name))
	except: pass

	# Set computational region
	os.system('g.region raster={r}'.format(r = dem_name))

	# Load rgb if not exists
	try:
		os.system('r.in.gdal input={i} output={o}'.format(i = rgb_file, o = rgb_name))
		os.system('r.composite -c red={r}.red green={r}.green blue={r}.blue levels=32 output={r}'.format(r = rgb_name))
	except: pass


def load_basin(tid, elev, basin_name, basin_name_elev, basin_dir, basin_pattern, elev_mult, color_table, **kwargs):
	'''
	Load the basin raster and apply the colour table. 
	Create an elevation layer to position it vertically.
	'''

	# Create an elevation model for the basin raster.
	basin_name = basin_name.format(id = tid)
	basin_name_elev = basin_name_elev.format(id = tid)
	e = int(round(elev * 100.)) * int(round(elev_mult / 100.))

	# Create the file path.
	f = os.path.join(basin_dir, basin_pattern.format(e = e))
	print('Loading basin', f, '-->', basin_name)

	# Load the file.
	os.system('r.in.gdal --overwrite input={i} output={o}'.format(i = f, o = basin_name))

	# Set the color table.
	os.system('r.colors map={r} rules={c}'.format(r = basin_name, c = color_table))

	# We create a raster with the elevation lower than the desired elevation 
	# everywhere except where there are valid pixels. There, we raise the pixels
	# above the desired elevation to avoid occlusion (somewhat).
	os.system('r.mapcalc --overwrite expression="{o} = ({elev} - {e}) + ({b} > 0) * {e}"'.format(
		elev = elev, e = 2., b = basin_name, o = basin_name_elev)
	)

def load_conn(tid, name, basins, elevation):
	'''
	Load the connection vectors for the basins.
	'''

	# Format the basin IDs into a comma-delimited string.
	basins = ','.join(list(map(str, basins)))

	# Temp view for spill segments.
	tmp = 'spill_tmp_{}'.format(tid)

	# Connect to the spill DB.
	conn = pg.connect('dbname={dbn} user={dbu} password={dbp} host={dbh}'.format(
		dbn = dbname, dbu = user, dbp = password, dbh = host
	))
	cur = conn.cursor()

	# Create a table on the relevant spill connections.
	# Gotta do this because the query parser in grass seems to be broke.
	# Find the shortest paths between two basins at or below the given elevation.
	sql = '''
		create table {tmp} as \
		with t as ( \
			select gid, bid1, bid2, elevation, max_elevation, st_length(geom) as len, geom \
			from {lyr} \
			where bid1 in ({ids}) and bid2 in ({ids}) \
				and (elevation - {prec}) <= {elev} \
		) \
		select distinct on(bid1, bid2, elevation, len) bid1, bid2, elevation, geom \
		from t \
		group by bid1, bid2, elevation, geom, len \
		order by elevation, len \
	'''.format(
		lyr = layer,
		tmp = tmp,
		ids = basins, elev = elevation, prec = 0.001
	)
	cur.execute('drop table if exists {tmp}'.format(tmp = tmp))
	cur.execute(sql)
	cur.execute('alter table {tmp} add column _gid serial primary key'.format(tmp = tmp))
	conn.commit()
	conn.close()

	# Create the vector in grass for rendering.
	os.system('''v.external -o --overwrite \
			input="PG:dbname={dbn} user={dbu} password={dbp} host={dbh}" \
			layer={tmp} \
			output={name} \
		'''.format(
			dbn = dbname, dbu = user, dbp = password, dbh = host, 
			tmp = tmp, name = name, 
		)
	)

	return True

def render_scene(tid, elev, frame, frames, basin_name, basin_name_elev, out_tpl, steps, 
	dem_name, rgb_name, position, focus, size, zexag, elev_mult, steps_per_elev, 
	dem_file, height, perspective, basin_ids, **kwargs):
	'''
	Render the 3D image in m.nviz.image.
	'''

	# Format filenames.
	basin_name = basin_name.format(id = tid)
	basin_name_elev = basin_name_elev.format(id = tid)
	e = int(round(elev * 100.)) * int(round(elev_mult / 100.))

	# Start, end position and focus.
	pos0 = position[0]
	pos1 = position[1]
	foc0 = get_focus(focus[0][0], focus[0][1], focus[0][2], dem_file)
	foc1 = get_focus(focus[1][0], focus[1][1], focus[1][2], dem_file)

	# Load the connection lines.
	spill_name = 'spill_{}'.format(tid)
	load_conn(tid, spill_name, basin_ids, elev)

	# Iterate over the frames for this elevation step.
	for s in range(steps_per_elev):
		
		# The completion proportion, 0 - 1.
		frame_p = float(frame - 1) / frames

		# Current position, focus and perspective.
		pos = (
			pos0[0] + (pos1[0] - pos0[0]) * frame_p, 
			pos0[1] + (pos1[1] - pos0[1]) * frame_p
		)
		foc = (
			foc0[0] + (foc1[0] - foc0[0]) * frame_p, 
			foc0[1] + (foc1[1] - foc0[1]) * frame_p, 
			foc0[2] + (foc1[2] - foc0[2]) * frame_p
		)

		# Configure the output file.
		out_file = out_tpl.format(e = e, n = frame)

		# Advance frame.
		frame += 1

		print('Rendering', out_file)
		print('Pos', pos)
		print('Foc', foc)

		# Main part of the command.
		cmd = '''m.nviz.image -a --overwrite \
				elevation_map={dem},{basin} \
				color_map={dem_color},{basin_color} \
				resolution_fine=1 \
				resolution_coarse=9 \
				mode={mode} \
				position={position} \
				height={height} \
				zexag={zexag} \
				focus={focus} \
				output={output} \
				size={size} \
				perspective={persp}'''

		params = {
			'dem' : dem_name, 
			'dem_color' : rgb_name,
			'basin' : basin_name_elev, 
			'basin_color' :basin_name, 
			'position' : '{0[0]},{0[1]}'.format(pos),
			'focus' : '{0[0]},{0[1]},{0[2]}'.format(foc),
			'size' : '{0[0]},{0[1]}'.format(size),
			'height' : height,
			'zexag' : zexag,
			'persp' : perspective,
			'output' : out_file,
			'mode' : nviz_mode
		}

		# Add vline stuff if necessary.
		# TODO: May be conditional in the future.
		if True:
			cmd += '''\
				vline={vline} \
				vline_width={vwidth} \
				vline_height={vheight} \
				vline_color=red \
				vline_mode=flat'''

			params.update({
				'vline' : spill_name,
				'vwidth' : 2,
				'vheight' : 2
			})

		# Execute nviz.
		os.system(cmd.format(**params))
		
# Run me!
if __name__ == '__main__':
	run()
