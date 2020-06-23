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

# Debug mode. Makes only the first and last images
# in coarse mode and always overwrites.
debug = False

# Set to true to overwrite existing outputs. Otherwise,
# set to false and merely delete the images you want redone.
overwrite = True

# Database connection parameters for vectors.
dbname = 'ec'
user = 'rob'
password = 'river'
host = 'localhost'
layer = 'spill_4m'

# Rendermode. Coarse is useful for first-run checking.
nviz_mode = 'fine' if not debug else 'coarse'

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
# region: The computational region.
# res: The grid resolution.
# steps_per_elev: How many frames produced for each elevation step.
# perspective: Angle of view.
# color_table: Colour table for basins.
# basin_ids: The basin IDs of interest. Determines the loading of spill paths.
sequences = [
        {
                'active': True,
                'name': 'Horseshoe Slough',
                'out_dir': '/home/rob/Desktop/ec/videos/horseshoe_slough',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (209.5, 211.),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((464641.750,6524696.567), (465237.964,6525004.502)),
                'height': 500,
                'zexag': 2,
                'focus': ((466261.683,6525412.351, 200.), (466261.683,6525412.351, 200.)),
                'region': (464632., 6530152., 474140., 6520368.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 40.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (16, 23)
        },
        {
                'active': True,
                'name': 'PAD 58/Lake 540',
                'out_dir': '/home/rob/Desktop/ec/videos/pad58_lake540',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (209.5, 212),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((467267.384,6519997.290), (466749.792,6519387.972)),
                'height': 500,
                'zexag': 2,
                'focus': ((468256., 6521197., 200.), (468656., 6521997., 200.)),
                'region': (464564., 6526792., 474072., 6517004.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 30.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (17, 18, 23)
        },
        {
                'active': True,
                'name': 'Lake 50',
                'out_dir': '/home/rob/Desktop/ec/videos/lake50',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (208.5, 211.5),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((482326.692,6522888.271), (481088.402,6522431.283)),
                'height': 500,
                'zexag': 2,
                'focus': ((482025.309,6521859.638, 200.), (482356.175,6521360.064, 200.)),
                'region': (476936., 6526864., 486444., 6517076.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 65.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (7, 21, 23),
        },
        {
                'active': True,
                'name': 'PAD 37 - Rocher Pond',
                'out_dir': '/home/rob/Desktop/ec/videos/pad37',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (209.0, 210.3),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((484745.944,6521523.859), (484647.667,6521530.410)),
                'height': 500,
                'zexag': 2,
                'focus': ((483671.449,6521432.133, 200), (483671.449,6521432.133, 200)),
                'region': (476936., 6526864., 486444., 6517076.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 45.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (7, 21, 23),
        },
        {
                'active': True,
                'name': 'PAD 5b - Mud Lake',
                'out_dir': '/home/rob/Desktop/ec/videos/pad5b',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (208.5, 210.5),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((478393.975,6516564.147), (478459.493,6516485.526)),
                'height': 750,
                'zexag': 2,
                'focus': ((480359.515,6517193.120, 200), (480359.515,6517193.120, 200)),
                'region': (474208., 6527604., 486584., 6512280.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 30.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (11, 14, 5, 15, 23),
        },
        {
                'active': True,
                'name': 'Lake 582',
                'out_dir': '/home/rob/Desktop/ec/videos/lake582',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (208.5, 210.0),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((478810.014,6520457.553), (478652.771,6520791.695)),
                'height': 500,
                'zexag': 2,
                'focus': ((479497.953,6519022.710, 200), (479347.262,6519350.300, 200)),
                'region': (471424., 6525612., 486348., 6511684.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 30.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (14, 22, 23),
        },
        {
                'active': True,
                'name': 'Edwards Lake',
                'out_dir': '/home/rob/Desktop/ec/videos/edwards_lk',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (209., 210.5),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((479897.613,6509537.344), (478705.186,6509275.272)),
                'height': 500,
                'zexag': 2,
                'focus': ((481470.045,6509170.443, 200), (481024.522,6509668.379, 200)),
                'region': (476596., 6513836., 490088., 6499824.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 45.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (9, 8, 22, 23),
        },
        {
                'active': True,
                'name': 'Flett Lake',
                'out_dir': '/home/rob/Desktop/ec/videos/flett_lk',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (209., 210.5),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((482053.155,6510726.495), (482990.062,6511499.607)),
                'height': 500,
                'zexag': 2,
                'focus': ((481679.702,6510333.387, 200), (481679.702,6510333.387, 200)),
                'region': (465612., 6514632., 483444., 6498432.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 30.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (9, 8, 22, 23),
        },
        {
                'active': True,
                'name': 'Lake 565',
                'out_dir': '/home/rob/Desktop/ec/videos/lake565',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (209.5, 210.5),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((479301.399,6515633.792), (478927.947,6513609.286)),
                'height': 500,
                'zexag': 2,
                'focus': ((480664.173,6514873.783, 200), (480408.653,6514729.644, 200)),
                'region': (473472., 6519164., 493212., 6499984.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 20.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (12, 11, 22, 23),
        },
        {
                'active': True,
                'name': 'Lake 577',
                'out_dir': '/home/rob/Desktop/ec/videos/lake577',
                'out_pattern': '{e}_{n}',
                'basin_dir': '/home/rob/Desktop/ec/videos/4m/r',
                'basin_pattern': '{e}.tif',
                'dem_file': '/home/rob/Desktop/ec/final_p10_4m_clipped_fill_min.tif',
                'dem_name': 'pad_dem_4m',
                'rgb_file': '/home/rob/Desktop/ec/videos/ortho.tif',
                'rgb_name': 'pad_ortho',
                'elev_range': (209.5, 210.5),
                'elev_step': 0.001,
                'elev_mult': 10000.,
                'size': (1920, 1080),
                'position': ((477833.796,6522251.108), (477689.657,6521759.723)),
                'height': 500,
                'zexag': 2,
                'focus': ((476988.614,6522631.113, 200), (477624.139,6521438.685, 200.)),
                'region': (471720., 6527580., 483780., 6514684.),
                'res': 4.,
                'steps_per_elev': 3,
                'perspective': 45.,
                'color_table': '/home/rob/Desktop/ec/videos/color.txt',
                'basin_ids': (1, 13, 15, 23),
        }
]

class Job(th.Thread):
    '''
    A processing job in the grass environment.
    '''

    def __init__(self, tid, queue):
        '''
        Create a job with the given thread id, from configs in the queu
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
            if overwrite or debug or not has_file(elev, frame, **seq):
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

        # Keep the first and last for debugging.
        if debug:
            jobs = jobs[:spe] + jobs[-spe:]

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

def e_to_filename(elev, elev_mult, elev_step):
    '''
    Convert the elevation to a filename-ready string. Need the multiplier
    and step value to fix the precision.
    '''
    es = int(1. / elev_step)
    return str(int(round(elev * es) * (elev_mult / es)))

def has_file(elev, frame, elev_mult, elev_step, steps_per_elev, out_dir, **seq):
    '''
    See if the sequence of files associated with the given elevation
    already exists. If any of the files is missing, all will be redone.
    '''
    e = e_to_filename(elev, elev_mult, elev_step)
    for i in range(frame, frame + steps_per_elev):
        f = os.path.join(out_dir, '{0}_{1}.ppm'.format(e, i))
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

def get_focus(x, y, z, region):
    '''
    Calculate the focus for nviz, which is metres in x/y/z, from the
    southwest corner of the image.
    '''
    blx = region[0]
    bly = region[3]
    return (x - blx, y - bly, z)

def vlen(v):
    return np.sqrt(np.sum(v**2))

def get_position(px, py, height, zexag, region):
    '''
    UV coordinates.
    gsd_real2model
    '''
    gs_unit_size = 1000.
    rng = 5. * gs_unit_size
    rng_off = 2. * gs_unit_size
    try:
        w = abs(region[2] - region[0])
        h = abs(region[3] - region[1])
        longdim = max(w, h)
        scale = gs_unit_size / longdim
        mx = (px - region[0]) * scale
        my = (region[1] - py) * scale
        dx = (mx + rng_off) / rng
        dy = (my + rng_off) / rng
        print(dx, dy)
    except:
        print(tb.format_exc())
    return dx, dy

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

def load_basemap(dem_file, dem_name, rgb_file, rgb_name, region, res, **kwargs):
    '''
    Load the dem and ortho, composite the ortho.
    '''

    print('Loading basemap:', dem_file, 'as', dem_name)

    # Load dem if not exists.
    try:
        os.system('r.in.gdal input={i} output={o}'.format(i = dem_file, o = dem_name))
    except: pass

    # Set computational region
    os.system('g.region raster={r}'.format(r = dem_name))

    print('Loading ortho:', rgb_file, 'as', rgb_name)

    # Load rgb if not exists
    try:
        os.system('r.in.gdal input={i} output={o}'.format(i = rgb_file, o = rgb_name))
        os.system('r.composite -c red={r}.red green={r}.green blue={r}.blue levels=32 output={r}'.format(r = rgb_name))
    except: pass

    print('Setting region to:', region)

    # Set the bounding region.
    cmd = 'g.region --overwrite --verbose align={r} n={n} w={w} e={e} s={s} res={res}'.format(
        r = dem_name, w = region[0], n = region[1], e = region[2], s = region[3], res = res
    )
    os.system(cmd)

def load_basin(tid, elev, basin_name, basin_name_elev, basin_dir, basin_pattern, elev_mult, elev_step, color_table, **kwargs):
    '''
    Load the basin raster and apply the colour table.
    Create an elevation layer to position it vertically.
    '''

    # Create an elevation model for the basin raster.
    basin_name = basin_name.format(id = tid)
    basin_name_elev = basin_name_elev.format(id = tid)
    e = e_to_filename(elev, elev_mult, elev_step)

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
        dem_name, rgb_name, position, focus, size, zexag, elev_mult, elev_step, steps_per_elev,
        dem_file, height, perspective, basin_ids, region, **kwargs):
    '''
    Render the 3D image in m.nviz.image.
    '''

    # Format filenames.
    basin_name = basin_name.format(id = tid)
    basin_name_elev = basin_name_elev.format(id = tid)
    e = e_to_filename(elev, elev_mult, elev_step)

    # Start, end position and focus.
    foc0 = get_focus(focus[0][0], focus[0][1], focus[0][2], region)
    foc1 = get_focus(focus[1][0], focus[1][1], focus[1][2], region)
    pos0 = get_position(*position[0], height, zexag, region)
    pos1 = get_position(*position[1], height, zexag, region)

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
        print('Pos', pos, '; Foc', foc)

        # Main part of the command.
        cmd = '''m.nviz.image -a --verbose --overwrite \
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
        cmd = cmd.format(**params)
        #print(cmd)
        os.system(cmd)

# Run me!
if __name__ == '__main__':
    run()
