#!/usr/bin/env python3

'''
This script adds labels to frames generated by the nviz.py script.
'''

import sys
import os
import math
import cairo
import json
import multiprocessing as mp

font = '/usr/share/fonts/type1/gsfonts/p052003l.pfb'

# Panel size.
height = 50
width = 300
pad = 10

font_size = 24

def load_colours(colour_file):
    '''
    Load the colours from the GRASS colour map.
    Convert to float values.
    '''
    colours = {}
    with open(colour_file, 'r') as f:
        try:
            while True:
                i, c = f.readline().strip().split()
                r, g, b = list(map(int, c.split(':')))
                colours[int(i)] = (r / 255., g / 255., b / 255.)
        except: pass
    return colours

def basin_panel(ctx, name, width, height, pad, colour):
    '''
    Print a legend item for a basin with the colour and name.
    '''

    # Background.
    ctx.move_to(0, 0)
    ctx.rectangle(0, 0, width, height)
    ctx.set_source_rgba(1, 1, 1, 0.75)
    ctx.fill()

    # Colour swatch.
    ctx.rectangle(pad, pad, height - pad * 2, height - pad * 2)
    ctx.set_source_rgb(*colour)
    ctx.fill();

    # Stroke border.
    ctx.rectangle(pad, pad, height - pad * 2, height - pad * 2)
    ctx.set_source_rgb(0, 0, 0)
    ctx.set_line_width(1)
    ctx.stroke()

    # Text
    ctx.move_to(height, height / 2 + pad)
    ctx.set_font_size(font_size)
    ctx.set_source_rgb(0, 0, 0)
    ctx.show_text(name)

def elev_panel(ctx, elev, width, height, pad):
    '''
    Add a legend item for the elevation.
    '''

    # Background.
    ctx.move_to(0, 0)
    ctx.rectangle(0, 0, width, height)
    ctx.set_source_rgba(1, 1, 1, 0.75)
    ctx.fill()

    # Format the label, get its extent.
    t = '{}: {:.2f}m'.format('Elevation', elev)
    te = ctx.text_extents(t)

    # Render text in center.
    ctx.move_to(width / 2 - te.width / 2, height / 2 + pad)
    ctx.set_font_size(font_size)
    ctx.set_source_rgb(0, 0, 0)
    ctx.show_text(t)

def text_width(ctx, ids, panels):
    '''
    Get the maximum text width for each of the panel names, given the id list.
    '''
    ctx.set_font_size(font_size)
    w = width
    for id_ in ids:
        w = max(w, ctx.text_extents(panels[id_]['name']).width)
    return w

def job(infiles, outdir, colour_file, ids, panels):
    '''
    Process the image.
    '''
    # Load the colours
    colours = load_colours(colour_file)

    for infile in infiles:

        print('Processing', infile)

        # Rearrange numbers to make ordering easier.
        parts = os.path.splitext(os.path.basename(infile))[0].split('_')
        outfile = os.path.join(outdir, '{}_{}.{}'.format(parts[0], parts[1], 'png'))
        if os.path.exists(outfile) and os.stat(outfile).st_size > 0:
            continue

        # Get the elevation from the filename.
        f = os.path.splitext(os.path.basename(infile))[0]
        parts = f.split('_') # TODO: Match the pattern from nviz.
        e = float(parts[1]) / 1000.

        # If the image is not a png convert it as a tmp.
        tmp = '/tmp/overlay_{}.png'.format(os.getpid())
        if not infile.lower().endswith('.png'):
            os.system('convert {} {}'.format(infile, tmp))
            infile = tmp

        # Initialize the image and context.
        inimg = cairo.ImageSurface.create_from_png(infile)
        ctx = cairo.Context(inimg)

        # Get the maximum text width.
        tw = text_width(ctx, ids, panels) + height + pad

        # Translate to start the panel drawing in the bottom left.
        ctx.translate(height / 2, inimg.get_height() - (len(ids) + 1) * height - height / 2)

        # Do the legend panels.
        for id_ in ids:
            name = panels[id_]['name']
            colour = colours[id_]
            basin_panel(ctx, name, tw, height, pad, colour)
            ctx.translate(0, height)

        # Do the elevation panel.
        elev_panel(ctx, e, tw, height, pad)

        # Write.
        inimg.write_to_png(outfile)

def run(indir, outdir, colour_file, bid):
    panels = {}
    ids = None
    with open('config.json', 'r') as f:
        config = json.loads(f.read())
        for conf in config['configs']:
            if conf['id'] == bid:
                ids = conf['basin_ids']
                break
        for conf in config['configs']:
            if conf['id'] in ids:
                for k, v in config.items():
                    if k != 'configs':
                        conf[k] = v
                panels[conf['id']] = conf

    print('Spawning processes')

    mp.set_start_method('spawn')

    files = []

    for f in [x for x in os.listdir(indir) if x.endswith('.ppm')]:
        files.append(os.path.join(indir, f))

    ts = 4
    procs = []
    for t in range(ts):
        files_ = files[t * int(len(files) / 4):(t+1) * int(len(files) / 4)]
        procs.append(mp.Process(target = job, args = (files_, outdir, colour_file, ids, panels)))
        procs[t].start()

    for t in range(ts):
        procs[t].join()

if __name__ == '__main__':
    try:
        in_file = sys.argv[1]
        out_file = sys.argv[2]
        colour_file = sys.argv[3]
        bid = int(sys.argv[4])

        run(in_file, out_file, colour_file, bid)

    except Exception as e:
        print(e)
        print('Usage: <in file> <out file> <colour file> <basin id,[basin id,[basin id,...]]>')
        sys.exit(1)
