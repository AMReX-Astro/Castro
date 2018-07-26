#!/usr/bin/env python3

import yt
from yt.units import cm
from yt import derived_field

import os
import sys
import glob
import argparse
import numpy as np

from collections import namedtuple

# Argument information
description = """Generates plots of datasets using a specified yt plot function. Works with any slice or projection
        plot, as well as ParticlePlot."""
name_help = "The name of the main module."
datasets_help = "A list of datasets to be loaded by yt. Will be sorted by plot number by default."
func_help = "The plotting function to use. SlicePlot by default."
out_help = "The desired output directory for the image files."
var_help = "The variable to plot. Set to 'Temp' by default."
bounds_help = "The bounds for the colorbar."
log_help = "If provided, sets the plot to a logarithmic scale."
time_help = "If provided, adds a timestamp to each plot with the given precision."
ext_help = "The extension of the file format to save to. PNG by default."
sort_help = """A floating point number specifying the digits to sort file names by. Digits preceding the decimal point
        give the starting index, digits following the decimal point give the number of characters. Make negative for
        descending order."""
quiver_help = """Overplots a vector field on each generated plot, taking the x and y components, the number of
        points to skip, and a scale factor for the arrows."""
contour_help = "Adds a contour map to the plot."
xlim_help = "The x-axis limits."
ylim_help = "The y-axis limits."
stream_help = "Adds streamlines to the plot, showing the given vector field and a number of points to skip."
stream_color_help = """Options for coloring the streamlines - will be ignored if streamlines themselves were not
        requested. The first argument should be the field to color by, and later options should be the option name (odd
        indices) followed by the value (even indices). Valid options are cmap (the colormap) and display_threshold (to
        turn off the streamlines below a certain value for field color), and cbar (display streamline colorbar, True by
        default)."""
flame_wave_help = """Populate settings with those for the flame wave project (colored streamlines with colorbar)."""

# Construct parser and parse
parser = argparse.ArgumentParser(description=description)
parser.add_argument('datasets', nargs='*', help=datasets_help)
parser.add_argument('-f', '--func', default='SlicePlot', help=func_help)
parser.add_argument('-o', '--out', default='', help=out_help)
parser.add_argument('-v', '--var', default='Temp', help=var_help)
parser.add_argument('-b', '--bounds', nargs=2, type=float, metavar=('LOWER', 'UPPER'), help=bounds_help)
parser.add_argument('--log', action='store_true', help=log_help)
parser.add_argument('-t', '--time', type=int, default=-1, metavar=('PRECISION',), help=time_help)
parser.add_argument('-e', '--ext', type=lambda s: s.lower(), default='png', help=ext_help)
parser.add_argument('-s', '--sort', type=float, default=0.0, help=sort_help)
parser.add_argument('-q', '--quiver', nargs=4, metavar=('XFIELD', 'YFIELD', 'FACTOR', 'SCALE'),
        help=quiver_help)
parser.add_argument('-c', '--contour', metavar=('FIELD',), help=contour_help)
parser.add_argument('-x', '--xlim', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=xlim_help)
parser.add_argument('-y', '--ylim', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=ylim_help)
parser.add_argument('-S', '--stream', nargs=3, metavar=('XFIELD', 'YFIELD', 'FACTOR'))
parser.add_argument('-Sc', '--stream_color', nargs='*', help=stream_color_help)
parser.add_argument('--flame_wave', action='store_true', help=flame_wave_help)

args = parser.parse_args(sys.argv[1:])

coloropts = ['field_color', 'cmap', 'display_threshold', 'cbar']
ColorOpt = namedtuple('ColorOpt', field_names=coloropts)
optdict = {'field_color': None, 'display_threshold': None, 'cmap': None, 'cbar': False}

contour_opt = {}
color_opt = None

if args.quiver is not None:

    args.quiver[2] = int(args.quiver[2])
    args.quiver[3] = float(args.quiver[3])

if args.stream is not None:

    args.stream[2] = int(args.stream[2])

if args.flame_wave:

    args.time = 4
    args.quiver = None
    args.stream = ['x_velocity', 'y_velocity', 16]
    color_opt = ColorOpt(field_color='transvel', cmap='kamae', display_threshold=1.5e7, cbar=True)
    args.contour = 'enuc'
    contour_opt = {'ncont': 3, 'clim': (1e16, 1e20), 'plot_args': {"colors": "0.7", "linewidths": 1}}
    args.ylim = (0.375e4, 1.5e4)
    args.xlim = (0.0, 5e4)

if args.stream_color is not None:

    opts = args.stream_color[1::2]
    vals = args.stream_color[2::2]
    for opt, val in zip(opts, vals):
        optdict[opt] = val

    optdict['field_color'] = args.stream_color[0]
    if isinstance(optdict['cbar'], str):
        optdict['cbar'] = eval(optdict['cbar'].capitalize())
    if optdict['display_threshold'] is not None:
        optdict['display_threshold'] = float(optdict['display_threshold'])

if color_opt is None: color_opt = ColorOpt(**optdict)

# Make output directory
if not args.out:
    args.out = os.getcwd()
if not os.path.exists(args.out):
    os.makedirs(args.out)

# Grab files from working directory if none were specified
ts = args.datasets
if not ts:
    ts = glob.glob('plt*')

# Exit if nothing could be loaded
if len(ts) < 1:
    sys.exit("No files were available to be loaded.")

# Sort and load files
desc = args.sort < 0
start = abs(int(args.sort))
nchars = int(str(args.sort).split('.')[1])

if nchars == 0:
    key = lambda fname: fname[start:]
else:
    key = lambda fname: fname[start:start + nchars]
ts.sort(key=key, reverse=desc)

tf = lambda file: yt.load(file.rstrip('/'))
ts = list(map(tf, ts))
print("Successfully loaded the following files: {}\n".format(ts))

# Generate plots
func = getattr(yt, args.func)
field = args.var

def get_width(ds, xlim=None, ylim=None, zlim=None):
    """ Get the width of the plot. """

    if xlim is None: xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
    else: xlim = xlim[0] * cm, xlim[1] * cm

    if ylim is None: ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
    else: ylim = ylim[0] * cm, ylim[1] * cm

    xwidth = (xlim[1] - xlim[0]).in_cgs()
    ywidth = (ylim[1] - ylim[0]).in_cgs()

    if ds.domain_dimensions[2] == 1:
        zwidth = 0.0
    else:
        if zlim is None: zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
        else: zlim = zlim[0] * cm, zlim[1] * cm

        zwidth = (zlim[1] - zlim[0]).in_cgs()

    return xwidth, ywidth, zwidth

def get_center(ds, xlim=None, ylim=None, zlim=None):
    """ Get the coordinates of the center of the plot. """

    if xlim is None: xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
    else: xlim = xlim[0] * cm, xlim[1] * cm

    if ylim is None: ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
    else: ylim = ylim[0] * cm, ylim[1] * cm

    xctr = 0.5 * (xlim[0] + xlim[1])
    yctr = 0.5 * (ylim[0] + ylim[1])

    if ds.domain_dimensions[2] == 1:
        zctr = 0.0
    else:
        if zlim is None: zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
        else: zlim = zlim[0] * cm, zlim[1] * cm

        zctr = 0.5 * (zlim[2] + zlim[2])

    return xctr, yctr, zctr

@derived_field(name='transvel', units='cm/s', sampling_type='cell')
def _transvel(field, data):
    """ Transverse velocity """

    return np.sqrt(data['x_velocity']**2 + data['y_velocity']**2)

print("Generating...")

# Loop and generate
for ds in ts:

    settings = {}
    settings['center'] = get_center(ds, args.xlim, args.ylim)
    settings['width'] = get_width(ds, args.xlim, args.ylim)
    if ds.geometry == 'cylindrical':
        settings['normal'] = 'theta'
        settings['origin'] = 'native'

    plot = func(ds, fields=field, **settings)

    if args.bounds is not None:
        plot.set_zlim(field, *args.bounds)
    plot.set_log(field, args.log)
    if args.time > 0:
        plot.annotate_timestamp(corner='upper_left', time_format='t = {{time:.{}f}}{{units}}'.format(args.time),
                time_unit='s', draw_inset_box=True, inset_box_args={'alpha': 0.0})

    if args.quiver is not None:
        plot.annotate_quiver(*args.quiver)

    if args.contour is not None:
        plot.annotate_contour(args.contour, **contour_opt)

    if args.stream is not None:

        plot.annotate_streamlines(args.stream[0], args.stream[1], factor=args.stream[2],
                field_color=color_opt.field_color, add_colorbar=color_opt.cbar,
                display_threshold=color_opt.display_threshold,
                plot_args={'cmap': color_opt.cmap, 'arrowstyle': '->'})

    suffix = args.func.replace('Plot', '').lower()
    plot.save(os.path.join(args.out, '{}_{}_{}.{}'.format(ds, field, suffix, args.ext)))
    print()

print("Task completed.")
