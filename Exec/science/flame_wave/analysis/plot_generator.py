#!/usr/bin/env python3

import yt
from yt.units import amu, cm

import os
import sys
import glob
import argparse
import numpy as np

from collections import namedtuple
from functools import reduce

# Argument information
description = """Generates plots of datasets using a specified yt plot function. Works with any slice or projection
        plot, as well as ParticlePlot."""
datasets_help = "A list of datasets to be loaded by yt. Will be sorted by plot number by default."
func_help = "The plotting function to use. SlicePlot by default."
out_help = "The desired output directory for the image files."
var_help = "The variable to plot. Set to 'Temp' by default."
bounds_help = "The bounds for the colorbar."
cmap_help = "The colormap for the variable to plot."
log_help = "If provided, sets the plot to a logarithmic scale."
time_help = "If provided, adds a timestamp to each plot with the given precision."
ext_help = "The extension of the file format to save to. PNG by default."
sort_help = """A floating point number specifying the digits to sort file names by. Digits preceding the decimal point
        give the starting index, digits following the decimal point give the number of characters. Make negative for
        descending order."""
quiver_help = """Overplots a vector field on each generated plot, taking the x and y components, the number of
        points to skip, and a scale factor for the arrows."""
contour_help = "Adds a contour map to the plot, with NCONT contours and limits specified by CLIM_LOW and CLIM_HI."
contour_opt_help = "Plot args to supply to the matplotlib contour function. Takes colors and linewidths."
xlim_help = "The x-axis limits."
ylim_help = "The y-axis limits."
stream_help = "Adds streamlines to the plot, showing the given vector field and a number of points to skip."
stream_color_help = """Options for coloring the streamlines - will be ignored if streamlines themselves were not
        requested. The first argument should be the field to color by, and later options should be the option name (odd
        indices) followed by the value (even indices). Valid options are cmap (the colormap) and display_threshold (to
        turn off the streamlines below a certain value for field color), and cbar (display streamline colorbar, True by
        default)."""
grid_help = """Add an overlay to the plot showing the hierarchical grid structure. May supply additional options to yt
        as whitespace-delimited pairs of the keyword and its value (e.g. --grid alpha 0.5 min_level 1 cmap gray)."""
cell_edges_help = "Overplot the edges of the grid cells."
flame_wave_help = """Populate settings with those for the flame wave project (colored streamlines with colorbar)."""

# Construct parser and parse
parser = argparse.ArgumentParser(description=description)
parser.add_argument('datasets', nargs='*', help=datasets_help)
parser.add_argument('-f', '--func', default='SlicePlot', help=func_help)
parser.add_argument('-o', '--out', default='', help=out_help)
parser.add_argument('-v', '--var', default='Temp', help=var_help)
parser.add_argument('-b', '--bounds', nargs=2, type=float, metavar=('LOWER', 'UPPER'), help=bounds_help)
parser.add_argument('-c', '--cmap', metavar=('NAME',), help=cmap_help)
parser.add_argument('--log', action='store_true', help=log_help)
parser.add_argument('-t', '--time', type=int, metavar=('PRECISION',), help=time_help)
parser.add_argument('-e', '--ext', type=lambda s: s.lower(), default='png', help=ext_help)
parser.add_argument('-s', '--sort', type=float, default=0.0, help=sort_help)
parser.add_argument('-q', '--quiver', nargs=4, metavar=('XFIELD', 'YFIELD', 'FACTOR', 'SCALE'),
        help=quiver_help)
parser.add_argument('-C', '--contour', nargs=4, metavar=('FIELD', 'NCONT', 'CLIM_LOW', 'CLIM_HI'), help=contour_help)
parser.add_argument('-Co', '--contour_opt', nargs=2, metavar=('COLORS', 'LINEWIDTHS'),
        help=contour_opt_help)
parser.add_argument('-x', '--xlim', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=xlim_help)
parser.add_argument('-y', '--ylim', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=ylim_help)
parser.add_argument('-S', '--stream', nargs=3, metavar=('XFIELD', 'YFIELD', 'FACTOR'))
parser.add_argument('-Sc', '--stream_color', nargs='*', help=stream_color_help)
parser.add_argument('--grid', nargs='*', default=None, help=grid_help)
parser.add_argument('--cell_edges', action='store_true', help=cell_edges_help)
parser.add_argument('--flame_wave', action='store_true', help=flame_wave_help)

args = parser.parse_args(sys.argv[1:])

coloropts = ['field_color', 'cmap', 'display_threshold', 'cbar']
ColorOpt = namedtuple('ColorOpt', field_names=coloropts)
optdict = dict(field_color=None, display_threshold=None, cmap=None, cbar=False)

if args.quiver is not None:

    args.quiver[2] = int(args.quiver[2])
    args.quiver[3] = float(args.quiver[3])

if args.stream is not None:

    args.stream[2] = int(args.stream[2])

if args.contour is not None:

    args.contour[1:] = list(map(float, args.contour[1:]))

    if args.contour_opt is not None:

        contour_opt = {'plot_args': {'colors': args.contour_opt[0], 'linewidths': int(args.contour_opt[1])}}

else:

    contour_opt = {}

if args.flame_wave:

    args.time = 4
    args.quiver = None
    args.stream = ['x_velocity', 'y_velocity', 16]
    optdict = dict(field_color='transvel', cmap='kamae', display_threshold=1.5e7, cbar=False)
    args.contour = ['enuc', 3, 1e16, 1e20]
    contour_opt = {'plot_args': {'colors': '0.7', 'linewidths': 1}}
    args.ylim = 0.375e4, 1.5e4
    args.xlim = 0.0, 5e4
    args.bounds = 0.0, 2e9

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

color_opt = ColorOpt(**optdict)

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
print(f"Successfully loaded the following files: {ts}\n")

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

        zctr = 0.5 * (zlim[0] + zlim[1])

    return xctr, yctr, zctr

# Derived fields
def _transvel(field, data):
    """ Transverse velocity """

    return np.sqrt(data['x_velocity']**2 + data['y_velocity']**2)

print("Generating...")

# Loop and generate
for ds in ts:

    if ("gas", "transvel") not in ds.field_list:

        try:
            ds.add_field(("gas", "transvel"), function=_transvel,
                    units="cm/s", sampling_type="cell")
        except: pass

    settings = {}
    settings['center'] = get_center(ds, args.xlim, args.ylim)
    settings['width'] = get_width(ds, args.xlim, args.ylim)
    if ds.geometry in {'cylindrical', 'spherical'}:
        settings['normal'] = 'theta'
        settings['origin'] = 'native'
    else:
        settings['normal'] = 'z'

    plot = func(ds, fields=field, **settings)
    if args.cmap: plot.set_cmap(field=field, cmap=args.cmap)

    if args.bounds is not None:
        plot.set_zlim(field, *args.bounds)

    plot.set_log(field, args.log)

    if args.time:

        time_format = f't = {{time:.{args.time}f}}{{units}}'

        plot.annotate_timestamp(corner='upper_left', time_format=time_format,
                time_unit='s', draw_inset_box=True, inset_box_args={'alpha': 0.0})

    if args.quiver is not None:
        plot.annotate_quiver(*args.quiver)

    if args.contour is not None:
        plot.annotate_contour(args.contour[0], ncont=args.contour[1], clim=args.contour[2:], **contour_opt)

    if args.stream is not None:

        kw = {'field_color': color_opt.field_color, 'display_threshold': color_opt.display_threshold}
        if color_opt.cbar: kw['add_colorbar'] = color_opt.cbar

        plot.annotate_streamlines(args.stream[0], args.stream[1], factor=args.stream[2],
                plot_args={'cmap': color_opt.cmap, 'arrowstyle': '->'}, **kw)

    if args.grid is not None:

        opts = dict(zip(args.grid[::2], args.grid[1::2]))
        plot.annotate_grids(**opts)

    if args.cell_edges:

        plot.annotate_cell_edges()

    suffix = args.func.replace('Plot', '').lower()
    plot.save(os.path.join(args.out, f'{ds}_{field}_{suffix}.{args.ext}'))
    print()

print("Task completed.")
