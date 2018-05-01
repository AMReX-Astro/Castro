#!/usr/bin/env python3

import yt

import os
import sys
import glob
import argparse

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

# Construct parser and parse
parser = argparse.ArgumentParser(description=description)
parser.add_argument('__name__', help=name_help)
parser.add_argument('datasets', nargs='*', help=datasets_help)
parser.add_argument('-f', '--func', default='SlicePlot', help=func_help)
parser.add_argument('-o', '--out', help=out_help)
parser.add_argument('-v', '--var', default='Temp', help=var_help)
parser.add_argument('-b', '--bounds', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=bounds_help)
parser.add_argument('--log', action='store_true', help=log_help)
parser.add_argument('-t', '--time', type=int, default=-1, help=time_help)
parser.add_argument('-e', '--ext', type=lambda s: s.lower(), default='png', help=ext_help)
parser.add_argument('-s', '--sort', type=float, default=3.0, help=sort_help)

args = parser.parse_args(sys.argv)

# Make output directory
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
    key = lambda img: img[start:]
else:
    key = lambda img: img[start:start + nchars]
ts.sort(key=key, reverse=desc)

tf = lambda file: yt.load(file.rstrip('/'))
ts = list(map(tf, ts))
print("Successfully loaded the following files: {}\n".format(ts))

# Generate plots
func = getattr(yt, args.func)
field = args.var

def get_width(ds):
    """ Get the width of the plot. """

    xwidth = (ds.domain_right_edge[0] - ds.domain_left_edge[0]).in_cgs()

    ywidth = (ds.domain_right_edge[1] - ds.domain_left_edge[1]).in_cgs()

    if ds.domain_dimensions[2] == 1:
        zwidth = 0.0
    else:
        zwidth = (ds.domain_right_edge[2] - ds.domain_left_edge[2]).in_cgs()

    return xwidth, ywidth, zwidth

def get_center(ds):
    """ Get the coordinates of the center of the plot. """

    xctr = 0.5 * (ds.domain_left_edge[0] + ds.domain_right_edge[0])
    yctr = 0.5 * (ds.domain_left_edge[1] + ds.domain_right_edge[1])

    if ds.domain_dimensions[2] == 1:
        zctr = 0.0
    else:
        zctr = 0.5 * (self.file_info.ds.domain_left_edge[2] + self.file_info.ds.domain_right_edge[2])

    return xctr, yctr, zctr

print("Generating...")

# Loop and generate
for ds in ts:

    settings = {}
    settings['center'] = get_center(ds)
    settings['width'] = get_width(ds)
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

    suffix = args.func.replace('Plot', '').lower()
    plot.save('{}/{}_{}_{}.{}'.format(args.out, ds, field, suffix, args.ext))
    print()

print("Task completed.")
