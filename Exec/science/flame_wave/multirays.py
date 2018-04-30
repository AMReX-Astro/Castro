#!/usr/bin/env python3

import yt
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
import glob
import argparse

# Argument information
description = "Generates 1D vertical slices of axisymmetric datasets."
name_help = "The name of the main module."
datasets_help = "A list of datasets to be loaded by yt. Will be sorted by plot number by default."
out_help = "The desired output directory for the image files."
var_help = "The variable to plot against vertical position. Set to 'density' by default."
bounds_help = "The bounds for the vertical axis."
ext_help = "The extension of the file format to save to. PNG by default."
sort_help = """A floating point number specifying the digits to sort file names by. Digits preceding the decimal point
        give the starting index, digits following the decimal point give the number of characters. Make negative for
        descending order."""

# Construct parser and parse
parser = argparse.ArgumentParser(description=description)
parser.add_argument('__name__', help=name_help)
parser.add_argument('datasets', nargs='*', help=datasets_help)
parser.add_argument('-o', '--out', default='.', help=out_help)
parser.add_argument('-v', '--var', default='density', help=var_help)
parser.add_argument('-b', '--bounds', nargs=2, type=float, metavar=('UPPER', 'LOWER'), help=bounds_help)
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
field = args.var
bounds = args.bounds
colors = ('dodgerblue', 'springgreen', 'orangered')
labels = ('Left Boundary', 'Center', 'Right Boundary')

print("Generating...")

# Loop and generate
for i in range(0, len(ts)):

    print("Progress: {}%\r".format(int(100 * i / len(ts))), end='', flush=True)
    ds = ts[i]
    ad = ds.all_data()
    rmean, rmax = ad.mean('r'), ad.max('r')

    rays = ds.ortho_ray(1, (0, 0)), ds.ortho_ray(1, (rmean, 0)), ds.ortho_ray(1, (rmax, 0))
    z = rays[0]['z'], rays[1]['z'], rays[2]['z']
    iz = tuple(map(np.argsort, z))
    for z, iz, ray, label, color in zip(z, iz, rays, labels, colors):
        plt.plot(z[iz].d, ray[field][iz].d, label=label, color=color)

    plt.xlabel('z (cm)')
    plt.ylabel('{} ({})'.format(field, ad[field].units))
    plt.legend()

    plt.gca().set_yscale('log')
    if bounds is not None:
        plt.ylim(bounds)

    plt.savefig('{}/{}_{}_ray.{}'.format(args.out, ds, field, args.ext))
    plt.clf()

print("Task completed.")
