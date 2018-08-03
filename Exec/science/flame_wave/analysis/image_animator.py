#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib import animation

import sys
import argparse

# Argument information
description = "Animate a list of images and save the result as a video file."
name_help = "The name of the main module."
images_help = "The images to be animated. Must be supplied in order if sort option is not specified."
out_help = "The name of the output file. movie.mp4 by default."
dpi_help = """The desired dpi for the animation. The default is 256. Beware that larger numbers may cause memory errors
    on certain systems."""
sort_help = """A floating point number specifying the digits to sort file names by. Digits preceding the decimal point
    give the starting index, digits following the decimal point give the number of characters. Make negative for
    descending order."""

# Construct parser and parse
parser = argparse.ArgumentParser(description=description)
parser.add_argument('images', nargs='*', help=images_help)
parser.add_argument('-o', '--out', default='movie.mp4', help=out_help)
parser.add_argument('-d', '--dpi', type=int, default=256, help=dpi_help)
parser.add_argument('-s', '--sort', type=float, help=sort_help)

args = parser.parse_args(sys.argv[1:])
images = args.images

if not images:
    sys.exit("No images supplied for animation.")

# Sort if necessary
if args.sort is not None:
    # Descending or ascending
    desc = args.sort < 0
    # Staring index
    start = abs(int(args.sort))
    # Number of characters
    nchars = int(str(args.sort).split('.')[1])

    if nchars == 0:
        key = lambda img: img[start:]
    else:
        key = lambda img: img[start:start + nchars]
    images.sort(key=key, reverse=desc)

# Load images and animate
print("Reading...")
images = list(map(imread, images))

print("Animating...")
fig = plt.figure()
ax = plt.gca()
ax.set_axis_off()

imobj = ax.imshow(np.zeros(images[0].shape), origin='lower', alpha=1.0, zorder=1, aspect=1)

def init():
    imobj.set_data(np.zeros(images[0].shape))
    return imobj,

def animate(i):
    imobj.set_data(images[i][-1::-1])
    return imobj

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(images))
print("Saving...")
anim.save(args.out, dpi=args.dpi)
