#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib import animation

import sys
import argparse

# Argument information
description = "Animate a list of images and save the result as a video file."
images_help = "The images to be animated. Must be supplied in order if sort option is not specified."
out_help = "The name of the output file. movie.mp4 by default."
dpi_help = """The desired dpi for the animation. The default is 256. Beware of memory errors when using larger values."""
stack_help = """Create a single plot with multiple stacked subplots. The first argument is the number of subplots, and
        the second determines how they are grouped (set to 1 to stack each NSUBPLOTS images, and 0 to split the image
        list into NSUBPLOTS parts, and cycle through them in parallel). The input will be sorted prior to splitting the
        the list in the first case, and afterward in the latter one."""
sort_help = """A floating point number specifying the digits to sort file names by. Digits preceding the decimal point
    give the starting index, digits following the decimal point give the number of characters. Make negative for
    descending order."""

# Construct parser and parse
parser = argparse.ArgumentParser(description=description)
parser.add_argument('images', nargs='*', help=images_help)
parser.add_argument('-o', '--out', default='movie.mp4', help=out_help)
parser.add_argument('-d', '--dpi', type=int, default=256, help=dpi_help)
parser.add_argument('--stack', nargs=2, type=int, metavar=('NSUBPLOTS', 'GROUPMODE'), help=stack_help)
parser.add_argument('-s', '--sort', type=float, help=sort_help)

args = parser.parse_args(sys.argv[1:])
images = args.images

if not images:
    sys.exit("No images supplied for animation.")

if args.stack is not None:

    nsubplots = args.stack[0]
    groupmode = args.stack[1]

else:

    nsubplots = 1
    groupmode = 0

sublist_size = len(images) // nsubplots

if nsubplots > 1 and groupmode == 0:

    step = sublist_size
    sublists = [images[i-step:i] for i in range(step, len(images) + 1, step)]

else:

    sublists = [images]

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
    for imlist in sublists:
        imlist.sort(key=key, reverse=desc)

if nsubplots > 1 and groupmode == 1:

    sublists = [images[i::nsubplots] for i in range(nsubplots)]

# Load images and animate
print("Reading...")
sublists = [list(map(imread, images)) for images in sublists]

print("Animating...")

if nsubplots == 1:

    fig = plt.figure()
    axes = [plt.gca()]

else:

    fig, axes = plt.subplots(nsubplots)

plt.margins(0, 0)
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)

imobj = []

for ax, sub in zip(axes, sublists):

    ax.set_axis_off()

    zeros = np.zeros(sub[0].shape)
    imobj.append(ax.imshow(zeros, origin='lower', alpha=1.0, zorder=1, aspect=1))

def init():

    for im, sub in zip(imobj, sublists):
        im.set_data(np.zeros(sub[0].shape))
    return imobj,

def animate(i):

    for im, sub in zip(imobj, sublists):
        im.set_data(sub[i][-1::-1])
    return imobj

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=sublist_size)
anim.save(args.out, dpi=args.dpi)
