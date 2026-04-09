#!/usr/bin/env python

import sys
import os
import argparse
import tempfile as tf

# Argument information
description = "Animate a list of images and save the result as a video file."
images_help = "The images to be animated. Must be supplied in order if sort option is not specified."
out_help = "The name of the output file. movie.mp4 by default."
sort_help = """A floating point number specifying the digits to sort file names by. Digits preceding
    the decimal point give the starting index, digits following the decimal point give the number of
    characters. Make negative for descending order."""
framerate_help = "Output fps. Set to 30 by default."

# Construct parser and parse
parser = argparse.ArgumentParser(description=description)
parser.add_argument('images', nargs='*', help=images_help)
parser.add_argument('-o', '--out', default='movie', help=out_help)
parser.add_argument('-r', '--framerate', type=int, default=30, help=framerate_help)
parser.add_argument('-s', '--sort', nargs="?", type=float, default=0.0, help=sort_help)

args = parser.parse_args(sys.argv[1:])
images = args.images

# Sort if necessary
if args.sort is not None:

    # Descending or ascending
    desc = args.sort < 0
    # Starting index
    start = abs(int(args.sort))
    # Number of characters
    nchars = int(str(args.sort).split('.')[1])

    if nchars == 0:
        key = lambda filename: filename[start:]
    else:
        key = lambda filename: filename[start:start+nchars]
    images.sort(key=key, reverse=desc)

if args.out.endswith(".mp4"):
    args.out = args.out[:-4]
if not args.out:
    sys.exit("Invalid output file!")

with tf.TemporaryDirectory() as dir:

    # The files are numbered
    # Padding with zeros preserves sort order
    ndigits = len(str(len(images)))
    baselink = "__fftemp_{:0%dd}.png" % ndigits

    for i, image in enumerate(map(os.path.abspath, images)):

        # Create a symlink for each file
        linkname = baselink.format(i)
        linkname = os.path.join(dir, linkname)
        os.symlink(image, linkname)

    cmd = 'ffmpeg -y -r {} -pattern_type glob -i \'{}\' -c:v libx264 -pix_fmt yuv420p \'{}.mp4\''
    cmd = cmd.format(args.framerate, os.path.join(dir, '__fftemp_*.png'), args.out)
    os.system(cmd)

print("Task completed.")
