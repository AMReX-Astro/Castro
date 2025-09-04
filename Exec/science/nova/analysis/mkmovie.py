#!/usr/bin/env python3

import argparse
import os.path
import re
import shlex
import subprocess
import tempfile
from collections import defaultdict
from typing import Union


# routines for natural/human sorting
def atoi(text: str) -> str | int:
    return int(text) if text.isdigit() else text


def natural_keys(text: str) -> tuple[str | int, ...]:
    # list.sort(key=natural_keys) will sort according
    # to human sorting (or natural sorting)
    return tuple(atoi(c) for c in re.split(r"(\d+)", text))


parser = argparse.ArgumentParser(
    description="Python driver for ffmpeg. Takes a list of files and makes AVI and MP4 movies."
)
parser.add_argument(
    "-f",
    "--fps",
    metavar="x",
    type=int,
    default=15,
    help="set frames per second to be x",
)
copy_parser = parser.add_mutually_exclusive_group()
copy_parser.add_argument(
    "--double",
    dest="ncopies",
    action="store_const",
    const=2,
    default=1,
    help="double each frame (effectively 2x slower)",
)
copy_parser.add_argument(
    "-N",
    metavar="num",
    dest="ncopies",
    type=int,
    default=1,
    help="create num copies of each frame (to really slow it down)",
)
parser.add_argument(
    "-o",
    "--prefix",
    default="movie",
    help="prefix for the movie files",
)
type_parser = parser.add_mutually_exclusive_group()
type_parser.add_argument(
    "--avi",
    dest="mp4",
    action="store_false",
    default=True,
    help="only generate an avi movie",
)
type_parser.add_argument(
    "--mp4",
    dest="avi",
    action="store_false",
    default=True,
    help="only generate an mp4 movie",
)
parser.add_argument(
    "--endframes",
    metavar="N",
    type=int,
    default=0,
    help="add N duplicate frames to the end, to make the last frame persist",
)
parser.add_argument(
    "--startframes",
    metavar="N",
    type=int,
    default=0,
    help="add N duplicate frames to the start",
)
parser.add_argument(
    "-y", dest="yes", action="store_true", help="overwrite output files"
)

parser.add_argument(
    "frames", nargs="+", help='list of PNG files (will be sorted by "natural" order)'
)

args = parser.parse_args()

if "/" in args.prefix:
    os.makedirs(os.path.dirname(args.prefix), exist_ok=True)

# print("Got order:")
# print("\n".join("  " + x for x in args.frames))
# print()
args.frames.sort(key=natural_keys)
# print("Sorted:")
# print("\n".join("  " + x for x in args.frames))


def run_ffmpeg(frame_file: str, suffix: str, codec_opts: list[str]) -> str:
    """Returns the output movie file name."""
    common_opts = [
        "-r",
        str(args.fps),
        "-f",
        "concat",
        "-i",
        frame_file,
        "-framerate",
        str(args.fps),
    ]
    if args.yes:
        common_opts.insert(0, "-y")
    output_name = f"{args.prefix}{suffix}"
    subprocess.run(
        ["ffmpeg", *common_opts, *codec_opts, output_name], check=True
    )
    return output_name


def format_frame(path: str) -> str:
    return f"file {shlex.quote(path)}\nduration {1.0 / args.fps}\n"


with tempfile.NamedTemporaryFile("w", dir=".") as frame_list:
    frame_list.write(format_frame(args.frames[0]) * args.startframes)
    for frame in args.frames:
        frame_list.write(format_frame(frame) * args.ncopies)
    frame_list.write(format_frame(args.frames[-1]) * args.endframes)
    frame_list.flush()

    # encoder options: `ffmpeg -h encoder=<-c:v argument>`
    movies = defaultdict(list)
    if args.avi:
        movies["avi"].append(
            run_ffmpeg(
                frame_list.name,
                ".avi",
                "-c:v msmpeg4v2 -b:v 3000K -mbd 2 -trellis 1".split(),
            )
        )
    if args.mp4:
        movies["mp4"].append(
            run_ffmpeg(
                frame_list.name,
                ".mp4",
                "-c:v libx264 -pix_fmt yuv420p -crf 20 -coder vlc -level 30 -flags global_header -threads 2".split(),
            )
        )
        movies["mp4"].append(
            run_ffmpeg(
                frame_list.name,
                "_hg.mp4",
                "-c:v libx264 -pix_fmt yuv420p -crf 10 -me_method umh -subq 9".split(),
            )
        )
        # movies["mp4"].append(
        #     run_ffmpeg(
        #         frame_list.name,
        #         "_nvenc.mp4",
        #         "-vsync 0 -c:v h264_nvenc -pix_fmt yuv420p -preset:v p7 -tune:v hq -rc:v vbr -cq:v 19 -b:v 0 -profile:v high".split(),
        #     )
        # )

for k, v in movies.items():
    for movie in v:
        print(f"\n {k.title()} Movie: {movie}")