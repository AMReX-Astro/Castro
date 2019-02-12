#!/usr/bin/env python

desc = """ 
python driver to mencoder.  Take a list of files and make avi and mp4 movies.
./mkmovie.py [options] <list of PNG files>
options:
  -h, --help    : print this help page
  -f x, --fps=x : set frames per second to be x
  --double      : double each frame (effectively 2x slower)
  -N num        : create N copies of each frame (to really slow it down)
 
  -o prefix     : prefix for the movies (default: movie)
  --avi         : only generate an avi movie
  --mp4         : only generate an mp4 movie
  --endframes=N : add N duplicate frames to the end, to make the last 
                   frame persist
Note: for the mp4 movies, we duplicate the frames at the start to get
around a bug (?) in mencoder that skips some frames at the start.
M. Zingale (2013-02-19)
R. Orvedahl (2014-07-18)
	-add more options
	-change calling sequences
"""

import sys
import os
import getopt
import re

# routines for natural/human sorting
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    # list.sort(key=natural_keys) will sort according 
    # to human sorting (or natural sorting)
    return [ atoi(c) for c in re.split('(\d+)', text) ]


if len(sys.argv) == 1:
    print desc
    sys.exit(2)


prefix = "movie"
endFrames = 1
double = False
ncopies = 1
avi = True
mp4 = True
fps = 15

try: opts, next = getopt.getopt(sys.argv[1:], "ho:N:f:", 
                                ["double", "endframes=", "avi", "mp4", 
                                 "fps=", "help"])
except getopt.GetoptError:
    sys.exit("invalid calling sequence")

for o, a in opts:

    if o == "-o":
        prefix = a

    elif (o == "--fps" or o == "-f"):
        fps = a

    elif (o == "--help" or o == "-h"):
        print desc
        sys.exit(2)

    elif o == "--double":
        double = True

    elif o == "-N":
        ncopies = int(a)

    elif o == "--endframes":
        endFrames = int(a)

    elif o == "--avi":
        mp4 = False

    elif o == "--mp4":
        avi = False

try: frames = next[0:]
except IndexError:
    sys.exit("no frames specified")

print "number of frames: ", len(frames)

# sort frames
frames.sort(key=natural_keys)

# create temporary files that list the movie frames -- these are inputs
# to mencoder

if (avi):
    f = open("_mkmovie1.list", "w")

    for img in frames:
        f.write("%s\n" % (img))
        if double:
            f.write("%s\n" % (img))
        elif ncopies > 1:
            for i in range(ncopies):
                f.write("%s\n" % (img))

    if endFrames > 1:
        n = 0
        while (n < endFrames-1):
            f.write("%s\n" % (frames[len(frames)-1]))        
            n += 1

    f.close()


if (mp4):
    # for mp4, we want some extra frames at the start
    f = open("_mkmovie2.list", "w")

    n = 0
    while (n < 28):
        f.write("%s\n" % (frames[0]))
        n += 1

    for img in frames:
        f.write("%s\n" % (img))
        if double:
            f.write("%s\n" % (img))
        elif ncopies > 1:
            for i in range(ncopies):
                f.write("%s\n" % (img))

    if endFrames > 1:
        n = 0
        while (n < endFrames-1):
            f.write("%s\n" % (frames[len(frames)-1]))        
            n += 1

    f.close()



# make the movies
if (avi):
    str1 = "mencoder mf://@_mkmovie1.list -ovc lavc -lavcopts "
    str2 = "vcodec=msmpeg4v2:vbitrate=3000:vhq:mbd=2:trell "
    str3 = "-mf type=png:fps="+str(fps)
    str4 = " -o %s.avi" % (prefix)
    os.system(str1+str2+str3+str4)

if (mp4):
    str1 = "mencoder mf://@_mkmovie2.list -of lavf -lavfopts format=mp4 "
    str2 = "-ss 1 -ovc x264 -x264encopts crf=20.0:nocabac:level_idc=30:"
    str3 = "global_header:threads=2 -fps "+str(fps)
    str4 = " -o %s.mp4" % (prefix)
    os.system(str1+str2+str3+str4)

    str1 = "mencoder mf://@_mkmovie2.list -ovc x264 -x264encopts crf=10:"
    str2 = "me=umh:subq=9:nr=100:global_header -of lavf -lavfopts format=mp4 "
    str3 = "-fps "+str(fps)
    str4 = " -o %s_hg.mp4" % (prefix)
    os.system(str1+str2+str3+str4)


# remove the files
if (avi):
    os.remove("_mkmovie1.list")
if (mp4):
    os.remove("_mkmovie2.list")

if (avi):
    print "\n Avi Movie: "+prefix+".avi"
if (mp4):
    print "\n Mp4 Movie: "+prefix+".mp4"
    print "\n Mp4 Movie: "+prefix+"_hg.mp4"
