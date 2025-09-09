#!/usr/bin/python

from numpy import *

def read_gnu_file(filenm):
    x = []
    y = []
    f = open(filenm)
    line = f.readline()
    t = float(line.split('"')[1].split('=')[2])
    for line in f.readlines():
        if not line[0] == ";":
            words = line.split()
            x.append(float(words[0]))
            y.append(float(words[1]))
    f.close()
    return array(y), array(x), t
