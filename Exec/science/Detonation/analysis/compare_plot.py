#!/usr/bin/env python3

from __future__ import print_function

import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
import sys
import numpy as np


def get_Te_profile(plotfile):

    ds = yt.load(plotfile)

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])
    x_coord = np.array(ad['x'][srt])
    temp = np.array(ad['Temp'][srt])
    enuc = np.array(ad['enuc'][srt])

    return time, x_coord, temp, enuc


def doit(plotfiles):

    f = plt.figure()
    f.set_size_inches(7.0, 9.0)

    ax_T = f.add_subplot(211)
    ax_e = f.add_subplot(212)
    
    for p in plotfiles:
        if plotfiles[p] is None:
            continue

        time, x, T, enuc = get_Te_profile(plotfiles[p])

        ax_T.plot(x, T, label=p)
        ax_e.plot(x, enuc)
        
    ax_T.legend(frameon=False)

    ax_e.set_yscale("log")
    
    f.savefig("det.png")
    

if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--strang", help="plotfile with Strang splitting", type=str, default=None)
    p.add_argument("--sdc", help="plotfile with SDC", type=str, default=None)
    p.add_argument("--mol", help="plotfile with MOL integration", type=str, default=None)

    args = p.parse_args()

    pfiles = {}
    pfiles["Strang"] = args.strang
    pfiles["SDC"] = args.sdc
    pfiles["MOL"] = args.mol
    
    doit(pfiles)
