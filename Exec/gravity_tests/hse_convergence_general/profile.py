#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import argparse

import numpy as np
import matplotlib.pyplot as plt

import yt


def doit(pfile=None):

    fig, ax = plt.subplots(3,1)

    for pf in pfile:
        print(f"about to read {pfile}")
        ds = yt.load(pf, hint="castro")

        label = pf.split("plt")[-1]

        ray = ds.ray((0, 0, 0), (2.56e3, 0, 0))
        isrt = np.argsort(ray[("boxlib", "x")])

        ax[0].plot(ray[("boxlib", "x")][isrt], ray[("boxlib", "density")][isrt], label=label)
        ax[0].set_ylabel("density")
        ax[0].set_yscale("log")

        ax[1].plot(ray[("boxlib", "x")][isrt], ray[("boxlib", "Temp")][isrt], label=label)
        ax[1].set_ylabel("temperature")
        ax[1].set_yscale("log")

        ax[2].plot(ray[("boxlib", "x")][isrt], ray[("boxlib", "magvel")][isrt], label=label)
        ax[2].set_ylabel("|U|")
        ax[2].set_xlabel("x")
        ax[2].set_yscale("log")

    ax[0].legend(frameon=False)
    fig.tight_layout()
    fig.set_size_inches(9, 12)
    fig.savefig(f"{pfile[-1]}.png")

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("plotfile", metavar="plotfile name", type=str, nargs="*")
    args = p.parse_args()
    print(args.plotfile)
    doit(pfile=args.plotfile)
