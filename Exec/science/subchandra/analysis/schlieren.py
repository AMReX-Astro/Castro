#!/usr/bin/env python3

import matplotlib
matplotlib.use('agg')

import argparse
import os
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce

from mpl_toolkits.axes_grid1 import ImageGrid

# assume that our data is in CGS
from yt.units import cm, amu
from yt.frontends.boxlib.api import CastroDataset

#times = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]
times = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
#times = [0.0, 0.15, 0.3, 0.45]

clip_val = -35.0
max_val = -19.0

domain_frac = 0.15

def find_files(plist):

    mask = np.zeros(len(times))
    files_to_plot = []
    for pfile in plist:
        for k, t in enumerate(times):
            if mask[k]:
                continue
            print(pfile)
            ds = CastroDataset(pfile)
            if ds.current_time >= t:
                files_to_plot.append(pfile)
                mask[k] = 1.0

    return files_to_plot

def schlieren(data, cutoff=-16):
    """
    Calculates the Schlieren S = ln|(del^2 rho) / rho|. The cutoff is used to
    filter out low amplitude noise.
    """
    dens = data.field_data[('gas', 'density')]

    nr, nz, nl = dens.shape
    dens = dens[:, :, 0]

    coords = data.fcoords

    r = coords[:, 0].reshape((nr, nz, nl))[:, :, 0]
    z = coords[:, 1].reshape((nr, nz, nl))[:, :, 0]

    lapl = np.zeros_like(dens)

    # assume these are constant, which is at least true for this data
    dr = r[5, 5] - r[4, 5]
    dz = z[5, 5] - z[5, 4]

    rl = r - 0.5 * dr
    rr = r + 0.5 * dr

    # calculate the laplacian

    # radial term -- this is 1/r d/dr (r  dA/dr)

    lapl[1:-1, :] = 1 / (r[1:-1, :] * dr**2) * (
        -(rl[1:-1,:] + rr[1:-1,:]) * dens[1:-1:, :] +
        rl[1:-1, :] * dens[:-2, :] + rr[1:-1, :] * dens[2:, :])

    lapl[:, 1:-1] += 1 / dz**2 * \
        (dens[:, 2:] + dens[:, :-2] - 2 * dens[:, 1:-1])

    lapl[:, :] /= dens

    S = np.log(np.abs(lapl))
    print(S.min(), S.max())
    S[S < cutoff] = cutoff

    return S

def doit(pfiles):

    print("looking to plot: ", pfiles)

    fig = plt.figure()

    if len(pfiles) > 4:
        nrows = 2
        ncols = (len(pfiles) + 1)//2
    else:
        nrows = 1
        ncols = len(pfiles)

    grid = ImageGrid(fig, 111, nrows_ncols=(nrows, ncols),
                     axes_pad=0.75, cbar_pad=0.1, label_mode="L", cbar_mode="single")

    for i in range(nrows * ncols):

        if i < len(pfiles):
            pf = pfiles[i]
        else:
            grid[i].remove()
            continue

        ds = CastroDataset(pf)

        level = ds.index.max_level
        ref = int(np.product(ds.ref_factors[0:level]))

        dims = ds.domain_dimensions * ref

        sq = ds.covering_grid(level, left_edge=[0, 0, 0],
                              dims=dims, fields=["density"])

        data = schlieren(sq, cutoff=clip_val).T

        print(type(data))
        nx = data.shape[0]
        ny = data.shape[1]

        # note that things are transposed from what we normally think
        iymin = 0
        iymax = int(domain_frac * ny)
        ixmin = nx//2 - int(0.5 * domain_frac * nx)
        ixmax = nx//2 + int(0.5 * domain_frac * nx)

        yctr = 0.5 * (ds.domain_left_edge[1] + ds.domain_right_edge[1])
        dy = ds.domain_right_edge[1] - ds.domain_left_edge[1]

        xmin = yctr - 0.5 * domain_frac * dy
        xmax = yctr + 0.5 * domain_frac * dy

        ymin = 0.0
        ymax = domain_frac * ds.domain_right_edge[0]

        im = grid[i].axes.imshow(data[ixmin:ixmax+1, iymin:iymax+1],
                                 vmin=clip_val, vmax=max_val,
                                 extent=[ymin, ymax, xmin, xmax],
                                 origin='lower', cmap=plt.cm.bone_r)

        if i == 0:
            grid.cbar_axes[0].colorbar(im)

        grid[i].axes.text(0.05, 0.05, f"time = {float(ds.current_time):8.3f} s",
                          transform=grid[i].axes.transAxes)

        grid[i].axes.set_xlabel(r'$r$')
        grid[i].axes.set_ylabel(r'$z$')


    fig.suptitle(r'$\mathcal{S} = \ln|(\nabla^2 \rho) / \rho|$')

    #cbar = fig.colorbar(cax, orientation='vertical')

    fig.set_size_inches(19.2, 10.8)

    fig.tight_layout()

    fig.savefig(f"subch_schlieren.png")


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("plotfiles", type=str, nargs="+",
                   help="list of plotfiles to plot")

    args = p.parse_args()
    print(args)
    plist = find_files(args.plotfiles)

    doit(plist)
