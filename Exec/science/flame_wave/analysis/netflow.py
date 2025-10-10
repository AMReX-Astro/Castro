import argparse
import importlib
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pynucastro as pyna
import yt
from pynucastro.yt_utils import get_point, to_conditions
from yt.units import cm

# simulation

def doit(plotfile, reference):

    # the reference file is used to find the location at which we are
    # going to examine the network flow

    ds_r = yt.load(reference)

    val, loc = ds_r.find_max("enuc")

    print("using location for network flow: ", loc)

    # now load the plotfile we will plot

    if plotfile == reference:
        ds = ds_r
    else:
        ds = yt.load(plotfile)

    field = "enuc"

    # get the thermodynamic data at our desired location
    # from the current (not reference) plotfile

    pt = get_point(ds, loc)
    rho, T, comp = to_conditions(pt)


    # network
    mp = os.getenv("MICROPHYSICS_HOME")

    moddir = Path(mp) / "networks" / "he-burn" / "cno-he-burn-34am"
    sys.path.insert(0, str(moddir))

    import cno_he_burn_34am as pynet

    net = pynet.create_network()
    net.summary()

    # setup the figure where we will host the image

    fig = plt.figure(figsize=(19.2, 10.8))

    inch = 0.3
    W, H = fig.get_size_inches()
    margin = (inch / W, inch / H)
    fig.set_layout_engine('constrained',
                          rect=[margin[0], margin[1],
                                1.0 - 2.0*margin[0], 1.0 - 2.0*margin[1]])

    gs = fig.add_gridspec(nrows=2, ncols=1,
                          height_ratios=[1.5, 5.0])

    ax_img = fig.add_subplot(gs[0, 0])
    gs_flow = gs[1, 0].subgridspec(2, 2, height_ratios=(20, 1))


    # setup domain size
    # we'll use a buffer size that is proportional to the number of zones
    # in the data we are plotting

    xmin = ds.domain_left_edge[0]
    xmax = ds.domain_right_edge[0]
    xctr = 0.5*(xmin + xmax)
    L_x = xmax - xmin

    ymin = 0.0*cm
    ymax = 1.2288e4*cm

    yctr = 0.5*(ymin + ymax)
    L_y = ymax - ymin

    nx, ny, nz = ds.domain_dimensions
    ref = int(np.prod(ds.ref_factors[0:ds.max_level]))

    dx = ds.domain_width[0] / nx / ref
    dy = ds.domain_width[1] / ny / ref

    thinning = 1

    pxls = int(L_x / dx / thinning), int(L_y / dy / thinning)

    # do the plotting

    p = yt.SlicePlot(ds, "theta", field,
                     center=[xctr, yctr, 0.0*cm], width=[L_x, L_y, 0.0*cm])
    p.set_buff_size(pxls)
    p.set_zlim("enuc", 1.e15, 3.e19)
    p.set_background_color("enuc", None)
    p.set_axes_unit("cm")


    # get the image data and replot it on our axes

    im = p[field].image

    # Recreate the image on our axes with same appearance
    data   = im.get_array()
    cmap   = im.get_cmap()
    norm   = im.norm
    extent = im.get_extent()
    origin = getattr(im, "origin", "lower")


    im_new = ax_img.imshow(data, cmap=cmap, norm=norm, extent=extent, origin=origin)
    ax_img.set_xlabel("r [cm]")
    ax_img.set_ylabel("z [cm]")

    # keep physical aspect but within the GridSpec cell
    ax_img.set_aspect("equal", adjustable="box")

    ax_img.scatter([loc[0]], [loc[1]], marker="x", color="red")

    # and add a colorbar
    cb = fig.colorbar(im_new, ax=ax_img, orientation="horizontal", shrink=0.8)
    cb.set_label(field)


    # network flow part

    net.plot(rho, T, comp,
             screen_func=pyna.screening.screen.screen5,
             rotated=True, hide_xp=True, hide_xalpha=True,
             curved_edges=True,
             color_nodes_by_abundance=True,
             node_size=500, node_font_size="9",
             grid_spec=gs_flow)

    # text

    fig.text(0.05, 0.02, f"time = {float(ds.current_time*1000):5.1f} ms",
             transform=fig.transFigure, fontsize="12")

    # finalize
    fig.savefig(f"{os.path.basename(plotfile)}_netflow.png")


if __name__ == "__main__":

    p = argparse.ArgumentParser(description="plot enuc and the flow through the reaction network")

    p.add_argument("plotfile",
                   type=lambda p: Path(p).resolve() if Path(p).is_dir() else parser.error(f"{p} is not a directory"),
                   nargs=1,
                   help="plotfile to visualize (enuc slice + network flow")
    p.add_argument("--reference", type=str,
                   default=None, help="reference file to use for finding location of max enuc")

    args = p.parse_args()
    ref = args.reference
    if args.reference is None:
        ref = args.plotfile[0]

    print("reference = ", ref)
    doit(args.plotfile[0], ref)
