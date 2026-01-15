#!/usr/bin/env python

# render a single species for the massive star problem setup
# pass in the plotfile name _and_ the name of the species

import sys

import matplotlib
import numpy as np

import yt
from yt.frontends.boxlib.api import CastroDataset
from yt.units import cm
from yt.visualization.volume_rendering.api import Scene, create_volume_source

matplotlib.use('agg')


def doit(plotfile, species):

    ds = CastroDataset(plotfile)
    ds._periodicity = (True, True, True)

    zoom = 1.5

    t_drive = 0.0
    if "[*] castro.drive_initial_convection_tmax" in ds.parameters:
        t_drive = ds.parameters["[*] castro.drive_initial_convection_tmax"]
    elif "castro.drive_initial_convection_tmax" in ds.parameters:
        t_drive = ds.parameters["castro.drive_initial_convection_tmax"]
    print(t_drive)

    sc = Scene()

    field = f"X({species})"

    vol = create_volume_source(ds.all_data(), field=field)
    sc.add_source(vol)


    # transfer function
    vals = [-2, -1.5, -1, -0.5, 0]
    alpha = [0.01, 0.05, 0.1, 0.2, 0.4]
    sigma = 0.01

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cmap = "viridis"

    for v, a in zip(vals, alpha):
        tf.sample_colormap(v, sigma**2, alpha=a, colormap=cmap)

    sc.get_source(0).transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")
    cam.resolution = (1920, 1280)

    # view 1

    cam.position = [ds.domain_right_edge[0],
                    ds.domain_right_edge[1],
                    ds.domain_right_edge[2]]

    # look toward the center
    center = 0.5 * (ds.domain_left_edge + ds.domain_right_edge)
    # set the center in the vertical direction to be the height of the underlying base layer

    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal, north_vector=[0., 0., 1.])
    cam.set_width(ds.domain_width)
    cam.zoom(zoom)
    sc.camera = cam

    sc.save_annotated(f"{plotfile}_X{species}_zoom{zoom}.png",
                      label_fontsize="18",
                      label_fmt="%.3f",
                      sigma_clip=3,
                      text_annotate=[[(0.05, 0.05),
                                      f"$t - \\tau_\\mathrm{{drive}}$ = {float(ds.current_time) - t_drive:6.1f} s",
                                      dict(horizontalalignment="left", fontsize="18")]])


if __name__ == "__main__":

    # Choose a field
    plotfile = ""
    species = ""

    try:
        plotfile = sys.argv[1]
    except:
        sys.exit("ERROR: no plotfile specified")

    try:
        species = sys.argv[2]
    except:
        sys.exit("ERROR: no species specified")

    doit(plotfile, species)
