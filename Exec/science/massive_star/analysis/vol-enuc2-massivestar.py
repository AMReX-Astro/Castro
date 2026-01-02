#!/usr/bin/env python

# render a enuc using 2 separate transfer functions to better
# show positive and negative

import sys

import matplotlib
import numpy as np

import yt
from yt.frontends.boxlib.api import CastroDataset
from yt.units import cm
from yt.visualization.volume_rendering.api import Scene, create_volume_source

matplotlib.use('agg')

def _enuc_symlog(field, data):
    f = np.log10(np.abs(data["boxlib", "enuc"]))
    f[f < 10] = 0.0
    return np.copysign(f, data["boxlib", "enuc"])

yt.add_field(
    name = ("boxlib", "enuc_symlog"),
    display_name = "", #r"\mathrm{log}_{10}(\epsilon_\mathrm{nuc})~ [\mathrm{erg/g/s}]",
    function = _enuc_symlog,
    sampling_type = "local",
    units = None
)


def doit(plotfile):

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

    field = ("boxlib", "enuc_symlog")
    ds._get_field_info(field).take_log = False

    vals = [12, 13, 14, 15, 15.5, 16, 16.5, 17]
    sigma = 0.15
    alpha = [0.025, 0.05, 0.1, 0.125, 0.15, 0.2, 0.3, 0.4]

    # negative values of enuc

    vol2 = create_volume_source(ds.all_data(), field=field)
    #vol2.use_ghost_zones = True

    vals_r = [-v for v in reversed(vals)]
    alpha_r = [a for a in reversed(alpha)]

    cmap = "Blues_r"

    tf =  yt.ColorTransferFunction((min(vals_r), max(vals_r)))
    tf.clear()

    for v, a in zip(vals_r, alpha_r):
        tf.sample_colormap(v, sigma**2, alpha=a, colormap=cmap)

    vol2.set_transfer_function(tf)

    sc.add_source(vol2)

    # positive values of enuc

    vol = create_volume_source(ds.all_data(), field=field)
    #vol.use_ghost_zones = True

    cmap = "Reds"

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))
    tf.clear()

    for v, a in zip(vals, alpha):
        tf.sample_colormap(v, sigma**2, alpha=a, colormap=cmap)

    vol.set_transfer_function(tf)

    sc.add_source(vol)


    # camera

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

    sc.save_annotated(f"{plotfile}_enuc_zoom{zoom}.png",
                      label_fontsize="18",
                      label_fmt="%.3f",
                      sigma_clip=3,
                      text_annotate=[[(0.05, 0.05),
                                      f"$t - \\tau_\\mathrm{{drive}}$ = {float(ds.current_time) - t_drive:6.1f} s",
                                      dict(horizontalalignment="left", fontsize="18")]])


if __name__ == "__main__":

    # Choose a field
    plotfile = ""

    try:
        plotfile = sys.argv[1]
    except:
        sys.exit("ERROR: no plotfile specified")

    doit(plotfile)
