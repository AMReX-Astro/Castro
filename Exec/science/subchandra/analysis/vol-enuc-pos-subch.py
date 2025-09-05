#!/usr/bin/env python

import sys

import matplotlib
import numpy as np

import yt
from yt.frontends.boxlib.api import CastroDataset
from yt.units import cm
#from yt.visualization.volume_rendering.render_source import VolumeSource
from yt.visualization.volume_rendering.api import Scene, create_volume_source

matplotlib.use('agg')




# this is for the wdconvect problem

def _enuc_symlog(field, data):
    f = np.log10(np.abs(data["boxlib", "enuc"]))
    f[f < 10] = 0.0
    return np.copysign(f, data["boxlib", "enuc"])

yt.add_field(
    name = ("boxlib", "enuc_symlog"),
    function = _enuc_symlog,
    sampling_type = "local",
    units = None
)


def doit(plotfile):

    ds = CastroDataset(plotfile)
    ds._periodicity = (True, True, True)

    sc = Scene()

    # enuc

    field = ('boxlib', 'enuc')
    ds._get_field_info(field).take_log = True

    vol = create_volume_source(ds.all_data(), field=field)
    #vol.use_ghost_zones = True

    vals = [18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0]
    sigma = 0.1

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cmap = "viridis"

    for v in vals:
        if v > 21.0:
            alpha = 0.6
        elif v > 19.0:
            alpha = 0.4
        else:
            alpha = 0.2
        tf.sample_colormap(v, sigma**2, alpha=alpha, colormap=cmap)

    vol.set_transfer_function(tf)

    sc.add_source(vol)

    # temp

    field = ('boxlib', 'Temp')
    ds._get_field_info(field).take_log = True

    vol2 = create_volume_source(ds.all_data(), field=field)
    vol2.use_ghost_zones = True

    vals = [np.log10(2.e8)]
    sigma = 0.05

    Tmin = 1.e7
    Tmax = 1.e9

    tf =  yt.ColorTransferFunction((np.log10(Tmin), np.log10(Tmax)))

    tf.clear()

    cmap = "Oranges"

    for v in vals:
        alpha = 0.2
        tf.sample_colormap(v, sigma**2, alpha=alpha, colormap=cmap)

    vol2.set_transfer_function(tf)

    sc.add_source(vol2)


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
    cam.set_width(0.5*ds.domain_width)
    cam.zoom(3.0)
    sc.camera = cam

    sc.save_annotated(f"{plotfile}_enuc_pos_annotated.png",
                      sigma_clip=5.0,
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d}",
                                      dict(horizontalalignment="left")],
                                     [(0.5,0.95),
                                      "Castro simulation of double detonation SN Ia",
                                      dict(color="y", fontsize="24",
                                           horizontalalignment="center")]])


if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)
