#!/usr/bin/env python

# render Temperature for the subchandra problem setup

import sys

import matplotlib
import numpy as np

import yt
from yt.frontends.boxlib.api import CastroDataset
from yt.units import cm
#from yt.visualization.volume_rendering.render_source import VolumeSource
from yt.visualization.volume_rendering.api import Scene, create_volume_source

matplotlib.use('agg')


def doit(plotfile):

    ds = CastroDataset(plotfile)
    ds._periodicity = (True, True, True)

    field = ('boxlib', 'Temp')
    ds._get_field_info(field).take_log = True

    sc = Scene()

    vol = create_volume_source(ds.all_data(), field=field)
    #vol.use_ghost_zones = True

    sc.add_source(vol)


    # transfer function
    _vals = [1.e8, 2.e8, 5.e8, 1.e9, 2.e9, 5.e9]
    vals = []
    for v in _vals:
        vals.append(np.log10(v))

    alpha = [0.05, 0.2, 0.3, 0.3, 0.4, 0.5]

    sigma = 0.05

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cmap = "Oranges"

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
    cam.zoom(3.0)
    sc.camera = cam

    sc.save_annotated("{}_Temp_annotated.png".format(plotfile),
                      label_fontsize="18",
                      sigma_clip=3,
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:6.3f}",
                                      dict(horizontalalignment="left", fontsize="18")],
                                     [(0.5, 0.95),
                                      "Castro simulation of double detonation SN Ia",
                                      dict(color="y", fontsize="24",
                                           horizontalalignment="center")]])


if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)
