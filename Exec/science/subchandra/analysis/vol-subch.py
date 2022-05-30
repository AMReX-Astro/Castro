#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import numpy as np

import sys

import yt
from yt.frontends.boxlib.api import CastroDataset
import numpy as np
#from yt.visualization.volume_rendering.render_source import VolumeSource
from yt.visualization.volume_rendering.api import create_volume_source, Scene
from yt.units import cm

# this is for the wdconvect problem

def doit(plotfile):

    ds = CastroDataset(plotfile)
    ds._periodicity = (True, True, True)

    field = ('boxlib', 'Temp')
    ds._get_field_info(field).take_log = True

    sc = Scene()


    # add a volume: select a sphere
    #center = (0, 0, 0)
    #R = (5.e8, 'cm')

    #dd = ds.sphere(center, R)

    vol = create_volume_source(ds.all_data(), field=field)
    sc.add_source(vol)


    # transfer function
    vals_tmp = [5.e7, 1.e8, 5.e8, 1.e9, 2.e9, 2.5e9, 3.e9]
    vals = []
    for v in vals_tmp:
        vals.append(np.log10(v))

    sigma = 0.05

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cmap = "viridis"

    for v in vals:
        if v < 9.01:
            alpha = 0.25
        else:
            alpha = 0.75

        tf.sample_colormap(v, sigma**2, alpha=alpha, colormap=cmap)

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
    #cam.zoom(3.0)
    sc.camera = cam

    sc.save_annotated("{}_Hnuc_annotated.png".format(plotfile),
                      text_annotate=[[(0.05, 0.05),
                                      "t = {}".format(ds.current_time.d),
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
