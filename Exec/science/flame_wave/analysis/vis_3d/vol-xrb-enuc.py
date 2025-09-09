#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

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

    field = ('boxlib', 'enuc')
    ds._get_field_info(field).take_log = True

    sc = Scene()


    # add a volume: select a sphere
    #center = (0, 0, 0)
    #R = (5.e8, 'cm')

    #dd = ds.sphere(center, R)

    vol = create_volume_source(ds.all_data(), field=field)
    sc.add_source(vol)


    # transfer function
    vals = [16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5]
    sigma = 0.1

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cmap = "viridis"

    for v in vals:
        if v < 19.0:
            alpha = 0.25
        else:
            alpha = 0.75

        tf.sample_colormap(v, sigma**2, alpha=alpha, colormap=cmap)

    sc.get_source(0).transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")
    cam.resolution = (1920, 1280)

    # view 1

    cam.position = [0.75*ds.domain_right_edge[0],
                    0.5*(ds.domain_left_edge[1] + ds.domain_right_edge[1]),
                    ds.domain_right_edge[2]] # + 0.25 * (ds.domain_right_edge[2] - ds.domain_left_edge[2])]

    # look toward the center
    center = 0.5 * (ds.domain_left_edge + ds.domain_right_edge)
    # set the center in the vertical direction to be the height of the underlying base layer
    center[-1] = 2000*cm

    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal, north_vector=[0., 0., 1.])
    cam.set_width(3*ds.domain_width)
    cam.zoom(3.0)
    sc.camera = cam

    sc.save_annotated(f"{plotfile}_Hnuc_annotated_side.png",
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:7.5f} s",
                                      dict(horizontalalignment="left")]])

    # view 2

    dx = ds.domain_right_edge[0] - ds.domain_left_edge[0]
    cam.position = [0.5*(ds.domain_left_edge[0] + ds.domain_right_edge[0]) + 0.0001 * dx,
                    0.5*(ds.domain_left_edge[1] + ds.domain_right_edge[1]),
                    ds.domain_right_edge[2]] # + 0.25 * (ds.domain_right_edge[2] - ds.domain_left_edge[2])]

    # look toward the center
    center = 0.5 * (ds.domain_left_edge + ds.domain_right_edge)
    # set the center in the vertical direction to be the height of the underlying base layer
    center[-1] = 2000*cm

    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal, north_vector=[0., 0., 1.])
    cam.set_width(3*ds.domain_width)
    cam.zoom(0.6)
    sc.camera = cam

    sc.save_annotated(f"{plotfile}_Hnuc_annotated_top.png",
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:7.5f} s",
                                      dict(horizontalalignment="left")]])

if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)
