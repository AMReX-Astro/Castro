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

resolution = (1920, 1080)

# this is for the wdconvect problem

def doit(plotfile):

    ds = CastroDataset(plotfile)
    ds._periodicity = (True, True, True)

    field = ('boxlib', 'z_velocity')
    ds._get_field_info(field).take_log = False

    sc = Scene()


    # add a volume: select a sphere
    #center = (0, 0, 0)
    #R = (5.e8, 'cm')

    #dd = ds.sphere(center, R)

    vol = create_volume_source(ds.all_data(), field=field)
    sc.add_source(vol)


    # transfer function
    vals = [-1.e7, -5.e6, -2.5e6, 2.5e6, 5.e6, 1.e7]
    sigma = 5.e5

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cmap = "coolwarm"

    for v in vals:
        tf.sample_colormap(v, sigma**2, alpha=0.2, colormap=cmap)

    sc.get_source(0).transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")
    cam.resolution = resolution

    # view 1

    center = 0.5*(ds.domain_left_edge + ds.domain_right_edge)

    cam.position = [2.5*ds.domain_right_edge[0],
                    2.5*ds.domain_right_edge[1],
                    center[2]+0.5*ds.domain_right_edge[2]]

    # look toward the center -- we are dealing with an octant
    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal,
                           north_vector=[0., 0., 1.])
    cam.set_width(ds.domain_width)
    cam.zoom(1.3)

    sc.camera = cam

    sc.save(f"{plotfile}_zvel.png", sigma_clip=4.0)

    sc.annotate_axes(alpha=0.005, thickness=6)
    sc.annotate_domain(ds, color=np.array([0.05, 0.05, 0.05, 0.05]))
    sc.save(f"{plotfile}_zvel_axes.png", sigma_clip=3.0)

    sc.save_annotated(f"{plotfile}_zvel_annotated.png",
                      sigma_clip=4.0, #label_fmt="%6.2g",
                      label_fontsize="16",
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:6.4f} s",
                                      dict(horizontalalignment="left", fontsize=16)]])


if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)
