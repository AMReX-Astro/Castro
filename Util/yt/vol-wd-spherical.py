#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import sys

import yt
import numpy as np
from yt.visualization.volume_rendering.api import \
    Scene, \
    VolumeSource


# this is for the wdconvect problem

def doit(plotfile):

    ds = yt.load(plotfile)
    ds.periodicity = (True, True, True)

    field = ('boxlib', 'density')
    ds._get_field_info(field).take_log = True

    sc = Scene()


    # add a volume: select a sphere
    vol = VolumeSource(ds, field=field)
    vol.use_ghost_zones = True

    sc.add_source(vol)


    # transfer function
    vals = [-1, 0, 1, 2, 3, 4, 5, 6, 7]
    #vals = [0.1, 1.0, 10, 100., 1.e4, 1.e5, 1.e6, 1.e7]
    sigma = 0.1

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()
    cm = "coolwarm"
    cm = "spectral"
    for v in vals:
        if v < 3:
            alpha = 0.1
        else:
            alpha = 0.5
        tf.sample_colormap(v, sigma**2, colormap=cm, alpha=alpha)

    sc.get_source(0).transfer_function = tf


    # for spherical, youtube recommends an "equirectangular" aspect ratio
    # (2:1), suggested resolution of 8192 x 4096
    # see: https://support.google.com/youtube/answer/6178631?hl=en
    #
    # also see: http://yt-project.org/docs/dev/cookbook/complex_plots.html#various-lens-types-for-volume-rendering
    # the 2:1 is 2*pi in phi and pi in theta
    cam = sc.add_camera(ds, lens_type="spherical")
    #cam.resolution = (8192, 4096)
    cam.resolution = (4096, 2048)

    # look toward the +x initially
    cam.focus = ds.arr(np.array([ds.domain_left_edge[0], 0.0, 0.0]), 'cm')

    # center of the domain -- eventually we might want to do the
    # center of mass
    cam.position = ds.arr(np.array([0.0, 0.0, 0.0]), 'cm')

    # define up
    cam.north_vector = np.array([0., 0., 1.])

    normal = (cam.focus - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal,
                           north_vector=[0., 0., 1.])

    # there is no such thing as a camera width -- the entire volume is rendered
    #cam.set_width(ds.domain_width)

    #sc.annotate_axes()
    #sc.annotate_domain(ds)

    pid = plotfile.split("plt")[1]
    sc.render()
    sc.save("wdmerger_{}_spherical.png".format(pid), sigma_clip=6.0)

if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)
