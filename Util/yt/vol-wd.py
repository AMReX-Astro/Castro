#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import sys

import yt
import numpy as np
from yt.visualization.volume_rendering.api import \
    Scene, \
    Camera, \
    VolumeSource


# this is for the wdconvect problem

def doit(plotfile):

    ds = yt.load(plotfile)
    ds.periodicity = (True, True, True)

    cm = "coolwarm"

    field = ('boxlib', 'density')
    ds._get_field_info(field).take_log = True
        
    sc = Scene()


    # add a volume: select a sphere
    vol = VolumeSource(ds, field=field)
    sc.add_source(vol)


    # transfer function
    vals = [-1, 0, 1, 2, 4, 5, 6, 7]
    #vals = [0.1, 1.0, 10, 100., 1.e4, 1.e5, 1.e6, 1.e7]
    sigma = 0.1

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()
    cm = "coolwarm"
    cm = "spectral"
    for v in vals:
        if v < 4:
            alpha = 0.1
        else:
            alpha = 0.5
        tf.sample_colormap(v, sigma**2, colormap=cm, alpha=alpha)

    sc.get_source(0).transfer_function = tf

        
    cam = Camera(ds, lens_type="perspective")
    cam.resolution = (1280, 720)
    cam.position = 1.5*ds.arr(np.array([0.0, 5.e9, 5.e9]), 'cm')
    
    # look toward the center -- we are dealing with an octant
    center = 0.5*(ds.domain_left_edge + ds.domain_right_edge)
    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal,
                           north_vector=[0., 0., 1.])
    cam.set_width(ds.domain_width)

    sc.camera = cam
    #sc.annotate_axes()
    #sc.annotate_domain(ds)

    pid = plotfile.split("plt")[1]
    sc.render()
    sc.save("wdmerger_{}.png".format(pid), sigma_clip=6.0)
    sc.save_annotated("wdmerger_annotated_{}.png".format(pid), 
                      text_annotate=[[(0.05, 0.05), 
                                      "t = {:.3f}".format(float(ds.current_time.d)),
                                      dict(horizontalalignment="left")],
                                     [(0.5,0.95), 
                                      "Castro simulation of merging white dwarfs",
                                      dict(color="y", fontsize="24",
                                           horizontalalignment="center")]])

if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)


        
