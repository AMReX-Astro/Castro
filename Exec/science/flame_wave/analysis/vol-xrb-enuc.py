#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import sys

import yt
from yt.frontends.boxlib.api import CastroDataset
import numpy as np
#from yt.visualization.volume_rendering.render_source import VolumeSource
from yt.visualization.volume_rendering.api import create_volume_source, Scene


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
    vals = [17.5, 18.0, 18.5, 19.0, 19.5]
    sigma = 0.1

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cm = "viridis"

    for v in vals:
        if v < 19.0:
            alpha = 0.25
        else:
            alpha = 0.75

        tf.sample_colormap(v, sigma**2, alpha=alpha, colormap=cm)

    sc.get_source(0).transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")        
    cam.resolution = (1280, 720)
    cam.position = 1.0*ds.domain_right_edge
    
    # look toward the center -- we are dealing with an octant
    center = 0.5 * (ds.domain_left_edge + ds.domain_right_edge)
    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal,
                           north_vector=[0., 0., 1.])
    cam.set_width(3*ds.domain_width)
    cam.zoom(2.0)
    sc.camera = cam
    #sc.annotate_axes(alpha=0.05)
    #sc.annotate_domain(ds, color=np.array([0.05, 0.05, 0.05, 0.05]))
    #sc.annotate_grids(ds, alpha=0.05)

    sc.render()
    clip = 1.5
    sc.save("{}_Hnuc".format(plotfile), sigma_clip=clip)
    sc.save_annotated("{}_Hnuc_annotated.png".format(plotfile), 
                      sigma_clip=clip,
                      text_annotate=[[(0.05, 0.05), 
                                      "t = {}".format(ds.current_time.d),
                                      dict(horizontalalignment="left")],
                                     [(0.5,0.95), 
                                      "Castro simulation of XRB flame spreading",
                                      dict(color="y", fontsize="24",
                                           horizontalalignment="center")]])

if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)


        
