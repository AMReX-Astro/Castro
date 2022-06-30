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

    field = ('boxlib', 'enuc_symlog')
    ds._get_field_info(field).take_log = False

    sc = Scene()


    # add a volume: select a sphere
    #center = (0, 0, 0)
    #R = (5.e8, 'cm')

    #dd = ds.sphere(center, R)

    vol = create_volume_source(ds.all_data(), field=field)
    sc.add_source(vol)


    # transfer function
    vals = [-22.0, -21.0, -20.0, -19.0, -18.0, 18.0, 19.0, 20.0, 21.0, 22.0]
    alpha = [0.5,  0.4,   0.3,   0.2,   0.2,   0.2,  0.2,  0.3,  0.4,  0.5]
    sigma = 0.1

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()

    cmap = "coolwarm"

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
    #cam.zoom(3.0)
    sc.camera = cam

    sc.save_annotated("{}_enuc_annotated.png".format(plotfile),
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
