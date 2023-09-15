#!/usr/bin/env python
# pylint: disable=wrong-import-position, wrong-import-order, protected-access
# pylint: disable=too-many-locals, too-many-statements, use-dict-literal

import matplotlib

matplotlib.use('agg')

import sys
from pathlib import Path

import numpy as np
import yt
from yt.units import cm
from yt.visualization.volume_rendering.api import (PointSource, Scene,
                                                   create_volume_source)

resolution = (1920, 1080)

# this is for the wdconvect problem


def doit(plotfile):

    ds = yt.load(plotfile)

    plotfile = plotfile.replace("run_", "analysis_")
    Path(plotfile).parent.mkdir(parents=True, exist_ok=True)
    ds._periodicity = (True, True, True)

    field = ('boxlib', 'enuc')
    ds._get_field_info(field).take_log = True

    sc = Scene()


    # add a volume: select a sphere
    #center = (0, 0, 0)
    #R = (5.e8, 'cm')

    #dd = ds.sphere(center, R)

    center = ds.domain_center
    # set the center in the vertical direction to be the height of the underlying base layer
    center[-1] = 2000*cm

    vol = create_volume_source(ds.all_data(), field=field)
    sc.add_source(vol)

    axes_point = min(ds.domain_center)
    ref_points = PointSource(
        positions=np.array([
            ds.domain_left_edge,
            [axes_point, ds.domain_left_edge[1], ds.domain_left_edge[2]],
            [ds.domain_left_edge[0], axes_point, ds.domain_left_edge[2]],
            [ds.domain_left_edge[0], ds.domain_left_edge[1], axes_point],
            # center,
        ]),
        colors=np.array([
            [1, 1, 1, 0.01],
            [1, 0, 0, 0.01],
            [0, 1, 0, 0.01],
            [0, 0, 1, 0.01],
            # [0, 1, 1, 0.01],
        ]),
        radii=np.array([
            5,
            3,
            3,
            3,
            # 5,
        ])
    )

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

    vol.transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")
    cam.resolution = resolution

    # view 1 (side)

    setup_side_camera(cam, ds, center)
    sc.camera = cam
    print(f"side view: {cam!r}", flush=True)

    sc.save(f"{plotfile}_Hnuc_noaxes_side.png")

    # set the alpha value for the annotations
    alpha = 0.02
    ref_points.colors[:, 3] = alpha
    sc.add_source(ref_points)
    # sc.annotate_axes(alpha=alpha, thickness=6)
    sc.annotate_domain(ds, color=[1, 1, 1, alpha])

    sc.save(f"{plotfile}_Hnuc_side.png")

    sc.save_annotated(f"{plotfile}_Hnuc_annotated_side.png",
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:7.5f} s",
                                      dict(horizontalalignment="left")]])

    # view 2 (top)

    # remove the annotation sources for now
    print(sc.sources, flush=True)
    for key, source in list(sc.sources.items()):
        if source is not vol:
            sc.sources.pop(key)

    setup_top_camera(cam, ds, center)
    sc.camera = cam
    print(f"top view: {cam!r}", flush=True)

    sc.save(f"{plotfile}_Hnuc_noaxes_top.png")

    # set the alpha value for the annotations
    alpha = 0.005
    ref_points.colors[:, 3] = alpha
    sc.add_source(ref_points)
    # sc.annotate_axes(alpha=alpha, thickness=6)
    sc.annotate_domain(ds, color=[1, 1, 1, alpha])

    sc.save(f"{plotfile}_Hnuc_top.png")

    sc.save_annotated(f"{plotfile}_Hnuc_annotated_top.png",
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:7.5f} s",
                                      dict(horizontalalignment="left")]])


def setup_side_camera(cam, ds, center):
    cam.set_position([center[0],
                      ds.domain_left_edge[1] - ds.domain_width[0] / 2,
                      center[2]]) # + 0.25 * (ds.domain_right_edge[2] - ds.domain_left_edge[2])]

    # look toward the center
    normal = center - cam.position
    normal /= np.sqrt(normal.dot(normal))

    cam.set_focus(center)
    cam.switch_orientation(normal_vector=normal, north_vector=[0., 0., 1.])
    # width[0] and width[1] are the length and height of the image plane
    # width[2] is the distance from the camera to the image plane
    # x is horizontal, z is vertical
    # put the image plane at the near side of the domain
    width = [ds.domain_width[0], ds.domain_width[2],
             min(abs(cam.position[1] - ds.domain_left_edge[1]),
                 abs(cam.position[1] - ds.domain_right_edge[1]))]
    cam.set_width(width)
    cam.zoom(1 / 1.1)


def setup_top_camera(cam, ds, center):
    cam.set_position([ds.domain_center[0],
                      ds.domain_center[1],
                      ds.domain_right_edge[2] + ds.domain_width[0] / 2],
                     north_vector=[0., 1., 0.])

    normal = center - cam.position
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal, north_vector=[0., 1., 0.])
    # width[0] and width[1] are the length and height of the image plane
    # width[2] is the distance from the camera to the image plane
    # x is horizontal, y is vertical
    # put the image plane at the top of the domain
    width = [ds.domain_width[0], ds.domain_width[1],
             min(abs(cam.position[2] - ds.domain_left_edge[2]),
                 abs(cam.position[2] - ds.domain_right_edge[2]))]
    cam.set_width(width)
    cam.zoom(1 / 1.1)


if __name__ == "__main__":

    if not sys.argv[1:]:
        sys.exit("ERROR: no plotfile specified")

    for file in sys.argv[1:]:
        doit(file)
