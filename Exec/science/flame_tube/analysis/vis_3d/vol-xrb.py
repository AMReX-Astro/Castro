#!/usr/bin/env python
# pylint: disable=wrong-import-position, wrong-import-order, protected-access, use-dict-literal

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


def doit(plotfile):

    ds = yt.load(plotfile)

    plotfile = plotfile.replace("run_", "analysis_")
    Path(plotfile).parent.mkdir(parents=True, exist_ok=True)
    ds._periodicity = (True, True, True)

    # center = ds.domain_center
    # # set the center in the vertical direction to be the height of the underlying base layer
    # center[-1] = 2000*cm

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

    ds._get_field_info(('boxlib', 'abar')).take_log = False
    render_field(ds, ('boxlib', 'abar'), plotfile,
                 ref_points=ref_points,
                 side_alpha=0.005,
                 mid_alpha=0.001,
                 top_alpha=0.001,
                 annotated_kwargs={"label_fmt": "%.2f"},
                 sigma_clip=3.0)

    ds._get_field_info(('boxlib', 'enuc')).take_log = True
    render_field(ds, ('boxlib', 'enuc'), plotfile,
                 ref_points=ref_points,
                 side_alpha=0.005,
                 mid_alpha=0.001,
                 top_alpha=0.0001,
                 sigma_clip=3.0)

    del ds


def make_tf(field):
    if field == ('boxlib', 'abar'):
        # transfer function
        vals = [4.25, 4.5, 4.75, 5, 5.25, 5.5]
        sigma = 0.05

        tf =  yt.ColorTransferFunction((min(vals), max(vals)))

        tf.clear()

        cmap = "plasma_r"

        for v in vals:
            if v < 5.0:
                alpha = 0.25
            else:
                alpha = 0.75

            tf.sample_colormap(v, sigma**2, alpha=alpha, colormap=cmap)

    elif field == ('boxlib', 'enuc'):
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

    return tf


def render_field(ds, field, output_prefix, ref_points, *,
                 side_alpha=0.01, mid_alpha=0.01, top_alpha=0.01,
                 annotated_kwargs=None, **kwargs):

    if annotated_kwargs is None:
        annotated_kwargs = {}

    output_prefix = f"{output_prefix}_{field[1]}"

    sc = Scene()


    # add a volume: select a sphere
    #center = (0, 0, 0)
    #R = (5.e8, 'cm')

    #dd = ds.sphere(center, R)

    vol = create_volume_source(ds.all_data(), field=field)
    sc.add_source(vol)

    vol.transfer_function = make_tf(field)

    cam = sc.add_camera(ds, lens_type="perspective")
    cam.resolution = resolution

    center = ds.domain_center
    # set the center in the vertical direction to be the height of the underlying base layer
    center[-1] = 2000*cm

    # view 1 (side)

    setup_side_camera(cam, ds, center)
    sc.camera = cam
    print(f"side view: {cam!r}", flush=True)

    sc.save(f"{output_prefix}_noaxes_side.png", **kwargs)

    # set the alpha value for the annotations
    ref_points.colors[:, 3] = side_alpha
    sc.add_source(ref_points)
    # sc.annotate_axes(alpha=side_alpha, thickness=6)
    sc.annotate_domain(ds, color=[1, 1, 1, side_alpha])

    sc.save(f"{output_prefix}_side.png", **kwargs)

    sc.save_annotated(f"{output_prefix}_annotated_side.png",
                      **kwargs, **annotated_kwargs,
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:7.5f} s",
                                      dict(horizontalalignment="left")]])

    remove_annotations(sc, vol)

    # view 2 (upper side)

    setup_mid_camera(cam, ds, center)
    sc.camera = cam
    print(f"mid view: {cam!r}", flush=True)

    sc.save(f"{output_prefix}_noaxes_mid.png", **kwargs)

    # set the alpha value for the annotations
    ref_points.colors[:, 3] = mid_alpha
    sc.add_source(ref_points)
    # sc.annotate_axes(alpha=mid_alpha, thickness=6)
    sc.annotate_domain(ds, color=[1, 1, 1, mid_alpha])

    sc.save(f"{output_prefix}_mid.png", **kwargs)

    sc.save_annotated(f"{output_prefix}_annotated_mid.png",
                      **kwargs, **annotated_kwargs,
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:7.5f} s",
                                      dict(horizontalalignment="left")]])

    remove_annotations(sc, vol)

    # view 3 (top)

    setup_top_camera(cam, ds, center)
    sc.camera = cam
    print(f"top view: {cam!r}", flush=True)

    sc.save(f"{output_prefix}_noaxes_top.png", **kwargs)

    # set the alpha value for the annotations
    # ref_points.colors[:, 3] = top_alpha
    # sc.add_source(ref_points)
    # sc.annotate_axes(alpha=top_alpha, thickness=6)
    sc.annotate_domain(ds, color=[1, 1, 1, top_alpha])

    sc.save(f"{output_prefix}_top.png", **kwargs)

    sc.save_annotated(f"{output_prefix}_annotated_top.png",
                      **kwargs, **annotated_kwargs,
                      text_annotate=[[(0.05, 0.05),
                                      f"t = {ds.current_time.d:7.5f} s",
                                      dict(horizontalalignment="left")]])

    del sc, vol


def remove_annotations(sc, vol):
    # print(sc.sources, flush=True)
    for key, source in list(sc.sources.items()):
        if source is not vol:
            sc.sources.pop(key)


def setup_side_camera(cam, ds, center):
    cam.set_position([center[0],
                      ds.domain_left_edge[1] - ds.domain_width[0],
                      center[2]],  # + 0.25 * (ds.domain_right_edge[2] - ds.domain_left_edge[2])]
                     north_vector=[0.0, 0.0, 1.0])
    cam.set_focus(center)

    # width[0] and width[1] are the length and height of the image plane
    # width[2] is the distance from the camera to the image plane
    # x is horizontal, z is vertical
    # put the image plane at the near side of the domain
    width = [ds.domain_width[0], ds.domain_width[2],
             min(abs(cam.position[1] - ds.domain_left_edge[1]),
                 abs(cam.position[1] - ds.domain_right_edge[1]))]
    cam.set_width(width)
    cam.zoom(1 / 1.1)


def setup_mid_camera(cam, ds, center):
    setup_side_camera(cam, ds, center)

    # rotate camera around the x-axis, keeping the distance from the focus constant
    cam.rotate(30 * 180.0 / np.pi, rot_center=center, rot_vector=np.array([1.0, 0.0, 0.0]))


def setup_top_camera(cam, ds, center):
    cam.set_position([ds.domain_center[0],
                      ds.domain_center[1],
                      ds.domain_right_edge[2] + ds.domain_width[0]],
                     north_vector=[0.0, 1.0, 0.0])
    cam.set_focus(center)

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
        sys.exit("ERROR: no plotfiles specified")

    for file in sys.argv[1:]:
        doit(file)
