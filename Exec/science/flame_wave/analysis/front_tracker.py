#!/usr/bin/env python3

import sys
import os
import yt
from yt.units import cm
import numpy as np

tf = lambda file: yt.load(file.rstrip('/'))
ts = list(map(tf, sys.argv[1:]))

output_dir = "front_tracking"

if not os.path.isdir(output_dir):
    
    os.mkdir(output_dir)
    
class Metrics:
    
    metric_names = ["one_percent", "one_tenth_percent"]
    
    global_max = None
    
    def __init__(self, frb):
        
        self.dim = frb.buff_size
        self.xlim = frb.bounds[:2]
        self.ylim = frb.bounds[2:]
        
        self.avg_burns = frb["enuc"].mean(axis=0)
        self.max_avg_burn = self.avg_burns.max()
        self.r = frb["r"].mean(axis=0)
        
        self._ten_percent = None
        self._five_percent = None
        self._one_percent = None
        self._one_tenth_percent = None
        
    @property
    def ten_percent(self):
        
        if self._ten_percent is None:
            
            ind = (self.avg_burns < 1e-1 * self.global_max).argmax()
            self._ten_percent = self.r[ind]
        
        return self._ten_percent
        
    @property
    def five_percent(self):
        
        if self._five_percent is None:
            
            ind = (self.avg_burns < 5e-2 * self.global_max).argmax()
            self._five_percent = self.r[ind]
        
        return self._five_percent
        
    @property
    def one_percent(self):
        
        if self._one_percent is None:
            
            ind = (self.avg_burns < 1e-2 * self.global_max).argmax()
            self._one_percent = self.r[ind]
        
        return self._one_percent
        
    @property
    def one_tenth_percent(self):
        
        if self._one_tenth_percent is None:
            
            ind = (self.avg_burns < 1e-3 * self.global_max).argmax()
            self._one_tenth_percent = self.r[ind]
        
        return self._one_tenth_percent
        
    @property
    def all(self):
        
        return map(lambda name: getattr(self, name), self.metric_names)
    
def get_window_parameters(ds, axis, width=None, center='c'):
    
    width = ds.coordinates.sanitize_width(axis, width, None)
    center, display_center = ds.coordinates.sanitize_center(center, axis)
    xax = ds.coordinates.x_axis[axis]
    yax = ds.coordinates.y_axis[axis]
    bounds = (display_center[xax]-width[0] / 2,
              display_center[xax]+width[0] / 2,
              display_center[yax]-width[1] / 2,
              display_center[yax]+width[1] / 2)
    return bounds, center, display_center
    
def get_width(ds, xlim=None, ylim=None, zlim=None):
    """ Get the width of the plot. """

    if xlim is None: xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
    else: xlim = xlim[0] * cm, xlim[1] * cm

    if ylim is None: ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
    else: ylim = ylim[0] * cm, ylim[1] * cm

    xwidth = (xlim[1] - xlim[0]).in_cgs()
    ywidth = (ylim[1] - ylim[0]).in_cgs()

    if ds.domain_dimensions[2] == 1:
        zwidth = 0.0
    else:
        if zlim is None: zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
        else: zlim = zlim[0] * cm, zlim[1] * cm

        zwidth = (zlim[1] - zlim[0]).in_cgs()

    return xwidth, ywidth, zwidth
    
def get_center(ds, xlim=None, ylim=None, zlim=None):
    """ Get the coordinates of the center of the plot. """

    if xlim is None: xlim = ds.domain_left_edge[0], ds.domain_right_edge[0]
    else: xlim = xlim[0] * cm, xlim[1] * cm

    if ylim is None: ylim = ds.domain_left_edge[1], ds.domain_right_edge[1]
    else: ylim = ylim[0] * cm, ylim[1] * cm

    xctr = 0.5 * (xlim[0] + xlim[1])
    yctr = 0.5 * (ylim[0] + ylim[1])

    if ds.domain_dimensions[2] == 1:
        zctr = 0.0
    else:
        if zlim is None: zlim = ds.domain_left_edge[2], ds.domain_right_edge[2]
        else: zlim = zlim[0] * cm, zlim[1] * cm

        zctr = 0.5 * (zlim[2] + zlim[2])

    return xctr, yctr, zctr
    
def make_metrics(ds):
    
    axis = "theta"
    origin = "native"
    fields = ["enuc", "r"]
    
    xlim = 0.0, 6e4
    ylim = 0.0, 1.5e4
    width = get_width(ds, xlim, ylim)
    center = get_center(ds, xlim, ylim)
    
    dim = 256, 1024
    
    axis = ds.coordinates.axis_id.get(axis, axis)
    bounds, center, display_center = get_window_parameters(ds, axis, width, center)
    slc = ds.slice(axis, center[axis], center=center)
    slc.get_data(fields)
    
    frb = yt.FixedResolutionBuffer(slc, bounds, dim)
    
    return Metrics(frb)
    
def set_max(metrics):
    
    Metrics.global_max = max(map(lambda obj: obj.max_avg_burn, metrics))
    
times = np.array([ds.current_time for ds in ts])
metrics = list(map(make_metrics, ts))
set_max(metrics)
radii = np.array([tuple(obj.all) for obj in metrics])

if len(radii) > 2:
    
    veloc = ((radii[1:] - radii[:-1]).T / (times[1:] - times[:-1])).T
    
else:
    
    veloc = None
    
with open("{}/radii.dat".format(output_dir), 'w') as file:
    
    print(*Metrics.metric_names, file=file)
    
    for seq in radii:
        
        print(*seq, file=file)

if veloc is None: sys.exit("Single dataset, no velocities to output.")

with open("{}/velocities.dat".format(output_dir), 'w') as file:
    
    print(*Metrics.metric_names, file=file)
    
    for seq in veloc:
        
        print(*seq, file=file)

with open("{}/times.dat".format(output_dir), 'w') as file:

    print("time", file=file)

    for time in times:

        print(time, file=file)

print("Task completed.")
