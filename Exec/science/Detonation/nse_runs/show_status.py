#!/bin/env python3

import glob
import tqdm

import detonation as dt

if __name__ == "__main__":

    # get all the data
    run_dirs = glob.glob("det_[rst]*")
    dets = []
    for run in tqdm.tqdm(run_dirs):
        dets.append(dt.Detonation(run))

    for det in sorted(dets):

        if det.crashed:
            print("{:55} : crashed".format(det.name))
        elif det.has_started:
            print("{:55} : t = {:8.5f}, # of steps = {:5}, v = {:15.8g} +/- {:15.8g}".format(det.name, det.end_time, det.nsteps, det.v, det.v_sigma))
        else:
            print("{:55} : has not started".format(det.name))
