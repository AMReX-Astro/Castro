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
            print(f"{det.name:55} : crashed")
        elif det.has_started:
            print(f"{det.name:55} : t = {det.end_time:8.5f}, # of steps = {det.nsteps:5}, v = {det.v:15.8g} +/- {det.v_sigma:15.8g}")
        else:
            print(f"{det.name:55} : has not started")
