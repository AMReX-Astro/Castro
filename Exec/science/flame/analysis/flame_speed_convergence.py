#!/usr/bin/env python

import flame_speed as fs

import argparse
import glob
import os

import matplotlib.pyplot as plt








class Flame:

    def __init__(self, name, dx, time, speed, width):
        self.name = name
        self.dx = dx
        self.time = time
        self.speed = speed
        self.width = width

    def __str__(self):
        return "{} @ {} s, dx = {} : {} cm/s , {} cm wide".format(
            self.name, self.time, self.dx, self.speed, self.width)


if __name__ == "__main__":

    # get the directories for the runs we are going to process through
    # the arg list
    p = argparse.ArgumentParser()

    p.add_argument("run_dirs", type=str, nargs="+",
                   help="list of directories containing simulation plotfiles")

    args = p.parse_args()

    # process each run and find its flame speed at the last time
    # interval

    flames = []

    for rdir in args.run_dirs:

        cwd = os.getcwd()
        os.chdir(rdir)

        # find all of the plotfiles -- the timestep number may be 5, 6, or 7 digits,
        # so we need to be careful with sorting so it is done numerically
        files = glob.glob("*plt?????") + glob.glob("*plt??????")
        prefix = files[0].split("plt")[0] + "plt"
        plot_nums = sorted([p.split("plt")[1] for p in files], key=int)
        sorted_files = [f"{prefix}{q}" for q in plot_nums]

        # we'll operate on the last 2 files.  Compute the width and
        # flame_speed from these
        f1 = fs.Profile(sorted_files[-1])
        f2 = fs.Profile(sorted_files[-2])

        x1 = f1.find_x_for_T(T_0=2.e9)
        x2 = f2.find_x_for_T(T_0=2.e9)

        v = (x1 - x2)/(f1.time - f2.time)
        width = 0.5*(f1.find_flame_width() + f2.find_flame_width())

        time = 0.5*(f1.time + f2.time)

        dx = f1.get_dx()

        flames.append(Flame(prefix, dx, time, v, width))

        os.chdir(cwd)

    # summary
    for flame in flames:
        print(flame)

    # compute the errors -- we will take the flame with the smallest dx to be the
    # exact
    v_ref = -1.e30
    dx_ref = 1.e30

    for flame in flames:
        if flame.dx < dx_ref:
            v_ref = flame.speed
            dx_ref = flame.dx

    dxs = []
    errors = []
    for flame in flames:
        if flame.dx == dx_ref:
            continue
        dxs.append(flame.dx)
        errors.append(abs(flame.speed - v_ref)/v_ref)
        print(dxs[-1], errors[-1])

    # make some plots
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(dxs, errors)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.savefig("flame_speed_conv.png")



