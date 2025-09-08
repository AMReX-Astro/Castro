#!/usr/bin/env python3


import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import yt
import sys
import numpy as np

class Profile:
    def __init__(self, x, T, rho, e, p, u):
        self.x = x
        self.T = T
        self.rho = rho
        self.e = e
        self.p = p
        self.u = u

def get_profile(plotfile):

    ds = yt.load(plotfile)

    time = float(ds.current_time)
    ad = ds.all_data()

    # Sort the ray values by 'x' so there are no discontinuities
    # in the line plot
    srt = np.argsort(ad['x'])

    x_coord = np.array(ad['x'][srt])
    rho = np.array(ad['density'][srt])
    temp = np.array(ad['Temp'][srt])
    e = np.array(ad['rho_e'][srt])/rho
    p = np.array(ad['pressure'][srt])
    u = np.array(ad['xmom'][srt])/rho

    return Profile(x_coord, temp, rho, e, p, u)


def doit(ppm_file, sdc_file, exact_file, outfile):

    # this gets the header information correct
    have_exact = False
    if exact_file != "":
        exact = np.genfromtxt(exact_file, skip_header=2, names=True)
        have_exact = True

    f = plt.figure()
    f.set_size_inches(7.0, 10.0)

    ppm_profile = get_profile(ppm_file)
    sdc_profile = get_profile(sdc_file)

    ax_rho = f.add_subplot(411)
    ax_rho.plot(ppm_profile.x, ppm_profile.rho, color="C0", marker="x", label="ppm")
    ax_rho.plot(sdc_profile.x, sdc_profile.rho, color="C1", marker="x", label="sdc")
    if have_exact:
        ax_rho.plot(exact["x"], exact["rho"], color="0.5", label="exact")

    ax_rho.legend(frameon=False)
    ax_rho.set_ylabel(r"$\rho~(\mathrm{g~cm^{-3}})$")
    ax_rho.set_yscale("log")

    ax_u = f.add_subplot(412)
    ax_u.plot(ppm_profile.x, ppm_profile.u, color="C0", marker="x")
    ax_u.plot(sdc_profile.x, sdc_profile.u, color="C1", marker="x")
    ax_u.set_ylabel(r"$u~(\mathrm{cm~s^{-1}})$")
    if have_exact:
        ax_u.plot(exact["x"], exact["u"], color="0.5", label="exact")
    #ax_u.set_yscale("log")
    ax_u.set_ylim(bottom=1.e-10)

    ax_p = f.add_subplot(413)
    ax_p.plot(ppm_profile.x, ppm_profile.p, color="C0", marker="x")
    ax_p.plot(sdc_profile.x, sdc_profile.p, color="C1", marker="x")
    if have_exact:
        ax_p.plot(exact["x"], exact["p"], color="0.5", label="exact")

    ax_p.set_ylabel(r"$p~(\mathrm{erg~cm^{-3}})$")
    ax_p.set_yscale("log")

    ax_T = f.add_subplot(414)
    ax_T.plot(ppm_profile.x, ppm_profile.T, color="C0", marker="x")
    ax_T.plot(sdc_profile.x, sdc_profile.T, color="C1", marker="x")
    ax_T.set_ylabel(r"$T~(\mathrm{K})$")
    ax_T.set_yscale("log")

    ax_T.set_xlabel(r"$x$ (cm)")

    f.savefig(outfile)


if __name__ == "__main__":

    p = argparse.ArgumentParser()

    p.add_argument("--ppm", type=str, nargs="+",
                   help="base name for PPM files")
    p.add_argument("--sdc", type=str, nargs="+",
                   help="base name for SDC files")
    p.add_argument("--exact", type=str, default="",
                   help="exact solution (ASCI columns)")
    p.add_argument("-o", type=str, default="plot.png",
                   help="output plot name")

    args = p.parse_args()

    ppm_prefix = args.ppm[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.ppm], key=int)
    last_ppm_file = f"{ppm_prefix}{plot_nums[-1]}"

    sdc_prefix = args.sdc[0].split("plt")[0] + "plt"
    plot_nums = sorted([p.split("plt")[1] for p in args.sdc], key=int)
    last_sdc_file = f"{sdc_prefix}{plot_nums[-1]}"

    doit(last_ppm_file, last_sdc_file, args.exact, args.o)

