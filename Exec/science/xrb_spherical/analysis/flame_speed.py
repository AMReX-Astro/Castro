#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from uncertainties import unumpy

# Set some fontsize
SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('xtick.major', size=7, width=2)
plt.rc('xtick.minor', size=5, width=1)
plt.rc('ytick.major', size=7, width=2)
plt.rc('ytick.minor', size=5, width=1)

# Use tanh + quadratic fit
def tanh_func(t, a0, v0, x0, a, b, c):
    return 0.5*a0*t**2 + v0*t + x0 + a*np.tanh(t/b + c)

# Another version for error propagation
def utanh_func(t, a0, v0, x0, a, b, c):
    return 0.5*a0*t**2 + v0*t + x0 + a*unumpy.tanh(t/b + c)

# def angular_velocity(t, a0, v0, x0, a, b, c):
#     return a0*t + v0 + a / (b * unumpy.cosh(t/b + c)**2)

def quadratic_func(t, a0, v0, x0):
    return 0.5*a0*t**2 + v0*t + x0

def angular_velocity(t, a0, v0, x0):
    return a0*t + v0

def fit_front(times, thetas, tmin, init_guess):
    # Helper function to fit the front

    # Now apply tmin so that we ignore the transient phase during fitting
    cond = times >= tmin
    fit_times = times[cond]
    fit_thetas = thetas[cond]

    # Do the fitting
    # error is given by the square root of the diagonal of the covariance matrix.
    popt, pcov = curve_fit(quadratic_func, fit_times, fit_thetas, p0=init_guess, method="lm")
    err = np.sqrt(np.diag(pcov))

    # Given the fitted parameters, recreate fitted curve along with error propagation
    # Error propagation is handled by the uncertainties package.
    # create fitted params with error
    fitted_params = unumpy.uarray(popt, err)
    theta_fit = quadratic_func(fit_times, *fitted_params)
    theta_nominal = unumpy.nominal_values(theta_fit)
    theta_err = unumpy.std_devs(theta_fit)

    # Get angular velocity in rad / s
    w_fit = angular_velocity(fit_times, *fitted_params) * 1e3
    w_nominal = unumpy.nominal_values(w_fit)
    w_err = unumpy.std_devs(w_fit)

    return {"fit_times": fit_times,
            "theta_nominal": theta_nominal, "theta_err": theta_err,
            "w_nominal": w_nominal, "w_err": w_err}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''This script uses output
    from front_tracker.py to plot the time evolution of flame front θ.''')

    parser.add_argument('tracking_fnames', nargs='+', type=str,
                        help='cvs file generated from front_tracker.py to track flame front position')
    parser.add_argument('--ash-tmin', default=0.0, type=float,
                        help='minimum time for ash front curve fitting. Note that this will not affect plotting')
    parser.add_argument('--flame-tmin', default=0.0, type=float,
                        help='minimum time for ash front curve fitting. Note that this will not affect plotting')
    parser.add_argument('--tmax', type=float,
                        help='maximum time for both plotting and curve fitting. This applies to both flame and ash')
    parser.add_argument('--initial-guess', nargs=3, type=float,
                        default=None, metavar="a0, v0, x0",
                        help="initial guess for the quadratic front fit")
    parser.add_argument('--plot-stride', type=int, default=1,
                        help="""Interval at which we plot the raw front position data.
                        Increasing it can make plot look nicer""")
    parser.add_argument("-o", "--output", default=None, type=str, metavar="FILENAME",
                        help="Output filename (PNG). If not set, shows interactive plot.")

    args = parser.parse_args()
    all_data = []

    # Loop over all tracking files and append them
    for fname in args.tracking_fnames:
        df = pd.read_csv(fname)
        all_data.append(df)

    # concatenate all files together
    # then sort by the time
    tracking_data = pd.concat(all_data, ignore_index=True)
    tracking_data = tracking_data.sort_values(by='time[ms]')

    # Columns available are
    # fname,time[ms],flame_theta,theta_max_avg,max_avg_enuc,ash_theta
    # Get time and theta, these should already be time-sorted
    times = tracking_data['time[ms]']
    flame_thetas = tracking_data['flame_theta']
    ash_thetas = tracking_data['ash_theta']
    t_nuc = tracking_data['t_nuc']
    D = tracking_data['D']
    ocean_height = tracking_data['ocean_height']
    coriolis_param = tracking_data['coriolis_param']
    ash_velocity = tracking_data['ash_velocity']

    # Only plot up to tmax
    if args.tmax is not None:
        cond = times < args.tmax
        times = times[cond]
        flame_thetas = flame_thetas[cond]
        ash_thetas = ash_thetas[cond]
        t_nuc = t_nuc[cond]
        D = D[cond]
        ocean_height = ocean_height[cond]
        coriolis_param = coriolis_param[cond]
        ash_velocity = ash_velocity[cond]

    # Now since t_nuc and diffusion coefficient, D, will be None when dataset is smallplt
    # since it doesn't have enough information. Do filtering on these dataset.
    none_mask = D.notna()
    times_vel = times[none_mask]
    t_nuc = t_nuc[none_mask]
    D = D[none_mask]
    ocean_height = ocean_height[none_mask]
    coriolis_param = coriolis_param[none_mask]
    ash_velocity = ash_velocity[none_mask]

    # Now get the fitted front data
    # init_guess = np.array([0.0, 0.0004, 0.0, 0.01, 5, -2.0])
    init_guess = args.initial_guess

    flame_fit = fit_front(times, flame_thetas, args.flame_tmin, init_guess)
    ash_fit = fit_front(times, ash_thetas, args.ash_tmin, init_guess)

    # Now do plotting
    fig, (ax_theta, ax_velocity) = plt.subplots(2, 1, figsize=(8, 8),
                                                sharex=True, constrained_layout=True)

    # Raw data, downsample the data so markers show up nicely
    stride = args.plot_stride
    ax_theta.plot(times[::stride], flame_thetas[::stride], '*', color='k', markersize=7, label='flame data')
    ax_theta.plot(times[::stride], ash_thetas[::stride], '^', color='k', markersize=7, label='ash data')

    # Fit data
    ax_theta.plot(flame_fit["fit_times"], flame_fit["theta_nominal"], linewidth=4,
                  color='tab:blue', linestyle='--', label='flame fit')

    ax_theta.plot(ash_fit["fit_times"], ash_fit["theta_nominal"], linewidth=4,
                  color='tab:green', linestyle='--', label='ash fit')

    # plot error of the fit
    ax_theta.fill_between(flame_fit["fit_times"],
                          flame_fit["theta_nominal"] - flame_fit["theta_err"],
                          flame_fit["theta_nominal"] + flame_fit["theta_err"],
                          alpha=0.5, color='tab:blue')

    ax_theta.fill_between(ash_fit["fit_times"],
                          ash_fit["theta_nominal"] - ash_fit["theta_err"],
                          ash_fit["theta_nominal"] + ash_fit["theta_err"],
                          alpha=0.5, color='tab:green')

    ax_theta.set_ylabel(r"$\theta$ [rad]")
    ax_theta.set_ylim(0.06, None)
    ax_theta.grid(linestyle=":")
    ax_theta.tick_params(top=True, bottom=True, left=True, right=True)
    ax_theta.legend(frameon=False)
    ax_theta.tick_params(axis="both",direction="in")

    # Then do velocity plotting.
    # Show fitting data and velocity from theoretical prediction.

    # Compute theoretical speed.
    # Conduction Speed (Landau)
    v1 = np.sqrt(D / t_nuc) * 1e-5 # convert to km/s

    # Ageotrophic Speed (Spitkovski 2002)
    g = 1.5e14
    L_R = np.sqrt(g * ocean_height) / coriolis_param
    v2 = L_R / t_nuc * 1e-5

    # Conduction + Ageotrophic (Cavecchi 2013)
    v3 = 2.5 * np.sqrt(D / t_nuc) * L_R / ocean_height * 1e-5

    # Plot theoretical speeds
    ax_velocity.plot(times_vel[::stride], v1[::stride], 'v', color='k', markersize=7, label=r'$\sqrt{D/t_n}$')
    ax_velocity.plot(times_vel[::stride], v2[::stride], 'p', color='k', markersize=7, label=r'$L_R/t_n$')
    ax_velocity.plot(times_vel[::stride], v3[::stride], 'X', color='k', markersize=7, label=r'$L_R/H \sqrt{D/t_n} $')

    # Assume neutron star of radius 11 km, so linear speed is R*omega
    R = 11
    ax_velocity.plot(flame_fit["fit_times"], R*flame_fit["w_nominal"], linewidth=3,
              color='tab:red', linestyle='-.', label='flame speed')

    ax_velocity.plot(ash_fit["fit_times"], R*ash_fit["w_nominal"], linewidth=3,
              color='tab:orange', linestyle='-.', label='ash speed')

    ax_velocity.fill_between(flame_fit["fit_times"],
                             R*(flame_fit["w_nominal"] - flame_fit["w_err"]),
                             R*(flame_fit["w_nominal"] + flame_fit["w_err"]),
                             alpha=0.5, color='tab:red')

    ax_velocity.fill_between(ash_fit["fit_times"],
                             R*(ash_fit["w_nominal"] - ash_fit["w_err"]),
                             R*(ash_fit["w_nominal"] + ash_fit["w_err"]),
                             alpha=0.5, color='tab:orange')

    ax_velocity.set_ylabel(r"R $\omega$ [km $s^{-1}$]")
    ax_velocity.grid(linestyle=":")
    ax_velocity.tick_params(top=True, bottom=True, left=True, right=True)
    ax_velocity.legend(frameon=False)
    ax_velocity.tick_params(axis="both",direction="in")

    ax_velocity.set_xlabel("time [ms]")

    fig.tight_layout()
    fig.set_size_inches(7, 9)
    # Store to output, otherwise show plot
    if args.output is not None:
        fig.savefig(args.output, format="png", bbox_inches="tight")
    else:
        plt.show()
