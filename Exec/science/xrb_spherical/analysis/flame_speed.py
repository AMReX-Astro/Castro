#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from uncertainties import unumpy

# Set some fontsize
SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22

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

parser = argparse.ArgumentParser(description='''This script uses output
from front_tracker.py to plot the time evolution of flame front θ.
''')

parser.add_argument('tracking_fname', type=str,
                    help='cvs file generated from front_tracker.py to track flame front position')
parser.add_argument('--tmin', default=0.0, type=float,
                    help='minimum time for curve fitting. Note that this will not affect plotting')
parser.add_argument('--tmax', type=float,
                    help='maximum time for both plotting and curve fitting')

args = parser.parse_args()

# Read in data
# data has columns: fname, time, front_theta, theta_max_avg, max_avg, theta_max, max_val.
tracking_data = pd.read_csv(args.tracking_fname)

# Get time and theta, these should already be time-sorted
times = tracking_data['time[ms]']
front_thetas = tracking_data['front_theta']

# Only plot up to tmax
if args.tmax is not None:
    cond = times < args.tmax

times = times[cond]
front_thetas = front_thetas[cond]

# Now do a curve fit to the data.
# Now apply tmin so that we ignore the transient phase during fitting
fit_times = times[args.tmin <= times]
fit_front_thetas = front_thetas[args.tmin <= times]

# Use tanh + quadratic fit
def tanh_func(t, a0, v0, x0, a, b, c):
    return 0.5*a0*t**2 + v0*t + x0 + a*np.tanh(t/b + c)

# Another version for error propagation
def utanh_func(t, a0, v0, x0, a, b, c):
    return 0.5*a0*t**2 + v0*t + x0 + a*unumpy.tanh(t/b + c)

# Give initial guess and solve for different parameters
# error is given by the square root of the diagonal of the covariance matrix.
init_guess = np.array([0.0, 0.0004, 0.0, 0.01, 5, -2.0])
popt, pcov = curve_fit(tanh_func, fit_times, fit_front_thetas, p0=init_guess, method="lm")
err = np.sqrt(np.diag(pcov))

# Given the fitted parameters, recreate fitted curve along with error propagation
# Error propagation is handled by the uncertainties package.
# create fitted params with error
fitted_params = unumpy.uarray(popt, err)
theta_fit = utanh_func(fit_times, *fitted_params)
theta_nominal = unumpy.nominal_values(theta_fit)
theta_err = unumpy.std_devs(theta_fit)

# Now use the fitted parameter to calculate angular velocity
# This is the derivative of utanh_func
def angular_velocity(t, a0, v0, x0, a, b, c):
    return a0*t + v0 + a / (b * unumpy.cosh(t/b + c)**2)

# Get angular velocity in rad / s
w_fit = angular_velocity(fit_times, *fitted_params) * 1e3
w_nominal = unumpy.nominal_values(w_fit)
w_err = unumpy.std_devs(w_fit)

# Now do plotting
fig, ax = plt.subplots()
ax.plot(times, front_thetas, 'x', color='k', label='θ: data')
ax.plot(fit_times, theta_nominal, linewidth=3, linestyle='--', label='θ: fit')
ax.fill_between(fit_times, theta_nominal - theta_err, theta_nominal + theta_err,
                alpha=0.5, color='b', label='θ: 1σ band')
ax.set_xlabel("time [ms]")
ax.set_ylabel("θ [rad]")
ax.set_ylim(0.06, None)

# Create twin ax to plot angular velocity
ax_twin = ax.twinx()
ax_twin.plot(fit_times, w_nominal, linewidth=3, color='r', linestyle='-.', label='ω fit')
ax_twin.fill_between(fit_times, w_nominal - w_err, w_nominal + w_err,
                     alpha=0.5, color='purple', label='ω 1σ band')
ax_twin.set_ylabel("ω [rad/s]")

# Combine legend
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax_twin.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

fig.tight_layout()
fig.set_size_inches(8, 8)
fig.savefig("flame_position.png", bbox_inches="tight")
