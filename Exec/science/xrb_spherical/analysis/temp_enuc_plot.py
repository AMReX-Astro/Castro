#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
from yt.frontends.boxlib.api import CastroDataset
import numpy as np
import yt

# Set some fontsize
SMALL_SIZE = 18
MEDIUM_SIZE = 24
BIGGER_SIZE = 28

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

def get_weight_fields(fnames, field_list=['Temp', 'enuc'],
                      weighted_field='density',
                      output='profile.dat',
                      cutoff_quantile=0.99):

    avg_weighted_fields = [[] for _ in range(len(field_list))]
    times = []

    # Loop over different data files, each represent different time
    for fname in fnames:
        ds = CastroDataset(fname)
        times.append(ds.current_time.in_units('ms'))

        # loop over different desired fields
        for i, f in enumerate(field_list):
            ad = ds.all_data()
            f_array = ad[f].to_ndarray()

            # Get the cutoff value corresponding to whatever quantile
            cutoff = np.quantile(f_array, cutoff_quantile)

            # Get yt data container that excludes value below the cutoff
            f_cutoff = ad.exclude_below(f, cutoff)

            # Get the weighted average using the weighted_field
            avg_weighted_field = f_cutoff.quantities.weighted_average_quantity(f, weighted_field)

            avg_weighted_fields[i].append(avg_weighted_field)

    # Sort by time just in case:
    times = np.array(times)
    sort_ind = np.argsort(times)

    avg_weighted_fields = np.array(avg_weighted_fields)
    avg_weighted_fields = avg_weighted_fields[:, sort_ind]

    # Append time array
    data = np.vstack((avg_weighted_fields, times))

    # Save data
    np.savetxt(output, data, delimiter=',')

    return data

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='''This script produces
    a plot showing the density weighted average Temperature and enuc
    from a series of data files''')

    parser.add_argument('fnames', nargs='+', type=str,
                        help='''data file names for generating the plot''')
    parser.add_argument('-q', '--quantile', type=float, default=0.99,
                        help='''Only consider zones that have values above
                        this quantile. Default is 0.99, i.e. only zones
                        above 99% quantile is considered for the average''')
    parser.add_argument('-o', '--output', type=str, default='profile.dat',
                        help='''output name containing relevant data''')
    # parser.add_argument('-f', '--fields', nargs='+', type=str,
    #                     default=['Temp', 'enuc'],
    #                     help='''Different field names for plotting''')
    # parser.add_argument('-w', '--weighted_field', type=str, default='density',
    #                     help='''Weight different fields by this field.''')

    args = parser.parse_args()

    # Find the data points
    # Right now just find average Temp and enuc weighted by density
    data = get_weight_fields(args.fnames, field_list=['Temp', 'enuc'],
                             weighted_field='density',
                             output=args.output,
                             cutoff_quantile=args.quantile)

    t_array = data[-1, :]

    # Do Plotting
    fig, ax = plt.subplots()
    ax.plot(t_array, data[0, :], color='blue', linewidth=3)
    ax.set_xlabel('time [ms]')
    ax.set_ylabel(r'$\left<T \right>$ [GK]')

    ax_twin = ax.twinx()
    ax_twin.plot(t_array, data[1, :] * 1e-19, color='red', linewidth=3)
    ax_twin.set_ylabel(r'$\left<\dot{e}_{nuc}\right>$ [$10^{19}$erg/g/s]')

    # Combine legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax_twin.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, frameon=False)

    fig.tight_layout()
    fig.set_size_inches(8, 8)
    fig.savefig("avg_weighted_temp_enuc.png", bbox_inches="tight")
