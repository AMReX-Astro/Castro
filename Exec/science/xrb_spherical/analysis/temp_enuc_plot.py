#!/usr/bin/env python3

import numpy as np
import yt
import argparse
import matplotlib.pyplot as plt
from yt.frontends.boxlib.api import CastroDataset
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed

# Set some fontsize settings
SMALL_SIZE = 18
MEDIUM_SIZE = 24
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)
plt.rc('xtick.major', size=7, width=2)
plt.rc('xtick.minor', size=5, width=1)
plt.rc('ytick.major', size=7, width=2)
plt.rc('ytick.minor', size=5, width=1)

def process_single_file(fname, field_list, weighted_field, cutoff_quantile):
    """Process a single file and return results for all fields.
    This calculates the field_list weighted by weighted_field for values
    above cutoff_quantile"""

    ds = CastroDataset(fname)
    time = ds.current_time.in_units('ms')

    field_results = []
    for f in field_list:
        ad = ds.all_data()
        f_array = ad[f].to_ndarray()

        # Get cutoff value
        cutoff = np.quantile(f_array, cutoff_quantile)

        # Find values above the cutoff value
        f_cutoff = ad.exclude_below(f, cutoff)

        # Get the weighted field
        avg_weighted_field = f_cutoff.quantities.weighted_average_quantity(f, weighted_field)
        field_results.append(avg_weighted_field)

    return time, field_results


def get_weight_fields_concurrent(fnames, field_list=['Temp', 'enuc'],
                                weighted_field='density',
                                output='profile.dat',
                                cutoff_quantile=0.99,
                                max_workers=None):
    """
    using concurrent.futures for parallelization over reading
    different data files. Then save and return the data.
    """

    process_func = partial(process_single_file,
                          field_list=field_list,
                          weighted_field=weighted_field,
                          cutoff_quantile=cutoff_quantile)

    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_fname = {executor.submit(process_func, fname): fname
                          for fname in fnames}

        # Collect results as they complete
        try:
            for future in as_completed(future_to_fname):
                fname = future_to_fname[future]
                try:
                    result = future.result()
                    results.append(result)
                except Exception as exc:
                    print(f'File {fname} generated an exception: {exc}')
                    raise
        except KeyboardInterrupt:
            print(
                "\n*** got ctrl-c, cancelling remaining tasks and waiting for existing ones to finish...\n",
                flush=True,
            )
            executor.shutdown(wait=True, cancel_futures=True)
            sys.exit(1)

    # Get the time and all fields as arrays
    times = [result[0] for result in results]
    all_field_results = [result[1] for result in results]

    avg_weighted_fields = [[] for _ in range(len(field_list))]
    for field_results in all_field_results:
        for i, field_result in enumerate(field_results):
            avg_weighted_fields[i].append(field_result)

    # Sort results by time
    times = np.array(times)
    sort_ind = np.argsort(times)
    avg_weighted_fields = np.array(avg_weighted_fields)
    avg_weighted_fields = avg_weighted_fields[:, sort_ind]
    times = times[sort_ind]

    # Append time to weighted_fields result
    data = np.vstack((avg_weighted_fields, times))

    # Save data and return
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
    parser.add_argument('-j', '--jobs', default=1, type=int,
                        help="""Number of workers to plot in parallel""")
    parser.add_argument('-d', '--data', type=str,
                        help='''If the average density weighted Temp and enuc
                        are already pre-computed, pass in the data file here.
                        Then the calculation is skipped.''')
    # parser.add_argument('-f', '--fields', nargs='+', type=str,
    #                     default=['Temp', 'enuc'],
    #                     help='''Different field names for plotting''')
    # parser.add_argument('-w', '--weighted_field', type=str, default='density',
    #                     help='''Weight different fields by this field.''')

    args = parser.parse_args()

    # Find the data points
    # Right now just find average Temp and enuc weighted by density
    # data = get_weight_fields(args.fnames, field_list=['Temp', 'enuc'],
    #                          weighted_field='density',
    #                          output=args.output,
    #                          cutoff_quantile=args.quantile)

    if args.data is None:
        data = get_weight_fields_concurrent(args.fnames, field_list=['Temp', 'enuc'],
                                            weighted_field='density',
                                            output=args.output,
                                            cutoff_quantile=args.quantile,
                                            max_workers=args.jobs)
    else:
        data = np.loadtxt(args.data, delimiter=',')

    t_array = data[-1, :]

    # Do Plotting
    fig, ax = plt.subplots()
    ax.plot(t_array, data[0, :], color='blue', linewidth=3, label=r'$\left<T\right>$')
    ax.set_xlabel('time [ms]')
    ax.set_ylabel(r'$\left<T \right>$ [GK]')

    ax_twin = ax.twinx()
    ax_twin.plot(t_array, data[1, :] * 1e-19, color='red', linewidth=3, label=r'$\left<\dot{e}_{nuc}\right>$')
    ax_twin.set_ylabel(r'$\left<\dot{e}_{nuc}\right>$ [$10^{19}$erg/g/s]')

    # Combine legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax_twin.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2,
              loc='lower right', frameon=False)

    fig.tight_layout()
    fig.set_size_inches(8, 8)
    fig.savefig("avg_weighted_temp_enuc.png", bbox_inches="tight")
