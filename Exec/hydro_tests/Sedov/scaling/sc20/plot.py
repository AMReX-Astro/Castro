import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

def plot(unnormalized, show_linear):

    results_dir = 'scaling_results/'

    results_files = [result for result in os.listdir(results_dir) if 'out' in result]

    n_gpus_per_node = 6

    # We are going to collect three scaling curves. The first is a fiducial case
    # that is a true weak scaling test: work per rank is held fixed. The second is
    # a best case: at each node count we try several different options for max_grid_size
    # and report the best one. This represents the optimal tradeoff between load balancing
    # and per-grid efficiency. The third is a worst case for each node count. It represents
    # the case where the user was not careful in their choice of max_grid_size, and results
    # in them having too many small boxes.

    fom_fiducial = []
    num_fiducial = []

    fom_best = []
    num_best = []

    fom_worst = []
    num_worst = []

    fiducial_grid_size = 64
    fiducial_ncell_size = 256
    fiducial_gpu_counts = [n_gpus_per_node,
                           n_gpus_per_node * 8,
                           n_gpus_per_node * 8**2,
                           n_gpus_per_node * 8**3]

    for results_file in results_files:

        names = results_file.split('.')

        if not 'ngpu' in names:
            continue

        n_gpu = int(names[names.index('ngpu') + 1])
        grid_size = int(names[names.index('grid') + 1])
        ncell_size = int(names[names.index('ncell') + 1])

        if n_gpu not in num_best:
            num_best.append(n_gpu)
            fom_best.append(0.0)

        if n_gpu not in num_worst:
            num_worst.append(n_gpu)
            fom_worst.append(0.0)

        for line in open(results_dir + results_file):

            if len(line.split()) > 0 and line.split()[0] == 'Average' and line.split()[6] == 'microsecond:':

                temp_fom = float(line.split()[7])

                if n_gpu in fiducial_gpu_counts and grid_size == fiducial_grid_size:

                    if n_gpu not in num_fiducial:
                        num_fiducial.append(n_gpu)
                        fom_fiducial.append(0.0)

                    idx = num_fiducial.index(n_gpu)
                    if n_gpu == fiducial_gpu_counts[0]:
                        scale_factor = 1
                    elif n_gpu == fiducial_gpu_counts[1]:
                        scale_factor = 2
                    elif n_gpu == fiducial_gpu_counts[2]:
                        scale_factor = 4
                    elif n_gpu == fiducial_gpu_counts[3]:
                        scale_factor = 8

                    if ncell_size == fiducial_ncell_size * scale_factor:
                        fom_fiducial[idx] = temp_fom

                idx = num_best.index(n_gpu)

                if fom_best[idx] < temp_fom or fom_best[idx] == 0.0:
                    fom_best[idx] = temp_fom

                idx = num_worst.index(n_gpu)

                if fom_worst[idx] > temp_fom or fom_worst[idx] == 0.0:
                    fom_worst[idx] = temp_fom

    fom_fiducial = np.array(fom_fiducial)
    num_fiducial = np.array(num_fiducial) / n_gpus_per_node

    fom_fiducial = np.array([x for _, x in sorted(zip(num_fiducial, fom_fiducial))])
    num_fiducial = sorted(num_fiducial)

    fom_linear = [fom_fiducial[0] * num for num in num_fiducial]
    num_linear = num_fiducial

    fom_best = np.array(fom_best)
    num_best = np.array(num_best) / n_gpus_per_node

    fom_best = np.array([x for _, x in sorted(zip(num_best, fom_best))])
    num_best = sorted(num_best)

    fom_worst = np.array(fom_worst)
    num_worst = np.array(num_worst) / n_gpus_per_node

    fom_worst = np.array([x for _, x in sorted(zip(num_worst, fom_worst))])
    num_worst = sorted(num_worst)

    if not unnormalized:
        fom_best = fom_best / fom_fiducial[0] / num_best
        fom_worst = fom_worst / fom_fiducial[0] / num_worst
        fom_linear = fom_linear / fom_fiducial[0] / num_linear
        fom_fiducial = fom_fiducial / fom_fiducial[0] / num_fiducial

    plt.plot(num_fiducial, fom_fiducial, marker='x', markersize=14, linestyle='-', lw=4, zorder=100, label='Weak scaling')
    if show_linear:
        plt.plot(num_linear, fom_linear, linestyle='-', color='k', lw=2, label='Linear scaling')
    plt.plot(num_best, fom_best, marker='o', markersize=12, linestyle='--', lw=4, label='Best case')
    plt.plot(num_worst, fom_worst, marker='s', markersize=12, linestyle='--', lw=4, label='Worst case')

    plt.xlim([0.9 * min(num_best), 1.1 * max(num_best)])
    if not unnormalized:
        plt.ylim([0, 1.25 * max(fom_best)])
    else:
        plt.ylim([0, 1.10 * max(fom_best)])

    if not unnormalized:
        plt.ylabel('Throughput (normalized)', fontsize=20)
    else:
        plt.ylabel('Throughput', fontsize=20)
    plt.xlabel('Number of nodes', fontsize=20)
    plt.title('Weak scaling of CASTRO Sedov', fontsize=20)
    plt.xscale('log', basex=2)
    ax = plt.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks([1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
    plt.tick_params(labelsize=14)
    plt.legend(loc='best', numpoints=3, markerscale=0.6, handlelength=5)
    plt.tight_layout()

    plt.savefig('scaling_results/scaling.eps')
    plt.savefig('scaling_results/scaling.png')

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--unnormalized', action='store_true',
                        help='Should we not normalize the figure of merit per GPU?')
    parser.add_argument('--show_linear', action='store_true',
                        help='Should we show a linear scaling curve?')

    args = parser.parse_args()

    plot(args.unnormalized, args.show_linear)

if __name__ == "__main__":

    main()
