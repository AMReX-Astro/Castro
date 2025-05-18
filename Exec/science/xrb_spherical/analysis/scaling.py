#!/usr/bin/env python3

import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

def process_outfiles(fnames):
    ''' Given a list of fnames, compute the average coarse timestep
    and number of MPI processors used.
    '''

    MPI_pattern = re.compile(r"MPI initialized with (\d+) MPI processes")
    CoarseTimestep_pattern = re.compile(r"Coarse TimeStep time:\s*([0-9]*\.[0-9]+)")

    numMPIs = []
    averageCoarseTimesteps = []

    for fname in fnames:
        coarseTimesteps = []
        with open(fname) as f:
            for line in f:

                # Find number of MPI
                mpiLine = MPI_pattern.search(line)
                if mpiLine:
                    numMPIs.append(int(mpiLine.group(1)))

                # Find coarse timestep
                timestepLine = CoarseTimestep_pattern.search(line)
                if timestepLine:
                    coarseTimesteps.append(float(timestepLine.group(1)))

        # Find the average
        averageCoarseTimestep = sum(coarseTimesteps) / len(coarseTimesteps)
        averageCoarseTimesteps.append(averageCoarseTimestep)

    numMPIs = np.array(numMPIs)
    averageCoarseTimesteps = np.array(averageCoarseTimesteps)

    # Sort
    index = np.argsort(numMPIs)
    numMPIs = numMPIs[index]
    averageCoarseTimesteps = averageCoarseTimesteps[index]

    return numMPIs, averageCoarseTimesteps


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This script is intended to
    Read in different output files, i.e. xxxxx.out that corresponds
    to runs that ran with different number of nodes and started with
    the same checkpoint file.

    Now within each xxxx.out files, we look for number of MPI processes used,
    this is number of CPU or GPU depending on system.

    To convert to number of nodes used, divide by the corresponding number
    of GPUs per node for specific machine.

    ***Note: Perlmutter has 4 GPUs and Frontier has 8 GPUs.***
    """)

    parser.add_argument("outFiles", nargs='+', type=str,
                        help="""xxxx.out files each representing runs used with
                        different MPI processors""")
    parser.add_argument("-r", "--reference", type=str,
                        help="""A xxxx.out file used as reference for the scaling test.
                        If no reference is provided, then the output file that used
                        the least number of processors will be used as the reference.""")

    args = parser.parse_args()

    numMPIs, averageCoarseTimesteps = process_outfiles(args.outFiles)
    if args.reference is not None:
        numMPIs_ref[0], averageCoarseTimesteps_ref[0] = process_outfiles(args.reference)
    else:
        # If no reference is provided, use the outfile that used the least number
        # of MPI processors. It's index 0 since its already sorted.
        numMPIs_ref = numMPIs[0]
        averageCoarseTimesteps_ref = averageCoarseTimesteps[0]

    # Speedup = t(0) / t(N)
    # ratio between the timestep used via least number of processors vs others.
    speedup = averageCoarseTimesteps_ref / averageCoarseTimesteps

    # ratio between the least number of processors used vs others
    MPIratio = numMPIs / numMPIs_ref

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(MPIratio, speedup, marker='x', color="k", label="simulation")
    ax.plot(MPIratio, MPIratio, linestyle="-.", label="theory")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xlabel("MPI Processor Ratio")
    ax.set_ylabel("Speedup")
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    fig.savefig("xrb_scaling.png")
