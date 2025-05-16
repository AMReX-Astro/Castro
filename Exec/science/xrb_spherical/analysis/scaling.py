#!/usr/bin/env python3


import sys
import re
import numpy as np
import matplotlib.pyplot as plt

'''
This script is intended to
Read in different output files, i.e. xxxxx.out that corresponds
to runs that ran with different number of nodes and started with
the same checkpoint file.

Now within each xxxx.out files, we look for number of MPI processes used,
this is number of CPU or GPU depending on system.

To convert to number of nodes used, divide by the corresponding number
of GPUs per node for specific machine.

***Note: Perlmutter has 4 GPUs and Frontier has 8 GPUs.***
'''

fnames = sys.argv[1:]

MPI_pattern = re.compile(r"MPI initialized with (\d+) MPI processes")
CoarseTimestep_pattern = re.compile(r"Coarse TimeStep time:\s*([0-9]*\.[0-9]+)")

numMPIs = []
averageCoarseTimesteps = []

for fname in fnames:
    coarseTimesteps = []
    with open(fname) as f:
        for line in f:
            mpiLine = MPI_pattern.search(line)
            if mpiLine:
                numMPIs.append(int(mpiLine.group(1)))

            timestepLine = CoarseTimestep_pattern.search(line)
            if timestepLine:
                coarseTimesteps.append(float(timestepLine.group(1)))

    averageCoarseTimestep = sum(coarseTimesteps) / len(coarseTimesteps)
    averageCoarseTimesteps.append(averageCoarseTimestep)

numMPIs = np.array(numMPIs)
averageCoarseTimesteps = np.array(averageCoarseTimesteps)

# Sort by number of MPI used, just in case.
index = np.argsort(numMPIs)
numMPIs = numMPIs[index]
averageCoarseTimesteps = averageCoarseTimesteps[index]

# Speedup = t(0) / t(N)
# ratio between the timestep used via least number of processors vs others.
speedup = averageCoarseTimesteps[0] / averageCoarseTimesteps

# ratio between the least number of processors used vs others
MPIratio = numMPIs / numMPIs[0]

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
