import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update({'xtick.labelsize': 9,
                     'ytick.labelsize': 9,
                     'font.size': 11})

plt.rc("axes", linewidth=1.5)
plt.rc("lines", markeredgewidth=1.5)

# we need this for normalization
nsteps = 25

class ScalingRun(object):
    def __init__(self, nodes=1, max_grid=1,
                 time=0.0, std=0.0, coarse_timestep=0.0,
                 react=0.0, parallelcopy=0.0, diffusion=0.0, ctu_hydro=0.0):
        self.nodes = nodes
        self.max_grid = max_grid
        self.time = time
        self.std = std
        self.coarse_timestep = coarse_timestep
        self.react = react
        self.parallelcopy = parallelcopy
        self.diffusion = diffusion
        self.ctu_hydro = ctu_hydro

def trend_line(c, t):
    cnew = np.array(sorted(list(set(c))))
    trend = t[0]*c[0]/cnew[:]
    return cnew, trend

truns = []

data = np.loadtxt("scaling_20200613.txt")

for row in data:
    truns.append(ScalingRun(nodes=row[0], max_grid=row[1],
                           time=row[2], std=row[3], coarse_timestep=row[4],
                           react=row[5], parallelcopy=row[6], diffusion=row[7], ctu_hydro=row[8]))


# for each node count, find the fastest case
runs = []

nodes = (q.nodes for q in truns)

for n in nodes:
    current_runs = [q for q in truns if q.nodes == n]
    fastest = None
    for run in current_runs:
        if fastest is None:
            fastest = run
        else:
            if run.time < fastest.time:
                fastest = run
    runs.append(fastest)

nodes = [q.nodes for q in runs]
times = [q.time  for q in runs]
std = [q.std for q in runs]

react = [q.react/nsteps for q in runs]
hydro = [q.ctu_hydro/nsteps for q in runs]
diffuse = [q.diffusion/nsteps for q in runs]
copy = [q.parallelcopy/nsteps for q in runs]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(nodes, times, label = "timestep", marker="s")

c, t = trend_line(nodes, times)
ax.plot(c, t)

ax.scatter(nodes, hydro, label = "hydro")
ax.scatter(nodes, react, label = "reactions")
ax.scatter(nodes, diffuse, label = "diffusion")
ax.scatter(nodes, copy, label = "MPI")

ax.set_yscale("log")
ax.set_xscale("log")

ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.yaxis.set_minor_formatter(mticker.FormatStrFormatter("%4.1f"))
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())

ax.set_xlim(200, 2500)
#ax.set_ylim(5, 50)

ax.set_ylabel("wallclock time / step")
ax.set_xlabel("number of nodes")

leg = ax.legend(facecolor="#555555", framealpha=0.1)

fig.savefig("flame_wave_scaling.png")

