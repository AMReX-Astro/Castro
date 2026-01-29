import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update({'xtick.labelsize': 10,
                     'ytick.labelsize': 10,
                     'font.size': 12})

plt.rc("axes", linewidth=1.5)
plt.rc("lines", markeredgewidth=1.5)

# we need this for normalization
nsteps = 25

class ScalingRun:
    def __init__(self, nodes=1, max_grid=1,
                 time=0.0, std=0.0):
        self.nodes = nodes
        self.max_grid = max_grid
        self.time = time
        self.std = std

def filter_fastest_runs(data):

    truns = []

    for row in data:
        truns.append(ScalingRun(nodes=row[0], max_grid=row[1],
                                time=row[2], std=row[3]))

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

    return runs

def trend_line(c, t):
    cnew = np.array(sorted(list(set(c))))
    cnew = np.linspace(cnew.min(), cnew.max(), 256, endpoint=True)
    trend = t[0]*c[0]/cnew[:]
    return cnew, trend

files = [("scaling_20200613.txt", "2020"),
         ("scaling_20210606.txt", "2021"),
         ("scaling_20220430.txt", "2022")]

# show the change through years

fig, ax = plt.subplots(1)

markers = {"2020": "s", "2021": "x", "2022": "o"}

for f, year in files:

    data = np.loadtxt(f)

    runs = filter_fastest_runs(data)

    nodes = [q.nodes for q in runs]
    times = [q.time  for q in runs]
    std = [q.std for q in runs]

    ax.scatter(nodes, times, label = year, marker=markers[year])

    if year == "2020":
        c, t = trend_line(nodes, times)
        ax.plot(c, t, alpha=0.5, linestyle=":")

ax.set_ylabel("wallclock time / step")
ax.set_xlabel("number of nodes")

leg = ax.legend(facecolor="#555555", framealpha=0.1)

fig.savefig("flame_wave_scaling_history.png")

ax.set_xscale("log")
ax.set_yscale("log")

fig.savefig("flame_wave_scaling_history_log.png")
