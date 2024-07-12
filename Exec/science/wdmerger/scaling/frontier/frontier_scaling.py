import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.rcParams.update({'xtick.labelsize': 10,
                     'ytick.labelsize': 10,
                     'font.size': 12})

plt.rc("axes", linewidth=1.5)
plt.rc("lines", markeredgewidth=1.5)

frontier_256base_data = np.loadtxt("frontier_256base_20240709.txt")
frontier_256base_nodes = frontier_256base_data[:, 0]
frontier_256base_times = frontier_256base_data[:, 4]

frontier_512base_data = np.loadtxt("frontier_512base_20240709.txt")
frontier_512base_nodes = frontier_512base_data[:, 0]
frontier_512base_times = frontier_512base_data[:, 4]

frontier_1024base_data = np.loadtxt("frontier_1024base_20240709.txt")
frontier_1024base_nodes = frontier_1024base_data[:, 0]
frontier_1024base_times = frontier_1024base_data[:, 4]


def trend_line(c, t):
    cnew = np.array(sorted(list(set(c))))
    cnew = np.linspace(cnew.min(), cnew.max(), 256, endpoint=True)
    trend = t[0]*c[0]/cnew[:]
    return cnew, trend


# first by nodes

fig, ax = plt.subplots(1)

ax.plot(frontier_256base_nodes, frontier_256base_times, ls="None", marker="x", label="$256^3$ coarse grid")
ax.plot(frontier_512base_nodes, frontier_512base_times, ls="None", marker="x", label="$512^3$ coarse grid")
ax.plot(frontier_1024base_nodes, frontier_1024base_times, ls="None", marker="x", label="$1024^3$ coarse grid")

c, t = trend_line(frontier_256base_nodes, frontier_256base_times)
ax.plot(c, t, alpha=0.5, linestyle=":")

ax.set_ylabel("wallclock time / step")
ax.set_xlabel("number of nodes")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_title("3D wedmerger scaling (Frontier)")

ax.legend()

fig.tight_layout()
fig.savefig("frontier_wdmerger.png")
