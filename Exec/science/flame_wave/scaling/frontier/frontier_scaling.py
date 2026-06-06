import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.rcParams.update({'xtick.labelsize': 10,
                     'ytick.labelsize': 10,
                     'font.size': 12})

plt.rc("axes", linewidth=1.5)
plt.rc("lines", markeredgewidth=1.5)

_data = np.genfromtxt("frontier-scaling-vode-iso7-20260606.txt")
frontier_iso7_vode_nodes = _data[:, 0]
frontier_iso7_vode_times = _data[:, 1]
frontier_iso7_vode_std = _data[:, 2]

_data = np.genfromtxt("frontier-scaling-vode-ase-20260606.txt")
frontier_ase_vode_nodes = _data[:, 0]
frontier_ase_vode_times = _data[:, 1]
frontier_ase_vode_std = _data[:, 2]

_data = np.genfromtxt("frontier-scaling-rosenbrock-iso7-20260606.txt")
frontier_iso7_rosenbrock_nodes = _data[:, 0]
frontier_iso7_rosenbrock_times = _data[:, 1]
frontier_iso7_rosenbrock_std = _data[:, 2]

_data = np.genfromtxt("frontier-scaling-rosenbrock-ase-20260606.txt")
frontier_ase_rosenbrock_nodes = _data[:, 0]
frontier_ase_rosenbrock_times = _data[:, 1]
frontier_ase_rosenbrock_std = _data[:, 2]


def trend_line(c, t):
    cnew = np.array(sorted(list(set(c))))
    cnew = np.linspace(cnew.min(), cnew.max(), 256, endpoint=True)
    trend = t[0]*c[0]/cnew[:]
    return cnew, trend


# first by nodes

fig, ax = plt.subplots(1)

ax.errorbar(frontier_iso7_vode_nodes, frontier_iso7_vode_times, yerr=frontier_iso7_vode_std,
            ls="None", marker="x", label="iso7 / VODE")

ax.errorbar(frontier_ase_vode_nodes, frontier_ase_vode_times, yerr=frontier_ase_vode_std,
            ls="None", marker="x", label="ase / VODE")

ax.errorbar(frontier_iso7_rosenbrock_nodes, frontier_iso7_rosenbrock_times, yerr=frontier_iso7_rosenbrock_std,
            ls="None", marker="x", label="iso7 / Rosenbrock")

ax.errorbar(frontier_ase_rosenbrock_nodes, frontier_ase_rosenbrock_times, yerr=frontier_ase_rosenbrock_std,
            ls="None", marker="x", label="ase / Rosenbrock")

c, t = trend_line(frontier_iso7_vode_nodes, frontier_iso7_vode_times)
ax.plot(c, t, alpha=0.5, linestyle=":", color="k")

c, t = trend_line(frontier_ase_vode_nodes, frontier_ase_vode_times)
ax.plot(c, t, alpha=0.5, linestyle=":", color="k")


ax.set_ylabel("wallclock time / step")
ax.set_xlabel("number of nodes")

ax.set_xscale("log")
ax.set_yscale("log")

ax.legend(fontsize="small")

ax.set_title("3D XRB flame scaling")

fig.tight_layout()
fig.savefig("frontier_flame_wave_scaling.png")
