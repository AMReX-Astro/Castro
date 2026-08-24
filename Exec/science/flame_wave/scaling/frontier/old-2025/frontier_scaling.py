import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.rcParams.update({'xtick.labelsize': 10,
                     'ytick.labelsize': 10,
                     'font.size': 12})

plt.rc("axes", linewidth=1.5)
plt.rc("lines", markeredgewidth=1.5)

frontier_data = np.genfromtxt("frontier-scaling-iso7-20250409.txt")

frontier_nodes = frontier_data[:, 0]
frontier_times = frontier_data[:, 3]
frontier_std = frontier_data[:, 4]


frontier_bignet_data = np.genfromtxt("frontier-scaling-he-burn-22a-20250409.txt")

frontier_bignet_nodes = frontier_bignet_data[:, 0]
frontier_bignet_times = frontier_bignet_data[:, 3]
frontier_bignet_std = frontier_bignet_data[:, 4]


def trend_line(c, t):
    cnew = np.array(sorted(list(set(c))))
    cnew = np.linspace(cnew.min(), cnew.max(), 256, endpoint=True)
    trend = t[0]*c[0]/cnew[:]
    return cnew, trend


# first by nodes

fig, ax = plt.subplots(1)

ax.errorbar(frontier_nodes, frontier_times, yerr=frontier_std,
            ls="None", marker="x", label="Frontier (ROCm 6.3.1)")
ax.errorbar(frontier_bignet_nodes, frontier_bignet_times, yerr=frontier_bignet_std,
            ls="None", marker="o", label="Frontier (ROCm 6.3.1; big network)")

c, t = trend_line(frontier_nodes, frontier_times)
ax.plot(c, t, alpha=0.5, linestyle=":", color="k")

c, t = trend_line(frontier_bignet_nodes, frontier_bignet_times)
ax.plot(c, t, alpha=0.5, linestyle=":", color="k")


ax.set_ylabel("wallclock time / step")
ax.set_xlabel("number of nodes")

ax.set_xscale("log")
ax.set_yscale("log")

ax.legend(fontsize="small")

ax.set_title("3D XRB flame scaling")

fig.tight_layout()
fig.savefig("frontier_flame_wave_scaling.png")
