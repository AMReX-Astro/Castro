import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.rcParams.update({'xtick.labelsize': 10,
                     'ytick.labelsize': 10,
                     'font.size': 12})

plt.rc("axes", linewidth=1.5)
plt.rc("lines", markeredgewidth=1.5)

frontier_data = np.loadtxt("frontier-scaling-2024-07-04.txt")

frontier_nodes = frontier_data[:, 0]
frontier_times = frontier_data[:, 3]
frontier_std = frontier_data[:, 4]

frontier_rkc_data = np.loadtxt("frontier-scaling-rkc-2024-07-04.txt")

frontier_rkc_nodes = frontier_rkc_data[:, 0]
frontier_rkc_times = frontier_rkc_data[:, 3]
frontier_rkc_std = frontier_rkc_data[:, 4]

summit_data = np.loadtxt("../summit/scaling_20230407.txt")

summit_nodes = summit_data[:, 0]
summit_times = summit_data[:, 2]
summit_std = summit_data[:, 3]

frontier_bignet_data = np.loadtxt("frontier-scaling-2024-07-04-subch_simple.txt")

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
            ls="None", marker="x", label="Frontier (ROCm 6.0)")
ax.errorbar(frontier_rkc_nodes, frontier_rkc_times, yerr=frontier_rkc_std,
            ls="None", marker="x", label="Frontier (ROCm 6.0; RKC integrator)")
ax.errorbar(summit_nodes, summit_times, yerr=summit_std,
            ls="None", marker="^", label="Summit (CUDA 11.4)")
ax.errorbar(frontier_bignet_nodes, frontier_bignet_times, yerr=frontier_bignet_std,
            ls="None", marker="o", label="Frontier (ROCm 6.0; big network)")

c, t = trend_line(frontier_nodes, frontier_times)
ax.plot(c, t, alpha=0.5, linestyle=":")

c, t = trend_line(frontier_bignet_nodes, frontier_bignet_times)
ax.plot(c, t, alpha=0.5, linestyle=":")


ax.set_ylabel("wallclock time / step")
ax.set_xlabel("number of nodes")

ax.set_xscale("log")
ax.set_yscale("log")

ax.legend()

ax.set_title("3D XRB flame scaling")

fig.savefig("frontier_flame_wave_scaling.png")
