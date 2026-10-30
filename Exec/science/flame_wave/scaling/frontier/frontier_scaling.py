import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

plt.rcParams.update({'xtick.labelsize': 10,
                     'ytick.labelsize': 10,
                     'font.size': 12})

plt.rc("axes", linewidth=1.5)
plt.rc("lines", markeredgewidth=1.5)

# iso7
_data = np.genfromtxt("frontier-scaling-vode-iso7-20260606.txt")
frontier_iso7_vode_nodes = _data[:, 0]
frontier_iso7_vode_times = _data[:, 1]
frontier_iso7_vode_std = _data[:, 2]

_data = np.genfromtxt("frontier-scaling-rosenbrock-iso7-20260606.txt")
frontier_iso7_rosenbrock_nodes = _data[:, 0]
frontier_iso7_rosenbrock_times = _data[:, 1]
frontier_iso7_rosenbrock_std = _data[:, 2]

_data = np.genfromtxt("frontier-scaling-rkc-iso7-20260606.txt")
frontier_iso7_rkc_nodes = _data[:, 0]
frontier_iso7_rkc_times = _data[:, 1]
frontier_iso7_rkc_std = _data[:, 2]

# ase
_data = np.genfromtxt("frontier-scaling-vode-ase-20260606.txt")
frontier_ase_vode_nodes = _data[:, 0]
frontier_ase_vode_times = _data[:, 1]
frontier_ase_vode_std = _data[:, 2]

_data = np.genfromtxt("frontier-scaling-rosenbrock-ase-20260606.txt")
frontier_ase_rosenbrock_nodes = _data[:, 0]
frontier_ase_rosenbrock_times = _data[:, 1]
frontier_ase_rosenbrock_std = _data[:, 2]

_data = np.genfromtxt("frontier-scaling-rkc-ase-20260606.txt")
frontier_ase_rkc_nodes = _data[:, 0]
frontier_ase_rkc_times = _data[:, 1]
frontier_ase_rkc_std = _data[:, 2]

# H/He
_data = np.genfromtxt("frontier-h-he-scaling-rosenbrock-cno-he-burn-34am.txt")
frontier_h_he_rosenbrock_nodes = _data[:, 0]
frontier_h_he_rosenbrock_times = _data[:, 1]
frontier_h_he_rosenbrock_std = _data[:, 2]

_data = np.genfromtxt("frontier-h-he-scaling-vode-cno-he-burn-34am.txt")
frontier_h_he_vode_nodes = _data[:, 0]
frontier_h_he_vode_times = _data[:, 1]
frontier_h_he_vode_std = _data[:, 2]

_data = np.genfromtxt("frontier-h-he-scaling-vode32bitjac-cno-he-burn-34am.txt")
frontier_h_he_vode32bitjac_nodes = _data[:, 0]
frontier_h_he_vode32bitjac_times = _data[:, 1]
frontier_h_he_vode32bitjac_std = _data[:, 2]


def trend_line(c, t):
    cnew = np.array(sorted(list(set(c))))
    cnew = np.linspace(cnew.min(), cnew.max(), 256, endpoint=True)
    trend = t[0]*c[0]/cnew[:]
    return cnew, trend


# first by nodes

fig, ax = plt.subplots(1)

plot_iso7 = False

# iso7
if plot_iso7:
    ax.errorbar(frontier_iso7_vode_nodes, frontier_iso7_vode_times, yerr=frontier_iso7_vode_std,
                ls="None", marker="^", color="C0", label="He / iso7 + VODE")

    ax.errorbar(frontier_iso7_rosenbrock_nodes, frontier_iso7_rosenbrock_times, yerr=frontier_iso7_rosenbrock_std,
                ls="None", marker="^", color="C1", label="He / iso7 + Rosenbrock")

    ax.errorbar(frontier_iso7_rkc_nodes, frontier_iso7_rkc_times, yerr=frontier_iso7_rkc_std,
                ls="None", marker="^", color="C2", label="He / iso7 + RKC")

    c, t = trend_line(frontier_iso7_vode_nodes, frontier_iso7_vode_times)
    ax.plot(c, t, alpha=0.5, linestyle=":", color="k")

# ase
ax.errorbar(frontier_ase_vode_nodes, frontier_ase_vode_times, yerr=frontier_ase_vode_std,
            ls="None", marker="o", color="C0", label="He / VODE")

ax.errorbar(frontier_ase_rosenbrock_nodes, frontier_ase_rosenbrock_times, yerr=frontier_ase_rosenbrock_std,
            ls="None", marker="o", color="C1", label="He / Rosenbrock")

ax.errorbar(frontier_ase_rkc_nodes, frontier_ase_rkc_times, yerr=frontier_ase_rkc_std,
            ls="None", marker="o", color="C2", label="He / RKC")

# h_he
ax.errorbar(frontier_h_he_vode_nodes, frontier_h_he_vode_times, yerr=frontier_h_he_vode_std,
            ls="None", marker="x", color="C0", label="H-He / VODE")

ax.errorbar(frontier_h_he_rosenbrock_nodes, frontier_h_he_rosenbrock_times, yerr=frontier_h_he_rosenbrock_std,
            ls="None", marker="x", color="C1", label="H-He / Rosenbrock")

ax.errorbar(frontier_h_he_vode32bitjac_nodes, frontier_h_he_vode32bitjac_times, yerr=frontier_h_he_vode32bitjac_std,
            ls="None", marker="x", color="C3", label="H-He / VODE (32-bit Jac)")


c, t = trend_line(frontier_ase_vode_nodes, frontier_ase_vode_times)
ax.plot(c, t, alpha=0.5, linestyle=":", color="k")

c, t = trend_line(frontier_h_he_vode_nodes, frontier_h_he_vode_times)
ax.plot(c, t, alpha=0.5, linestyle=":", color="k")


ax.set_ylabel("wallclock time / step")
ax.set_xlabel("number of nodes")

ax.set_xscale("log")
ax.set_yscale("log")

ax.legend(fontsize="9", ncol=2)

ax.set_title("3D XRB flame scaling")

ax.grid(ls=":", which="both")

fig.tight_layout()
fig.savefig("frontier_flame_wave_scaling.pdf")
