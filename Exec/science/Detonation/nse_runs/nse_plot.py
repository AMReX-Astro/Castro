import numpy as np
import matplotlib.pyplot as plt

strang = "strang/det_x_plt00741.slice"
strang_cfl025 = "strang_cfl025/det_x_plt01774.slice"
strang_cfl0083 = "strang_cfl0083/det_x_plt05262.slice"

sdc = "SDC/det_x_plt00634.slice"
sdc_cfl025 = "SDC_cfl025/det_x_plt01768.slice"
sdc_cfl0083 = "SDC_cfl0083/det_x_plt05272.slice"

sdc3 = "SDC3/det_x_plt00606.slice"
sdc4 = "SDC4/det_x_plt00607.slice"

datasets = {}

datasets["Strang"] = np.loadtxt(strang)
datasets["Strang (CFL = 0.25)"] = np.loadtxt(strang_cfl025)
datasets["Strang (CFL = 0.083)"] = np.loadtxt(strang_cfl0083)

datasets["SDC"] = np.loadtxt(sdc)
datasets["SDC (CFL = 0.25)"] = np.loadtxt(sdc_cfl025)
datasets["SDC (CFL = 0.083)"] = np.loadtxt(sdc_cfl0083)

datasets["SDC (3 iters)"] =  np.loadtxt(sdc3)
datasets["SDC (4 iters)"] = np.loadtxt(sdc4)

plots = []
plots.append([("Strang", "Strang (CFL = 0.25)", "Strang (CFL = 0.083)"), "strang_timestep.png"])
plots.append([("SDC", "SDC (CFL = 0.25)", "SDC (CFL = 0.083)"), "sdc_timestep.png"])
plots.append([("Strang", "SDC", "SDC (3 iters)"), "strang_sdc_compare.png"])
plots.append([("SDC", "SDC (3 iters)", "SDC (4 iters)"), "sdc_iters.png"])

for p in plots:

    fig = plt.figure(1)
    fig.clear()

    ax = fig.add_subplot(211)

    dsets, fname = p
    for d in dsets:
        ax.plot(datasets[d][:,0], datasets[d][:,1], label=d)

    ax.set_yscale("log")
    ax.set_ylabel("T [K]")

    ax.legend(frameon=False)


    ax2 = fig.add_subplot(212)

    dsets, fname = p
    for d in dsets:
        ax2.plot(datasets[d][:,0], abs(datasets[d][:,2]), label=d)

    ax2.set_yscale("log")
    ax2.set_ylabel("nuclear energy generation rate [erg/g/s]")
    ax2.set_xlabel("x [cm]")

    ax2.set_ylim(1.e14)

    fig.set_size_inches(7.0, 10.0)

    fig.tight_layout()
    fig.savefig(fname)
