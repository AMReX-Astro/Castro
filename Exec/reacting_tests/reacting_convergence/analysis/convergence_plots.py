import yt
import numpy as np
import matplotlib.pyplot as plt
import itertools

files = ["react_converge_64_sdc4_plt00301",
         "react_converge_64_sdc_plt00301",
         "react_converge_64_strang_plt00301",
         "react_converge_128_sdc4_plt00601",
         "react_converge_128_sdc_plt00601",
         "react_converge_128_strang_plt00601",
         "react_converge_256_sdc4_plt01201",
         "react_converge_256_sdc_plt01201",
         "react_converge_256_strang_plt01201",
         "react_converge_512_sdc_plt02401",
         "react_converge_512_strang_plt02401"]

varnames = ["density", "Temp", "X(C12)", "magvel"]

class Profile:

    def __init__(self, filename):

        # get the simulation params from the file name
        _tmp = filename.replace("react_converge_", "").replace("plt", "")
        fields = _tmp.split("_")
        self.nx = int(fields[0])
        self.integrator = fields[1]
        self.plotnum = int(fields[2])

        # read in the data
        ds = yt.load(filename)

        xctr = 0.5*(ds.domain_left_edge[0] + ds.domain_right_edge[0])
        L_x = ds.domain_right_edge[0] - ds.domain_left_edge[0]

        ymin = ds.domain_left_edge[1]
        ymax = ds.domain_right_edge[1]

        ray = ds.ray([xctr, ymin, 0], [xctr, ymax, 0])

        isrt = np.argsort(ray["t"])

        self.y = ray["y"][isrt]

        self.data = {}
        for f in varnames:
            self.data[f] = ray[f][isrt]

if __name__ == "__main__":

    profiles = []

    for f in files:
        profiles.append(Profile(f))


    # plot 64 and 256
    d64 = [q for q in profiles if q.nx == 64]
    d256 = [q for q in profiles if q.nx == 256]

    fig, _ax = plt.subplots(2,2)
    axes = list(itertools.chain(*_ax))

    fig.set_size_inches(7.0, 8.0)

    for i, f in enumerate(varnames):
        for itype, c in [("sdc", "C1"), ("sdc4", "C2"), ("strang", "C0")]:

            d1 = [q for q in d64 if q.integrator == itype][-1]
            d2 = [q for q in d256 if q.integrator == itype][-1]

            axes[i].plot(d2.y, d2.data[f], label=itype, ls="-", color=c)
            axes[i].plot(d1.y, d1.data[f], ls=":", color=c)

        axes[i].set_ylabel(f)
        axes[i].set_xlabel("y (cm)")

    plt.legend()
    plt.tight_layout()
    plt.savefig("reacting_convergence.png")

