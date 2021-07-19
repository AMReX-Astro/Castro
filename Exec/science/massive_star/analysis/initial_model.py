# take a raw inputs file for aprox21 and convert to the aprox19
# nuclei.  This means getting rid of Cr56 and Fe56 (by lumping them
# into Ni56)

import numpy as np
import matplotlib.pyplot as plt

def find_r_for_rho(r, rho, rho_want):
    idx = np.where(rho < rho_want)[0][0]
    return r[idx]


#file = "../15m_500_sec.aprox19.hse.6400"
file = "../15m_500_sec.aprox19.hse.20km"

data = np.loadtxt(file)

print(data.shape)

# now manually read to get the variable names

# the first column is position

names = ["r"]

with open(file) as f:
    for n, line in enumerate(f):
        if line.startswith("# num"):
            continue

        if line.startswith("# npts"):
            continue

        if not line.startswith("#"):
            break

        names.append(line.split()[-1].strip())

# now make plots

idens = names.index("density")
itemp = names.index("temperature")
iye = names.index("Ye")

fig = plt.figure()
ax = fig.add_subplot(211)

l1 = ax.plot(data[:,0], data[:,idens], label="density")
l2 = ax.plot(data[:,0], data[:,itemp], label="temperature")

# show where the refinement kicks in 
rho_refine = 1.e4

r_refine = find_r_for_rho(data[:,0], data[:,idens], rho_refine)

print(r_refine)

ax.axvline(r_refine, color="0.25", ls=":")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("r [cm]")
ax.set_ylabel(r"$\rho~[\rm{g/cm^3}]$, $T~[K]$")
ax2 = ax.twinx()

ax2.set_ylabel(r"$Y_e$")

l3 = ax2.plot(data[:,0], data[:,iye], color="C2", label="Ye")

lns = l1 + l2 + l3
labs = [l.get_label() for l in lns]

ax.legend(lns, labs, frameon=False)


ax = fig.add_subplot(212)

threshold = 0.1


for n, var in enumerate(names):
    print(n, var)

    if var in ["r", "density", "temperature", "pressure", "Ye"]:
        continue

    Xmax = data[:,n].max()

    if Xmax > threshold:
        if Xmax > 0.5:
            lw = 2
        else:
            lw = 1
        ax.plot(data[:,0], data[:,n], label=var, lw=lw)


ax.axvline(r_refine, color="0.25", ls=":")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("r [cm]")
ax.set_ylabel("mass fraction")

ax.legend(frameon=False, fontsize="small", ncol=2)

fig.set_size_inches((8, 12))

fig.tight_layout()

fig.savefig("initial_model.png")

