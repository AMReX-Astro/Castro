#!/usr/bin/env python

import matplotlib.pyplot as plt

import detonation

runs = [("res12.288km", 12.288),
        ("res1.536km", 1.536),
        ("res0.192km", 0.192),
        ("res0.024km", 0.024)] #,
#("res0.003km", 0.003)]

nsb_runs = [("res12.288km_noshockburn", 12.288),
            ("res1.536km_noshockburn", 1.536),
            ("res0.192km_noshockburn", 0.192),
            ("res0.024km_noshockburn", 0.024)] #,

res = []
v = []
dv = []

for ddir, dx in runs:
    res.append(dx)
    d = detonation.Detonation(ddir)
    v.append(d.v)
    dv.append(d.v_sigma)

nsb_res = []
nsb_v = []
nsb_dv = []

for ddir, dx in nsb_runs:
    nsb_res.append(dx)
    d = detonation.Detonation(ddir)
    nsb_v.append(d.v)
    nsb_dv.append(d.v_sigma)


fig, ax = plt.subplots()

ax.errorbar(res, v, yerr=dv, fmt="o", label="burning in shocks allowed")
ax.errorbar(nsb_res, nsb_v, yerr=nsb_dv, fmt="d", label="no shock burning")

ax.set_xlabel(r"$\Delta x$ (km)")
ax.set_ylabel(r"$v_\mathrm{det}$ (cm / s)")

ax.legend()

ax.set_xscale("log")

fig.savefig("det_speeds.png")
