#!/usr/bin/env python

import matplotlib.pyplot as plt

import detonation

runs = [("res24.0km", 24),
        ("res12.0km", 12),
        ("res6.0km", 6),
        ("res3.0km", 3),
        ("res1.5km", 1.5),
        ("res0.1875km", 0.1875)] #,
       #("res0.024km", 0.024)] #,
#("res0.003km", 0.003)]

nsb1_runs = [("res24.0km_noshockburn_1", 24),
             ("res12.0km_noshockburn_1", 12),
             ("res6.0km_noshockburn_1", 6),
             ("res3.0km_noshockburn_1", 3),
             ("res1.5km_noshockburn_1", 1.5),
             ("res0.1875km_noshockburn_1", 0.1875)] #,

nsb23_runs = [("res24.0km_noshockburn_0.666", 24),
              ("res12.0km_noshockburn_0.666", 12),
              ("res6.0km_noshockburn_0.666", 6),
              ("res3.0km_noshockburn_0.666", 3),
              ("res1.5km_noshockburn_0.666", 1.5),
              ("res0.1875km_noshockburn_0.666", 0.1875)] #,

res = []
v = []
dv = []

for ddir, dx in runs:
    res.append(dx)
    d = detonation.Detonation(ddir)
    v.append(d.v)
    dv.append(d.v_sigma)

nsb23_res = []
nsb23_v = []
nsb23_dv = []

for ddir, dx in nsb23_runs:
    nsb23_res.append(dx)
    d = detonation.Detonation(ddir)
    nsb23_v.append(d.v)
    nsb23_dv.append(d.v_sigma)

nsb1_res = []
nsb1_v = []
nsb1_dv = []

for ddir, dx in nsb1_runs:
    nsb1_res.append(dx)
    d = detonation.Detonation(ddir)
    nsb1_v.append(d.v)
    nsb1_dv.append(d.v_sigma)


fig, ax = plt.subplots()

ax.errorbar(res, v, yerr=dv, fmt="o", label="burning in shocks allowed")
ax.errorbar(nsb23_res, nsb23_v, yerr=nsb23_dv, fmt="d", label="no shock burning (f=2/3)")
ax.errorbar(nsb1_res, nsb1_v, yerr=nsb1_dv, fmt="d", label="no shock burning (f=1)")

ax.set_xlabel(r"$\Delta x$ (km)")
ax.set_ylabel(r"$v_\mathrm{det}$ (cm / s)")

ax.legend()

ax.set_xscale("log")

fig.savefig("det_speeds.png")
