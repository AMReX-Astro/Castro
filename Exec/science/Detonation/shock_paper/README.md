# Shock Burning Experiments

This directory is meant to explore shock burning with detonation.  Compile as:

```
make USE_SIMPLIFIED_SDC=TRUE USE_SHOCK_VAR=TRUE NETWORK_DIR=aprox13 -j 4
```

Then the script `setup_runs.py` will setup a suite of simulations with
the following resolutions into separate directories (using the
`inputs-shock-burn.template`):


| resolution   |  base grid  |  levels (4x jumps)  |
| ------------ | ----------- | ------------------- |
|       24 km  |       48    |         1           |
|       12 km  |       96    |         1           |
|        6 km  |      192    |         1           |
|        3 km  |      384    |         1           |
|      1.5 km  |      768    |         1           |
|   0.1875 km  |     6144    |         1           |
|  2343.74 cm  |    12288    |         2           |

you can set the value of the shock detection threshold there
and the directory names will reflect that setting.

## plotting

The following scripts can make useful plots (some use the
`detonation.py` module as support):

* `det_speed_comp.py` : plot detonation speed vs. resolution, using
  simple differencing to estimate the detonation speed from the last 3
  plotfiles.

* `profile_compare.py` : given a list of pruns (from different
  resolutions), make a plot showing all of their profiles together.

* `profiles.py` : for a single run, plot profiles of T and enuc for
  several times.

* `show_shock_flag.py` : simply plot the shock variable on top of T
  and enuc profiles.

* `zoom_summary.py` : given a list of runs (from different
  resolutions), plot the last plotfile from each run, zoomed in on
  where the peak energy generation is.
