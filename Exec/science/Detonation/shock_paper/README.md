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



