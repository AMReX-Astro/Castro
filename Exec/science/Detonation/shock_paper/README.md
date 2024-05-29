# Shock Burning Experiments

This directory is meant to explore shock burning with detonation.  Compile as:

```
make USE_SIMPLIFIED_SDC=TRUE USE_SHOCK_VAR=TRUE NETWORK_DIR=aprox13 -j 4
```

Then setup a suite of simulations with the following resolutions (the
`inputs-det-x.subch_base`) here is for the coarsest run:


| resolution   |  base grid  |  levels (4x jumps)  |
| ------------ | ----------- | ------------------- |
|   12.288 km  |       48    |         1           |
|    1.536 km  |      384    |         1           |
|    0.192 km  |     3072    |         1           |
|     2400 cm  |     6144    |         2           |
|      300 cm  |    12288    |         3           |




