# Input parameters to simulate WD mergers

### Problem set up
```
problem.problem = 1
```
### Easiest way to make the two stars merge

```
problem.roche_radius_factor = 1.0e0
```
This is used to calculate the orbital distance and the rotational period. (Ref: https://articles.adsabs.harvard.edu/pdf/1983ApJ...268..368E)

### Desired masses or densities

Castro either uses the initial WD masses or the initial central densities to generate desired stars. Users can specify one of these and leave the other parameter negative. If both are provided, it will use the central densities to generate the model and assume masses as an estimate for the purpose of determining the envelope mass boundary. This is established [here](https://github.com/AMReX-Astro/Castro/blob/10534968cf09ba73f1044186389740559adc4188/Util/model_parser/model_parser.H#L208).

**Primary**
```
problem.mass_P = 1.0e0
problem.central_density_P = -1.0e0
```

**Secondary**
```
problem.mass_S = 1.0e0
problem.central_density_S = -1.0e0
```

### Desired compositions

Inspired by Dan et al. 2012, we classify the WDs using their masses and the maximum values that they can hold, and set their compositions accordingly.

Stars with masses less than 0.45 $M_{\odot}$ are given a pure He composition. This maximum limit can be changed using:
```
problem.max_he_wd_mass =  0.45e0
```

To simulate stars with pure CO cores and a thin He shell between 0.45 $M_{\odot}$ and 0.6 $M_{\odot}$ we use:

```
problem.max_he_wd_mass =  0.45e0
problem.max_hybrid_wd_mass = 0.6e0

problem.hybrid_wd_c_frac = 0.50e0
problem.hybrid_wd_o_frac = 0.50e0
problem.hybrid_wd_he_shell_mass = 0.10e0
```
On the other hand, to simlute stars with pure CO cores and a He shell between 0.6 $M_{\odot}$ and 1.05 $M_{\odot}$, and with different CO mass fractions, we use:
```
problem.max_hybrid_wd_mass = 0.6e0
problem.max_co_wd_mass = 1.05e0

problem.co_wd_c_frac = 0.40e0
problem.co_wd_o_frac = 0.60e0
problem.co_wd_he_shell_mass = 0.0e0
```
**NOTE:** To create stars with masses higher than the default maximum values, we need too increase the limits otherwise the conditions will not be satisfied and hence, we won't get the desired compositions.

### Stellar temperature
We can set the initial stellar temperature to the desired value using
```
problem.stellar_temp = 1.0e7
```

### Reactions/Network
In general, we want to disable burning in shocks, especially for the case of double-detonations (Ref: https://iopscience.iop.org/article/10.3847/2515-5172/add686). We use the make variable `USE_SHOCK_VAR = TRUE` and additionally the input parameter
```
castro.disable_shock_burning = 1
```