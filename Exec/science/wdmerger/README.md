# Input parameters to simulate WD mergers

### Easiest way to start a merger:

```
problem.roche_radius_factor = 1.0e0
```
### This is used to calculate the orbital distance and the rotational period. (Ref: https://articles.adsabs.harvard.edu/pdf/1983ApJ...268..368E)

### Castro either uses the inital WD masses or the initial central densities to generate desired stars. Users can specify one of these and leave the other parameter negative. If both are provided, it will use the central densities to generate the model and assume masses as an estimate for the purpose of determining the envelope mass boundary. This is established [here](https://github.com/AMReX-Astro/Castro/blob/10534968cf09ba73f1044186389740559adc4188/Util/model_parser/model_parser.H#L208).

### To set the CO mass fractions and He shell mass for WDs with  masses between 0.45 M_sun and 0.60 M_sun, use:
```
problem.hybrid_wd_c_frac = 0.00e0
problem.hybrid_wd_o_frac = 0.00e0
problem.hybrid_wd_he_shell_mass = 0.00e0
```
### To set the CO mass fractions and He shell mass for WDs with masses between 0.6 M_sun and 1.05 M_sun, use:
```
problem.co_wd_c_frac = 0.00e0
problem.co_wd_o_frac = 0.00e0
problem.co_wd_he_shell_mass = 0.00e0
```

