# Citing Castro

There are a number of Castro papers that describe parts of the
algorithm.  We ask that you cite all of the appropriate papers
describing the capabilities you used.

## General code use

If you use Castro, we appreciate you citing the most recent code paper:

```
@article{Almgren2020,
  doi = {10.21105/joss.02513},
  url = {https://doi.org/10.21105/joss.02513},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {54},
  pages = {2513},
  author = {Ann Almgren and Maria Barrios Sazo and John Bell and Alice Harpole and Max Katz and Jean Sexton and Donald Willcox and Weiqun Zhang and Michael Zingale},
  title = {CASTRO: A Massively Parallel Compressible Astrophysics Simulation Code},
  journal = {Journal of Open Source Software}
}
```

You are welcome to cite the original code paper as well (which
provides some details on the algorithmic implementations):

```
    @ARTICLE{2010ApJ...715.1221A,
       author = {{Almgren}, A.~S. and {Beckner}, V.~E. and {Bell},
                      J.~B. and {Day}, M.~S. and {Howell}, L.~H. and
                      {Joggerst}, C.~C. and {Lijewski}, M.~J. and
                      {Nonaka}, A. and {Singer}, M. and {Zingale}, M.},
        title = "{CASTRO: A New Compressible Astrophysical
                      Solver. I. Hydrodynamics and Self-gravity}",
      journal = {\apj},
    archivePrefix = "arXiv",
       eprint = {1005.0114},
     primaryClass = "astro-ph.IM",
     keywords = {equation of state, gravitation, hydrodynamics, methods:
                      numerical, nuclear reactions, nucleosynthesis,
                      abundances},
         year = 2010,
        month = jun,
       volume = 715,
        pages = {1221-1238},
          doi = {10.1088/0004-637X/715/2/1221},
       adsurl = {http://adsabs.harvard.edu/abs/2010ApJ...715.1221A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

You should also cite the zenodo DOI for the code release.  A bibtex
entry for the latest release can be found here:

https://doi.org/10.5281/zenodo.2301848

## Radiation hydrodynamics

If you use the radiation hydrodynamics capabilities, please additionally
cite the following:

```
    @ARTICLE{2011ApJS..196...20Z,
       author = {{Zhang}, W. and {Howell}, L. and {Almgren}, A. and
                      {Burrows}, A. and {Bell}, J.},
        title = "{CASTRO: A New Compressible Astrophysical
                      Solver. II. Gray Radiation Hydrodynamics}",
      journal = {\apjs},
    archivePrefix = "arXiv",
       eprint = {1105.2466},
     primaryClass = "astro-ph.IM",
     keywords = {diffusion, hydrodynamics, methods: numerical, radiative
                      transfer},
         year = 2011,
        month = oct,
       volume = 196,
          eid = {20},
        pages = {20},
          doi = {10.1088/0067-0049/196/2/20},
       adsurl = {http://adsabs.harvard.edu/abs/2011ApJS..196...20Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

```
    @ARTICLE{2013ApJS..204....7Z,
       author = {{Zhang}, W. and {Howell}, L. and {Almgren}, A. and
                      {Burrows}, A. and {Dolence}, J. and {Bell}, J.},
        title = "{CASTRO: A New Compressible Astrophysical
                      Solver. III. Multigroup Radiation Hydrodynamics}",
      journal = {\apjs},
    archivePrefix = "arXiv",
       eprint = {1207.3845},
     primaryClass = "astro-ph.IM",
     keywords = {diffusion, hydrodynamics, methods: numerical, radiative
                      transfer },
         year = 2013,
        month = jan,
       volume = 204,
          eid = {7},
        pages = {7},
          doi = {10.1088/0067-0049/204/1/7},
       adsurl = {http://adsabs.harvard.edu/abs/2013ApJS..204....7Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

## Gravity and rotation

If you use the Poisson gravity solver or rotation, please cite the
following, which included a lot of improvements to the coupling of
hydro and gravity:

```
    @ARTICLE{2016ApJ...819...94K,
       author = {{Katz}, M.~P. and {Zingale}, M. and {Calder}, A.~C. and
                      {Swesty}, F.~D. and {Almgren}, A.~S. and {Zhang},
                      W.},
        title = "{White Dwarf Mergers on Adaptive Meshes. I. Methodology
                      and Code Verification}",
      journal = {\apj},
    archivePrefix = "arXiv",
       eprint = {1512.06099},
     primaryClass = "astro-ph.HE",
     keywords = {hydrodynamics, methods: numerical, supernovae: general,
                      white dwarfs},
         year = 2016,
        month = mar,
       volume = 819,
          eid = {94},
        pages = {94},
          doi = {10.3847/0004-637X/819/2/94},
       adsurl = {http://adsabs.harvard.edu/abs/2016ApJ...819...94K},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

## GPUs and scaling

For CPU performance numbers, please cite:

```
    @INPROCEEDINGS{2018JPhCS1031a2024Z,
           author = {{Zingale}, M. and {Almgren}, A.~S. and {Barrios Sazo}, M.~G. and
             {Beckner}, V.~E. and {Bell}, J.~B. and {Friesen}, B. and
             {Jacobs}, A.~M. and {Katz}, M.~P. and {Malone}, C.~M. and
             {Nonaka}, A.~J. and {Willcox}, D.~E. and {Zhang}, W.},
            title = "{Meeting the Challenges of Modeling Astrophysical Thermonuclear Explosions: Castro, Maestro, and the AMReX Astrophysics Suite}",
         keywords = {Astrophysics - Instrumentation and Methods for Astrophysics},
        booktitle = {Journal of Physics Conference Series},
             year = 2018,
           series = {Journal of Physics Conference Series},
           volume = {1031},
            month = may,
              eid = {012024},
            pages = {012024},
              doi = {10.1088/1742-6596/1031/1/012024},
    archivePrefix = {arXiv},
           eprint = {1711.06203},
     primaryClass = {astro-ph.IM},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2018JPhCS1031a2024Z},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

For GPU performance, please cite:

```
    @ARTICLE{2020arXiv200705218K,
           author = {{Katz}, Max P. and {Almgren}, Ann and {Barrios Sazo}, Maria and
             {Eiden}, Kiran and {Gott}, Kevin and {Harpole}, Alice and
             {Sexton}, Jean M. and {Willcox}, Don E. and {Zhang}, Weiqun and
             {Zingale}, Michael},
            title = "{Preparing Nuclear Astrophysics for Exascale}",
          journal = {arXiv e-prints},
         keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - High Energy Astrophysical Phenomena},
             year = 2020,
            month = jul,
              eid = {arXiv:2007.05218},
            pages = {arXiv:2007.05218},
    archivePrefix = {arXiv},
           eprint = {2007.05218},
     primaryClass = {astro-ph.IM},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200705218K},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

## Spectral deferred corrections

For the 2nd and 4th order SDC coupling of hydro and reactions, please cite:

```
    @ARTICLE{2019ApJ...886..105Z,
           author = {{Zingale}, M. and {Katz}, M.~P. and {Bell}, J.~B. and {Minion}, M.~L. and
             {Nonaka}, A.~J. and {Zhang}, W.},
            title = "{Improved Coupling of Hydrodynamics and Nuclear Reactions via Spectral Deferred Corrections}",
          journal = {\apj},
         keywords = {Hydrodynamics, Astrophysical fluid dynamics, Computational methods, Computational astronomy, Astronomy software, Nuclear astrophysics, Nucleosynthesis, Stellar nucleosynthesis, Physics - Computational Physics, Astrophysics - Instrumentation and Methods for Astrophysics},
             year = 2019,
            month = dec,
           volume = {886},
           number = {2},
              eid = {105},
            pages = {105},
              doi = {10.3847/1538-4357/ab4e1d},
    archivePrefix = {arXiv},
           eprint = {1908.03661},
     primaryClass = {physics.comp-ph},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2019ApJ...886..105Z},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

