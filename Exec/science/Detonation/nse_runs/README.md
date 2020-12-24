# `nse_detonations`

This directory contains a set of scripts that setup, run, and analyze
a set of detonations, comparing Strang splitting to simplified SDC.

The basic usage is:

#. in a `screen` session of similar do:
   ```
   ./setup_runs.py
   ```
   this will create the run directories, copy the needed files, and
   run the jobs in parallel using the python Pool mechanism.

#. while the jobs are running, you can check the status by doing:
   ```
   ./show_status.py
   ```

#. once the runs are finised, you can make the suite of plots by
   ```
   ./make_plots.py
   ```


