# `nse_test`

This is a simple test problem designed to explore how well hydro and
reactions are coupled when a system enters NSE.

This version is based on ``reacting_convergence`` (which is in turn
based on ``acoustic_pulse_general``), but using the ``aprox19``
network with the NSE table enabled.

You can run the simplified-SDC convergence test with the script
`convergence_simplified_sdc_w_vel.sh`

The script `create_pretty_tables.py` will take the 2 (or 3) output
files and make a single LaTeX-formatted table of the results.  Use
the `--simple` argument to format in plaintext.

