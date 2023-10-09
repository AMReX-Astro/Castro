This is an exact Riemann solver for a general equation of state.  It
follows the outline for Colella & Glaz 1985, section 1.  This is too
slow to use in an actual hydro run, but instead is intended to
generate exact solutions to the Riemann problem for comparison with
Castro shocktube output.  Several inputs files for Helmholtz EOS-based
shocktubes are provided.

This solver is used in Zingale & Katz (2015):

https://ui.adsabs.harvard.edu/abs/2015ApJS..216...31Z/abstract

and more details are given there.

To build the solver, simply type 'make' in this directory.
