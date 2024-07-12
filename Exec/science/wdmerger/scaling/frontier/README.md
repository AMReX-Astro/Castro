# wdmerger scaling on Frontier

This explores a 12.5 km resolution wdmerger simulation using the
Pakmor initial conditions.

We consider 3 different gridding strategies:

* 256^3 base + 3 AMR levels, each a jump of 4

* 512^3 base + 3 AMR levels with jumps of 4, 4, 2

* 1024^3 base + 2 AMR levels with jumps of 4, 4

The inputs file here is setup for the 256^3 base.

Some general observations:

* We seem to do well with `max_grid_size` set to 64 or 128, but not 96

* At large node counts, it really doesn't matter which of the gridding
  strategies we use, since there is plenty of work to go around.  The
  main consideration would be that the larger coarse grid would make
  the plotfiles bigger.

* We seem to benefit from using `castro.hydro_memory_footprint_ratio=3`

* There really is no burning yet, since this is early in the
  evolution, so we would expect scaling to improve as the stars
  interact (more grids) and burning begins (more local work).
