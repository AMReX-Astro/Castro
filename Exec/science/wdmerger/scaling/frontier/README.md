# wdmerger scaling on Frontier

This explores a 12.5 km resolution wdmerger simulation using the
Pakmor initial conditions.

We consider 3 different gridding strategies:

* 256^3 base + 3 AMR levels, each a jump of 4

* 512^3 base + 3 AMR levels with jumps of 4, 4, 2

* 1024^3 base + 2 AMR levels with jumps of 4, 4

The inputs file here is setup for the 256^3 base.

We report the total evolution time excluding initialization that is
output by Castro at the end of the run.

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

Note that for the 256^3 base grid, on 64 nodes, the grid structure is:

INITIAL GRIDS
  Level 0   512 grids  16777216 cells  100 % of domain
            smallest grid: 32 x 32 x 32  biggest grid: 32 x 32 x 32
  Level 1   96 grids  3145728 cells  0.29296875 % of domain
            smallest grid: 32 x 32 x 32  biggest grid: 32 x 32 x 32
  Level 2   674 grids  38797312 cells  0.05645751953 % of domain
            smallest grid: 32 x 32 x 32  biggest grid: 64 x 32 x 32
  Level 3   7247 grids  1428029440 cells  0.03246963024 % of domain
            smallest grid: 32 x 32 x 32  biggest grid: 64 x 64 x 64

So only a small amount of the finest grid is refined in this problem.
