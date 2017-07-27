# Castro Coding Standards


## C++ to Fortran Calling

These are the general conventions that should be followed when calling
a Fortran routine from C++

 * The Fortran routines called by C++ should begin with `ca_`. This is true
   even if these routines are also called by other Fortran routines.
 
 * Passing a FAB to Fortran should be done using the macro BL_TO_FORTRAN_ANYD(fab),
   which passes in the FAB data pointer, and a three-dimensional version of lo and
   hi (in 1D and 2D, the extra dimensions are filled with zeros). In the C++ function
   declaration, you can either use BL_FORT_FAB_ARG_ANYD(fab), or you can explicitly
   list out Real*, const int*, const int*.

 * A separate tile `lo`/`hi` should be passed in to each subroutine, with the expectation
   that the Fortran subroutine will act exactly over the domain (lo, hi). If you want to
   include ghost zones in the calculation, include them in your box using mfi.growntilebox(ng)
   instead of mfi.validbox().

 * Use `bind(C)` on all Fortran routines that are called from C++, with an explicit `name=` clause.
 

## Fortran style

 * use `implicit none` in all routines, even those in a module already
   covered by a blanket `implicit none`

 * use `only` clause when using from a module

 * use `intent` on all arguments


