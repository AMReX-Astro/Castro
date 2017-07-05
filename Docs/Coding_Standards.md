# Castro Coding Standards


## C++ to Fortran Calling

These are the general conventions that should be followed when calling
a Fortran routine from C++

 * the Fortran routines called by C++ should begin with `ca_`
 
 * for each FAB, pass in array name, `lo(3)`, `hi(3)`, regardless of dimension

 * use `BL_TO_FORTRAN_3D` when calling Fortran from C++

 * a separate tile `lo`/`hi` should be passed in too

 * use `bind(C)` on the Fortran routines, with an explicit `name=` clause
 

## Fortran style

 * use `implicit none` in all routines, even those in a module already
   covered by a blanket `implicit none`

 * use `only` clause when using from a module

 * use `intent` on all arguments


