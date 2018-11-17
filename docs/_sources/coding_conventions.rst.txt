******************
Coding Conventions
******************


C++
===


C++ to Fortran
==============

* Fortran routines that are called from C++ should being with ``ca_``.
  This is true even if these routines are also called by other
  Fortran.

* A separate tile ``lo`` and ``hi`` should be passed in to each
  subroutine, with the expectation that the Fortran subroutine will
  act exactly over the domain (lo, hi). If you want to include ghost
  zones in the calculation, include them in your box using
  ``mfi.growntilebox(ng)`` instead of ``mfi.validbox()``.

* macros?


Fortran
=======

* Put all routines in a module—this ensures that argument lists are
  checked.

* Use the ``only`` clause in module ``use`` statements to explicitly
  make clear what is being accessed.

