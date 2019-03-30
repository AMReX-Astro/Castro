******************
Coding Conventions
******************

Castro development should do its best to adhere to the following coding
style guidelines.

General
=======

* Indentation should be spaces, not tabs, with 4 spaces preferred.


C++
===

* Conditional should appear as:

  .. code:: c++

     if (condition)
     {
         ...
     }
     else if (condition)
     {
         ...
     }
     else
     {
         ...
     }


C++ to Fortran
==============

* All C to Fortran interfaces should use the ISO C binding.  The
  Fortran subroutine should use

  .. code:: fortran

     subroutine subroutine_name(args) bind(C, name="subroutine_name)

* Data passed by reference from C++ should use ``*`` and not ``&``.

* Scalars should be passed by value, using the Fortran ``value`` attribute.

* Fortran routines that are called from C++ should being with ``ca_``.
  This is true even if these routines are also called by other
  Fortran.

* A separate tile ``lo`` and ``hi`` should be passed in to each
  subroutine, with the expectation that the Fortran subroutine will
  act exactly over the domain (lo, hi). If you want to include ghost
  zones in the calculation, include them in your box using
  ``mfi.growntilebox(ng)`` instead of ``mfi.validbox()``.


Fortran
=======

* Put all routines in a module—this ensures that argument lists are
  checked.

* Use the ``only`` clause in module ``use`` statements to explicitly
  make clear what is being accessed.

* In a module, there should be no "top-level" ``use`` statements (with
  the exception of getting access to the ``rt`` type).  Instead each
  function / subroutine in the module  should use what it needs directly.

* New Fortran files should have the .F90 file extension, not the .f90
  file extension, so that they can be preprocessed.


Documentation
=============

C++ routines should use Doxygen style comments.

* C++ functions should be documented in the header file using the style:

  .. code:: c++

     ///
     /// Description of the function
     ///
     /// @param bar       Brief description of the variable
     ///
     void foo(int bar) { ...

* Member variables can either be documented using the above style of comment block or
  with a brief inline description:

  .. code:: c++

     int var; ///< Brief description after the variable

Fortran functions should use Sphinx style comments

* Fortran functions should be documented by placing a comment block
  immediately after their prototype (i.e. `without` a line in betwen ) using the style:

  .. code:: fortran

     subroutine foo(bar)
       ! Description of the function

       use some_module

       implicit none

       integer, intent(inout) :: bar   ! Brief description of bar
       ...

  Documentation for modules should be similarly formatted, with the comment block again
  coming `immediately` after the module definition.
