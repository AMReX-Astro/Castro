******************
Coding Conventions
******************


C++
===



Fortran
=======

* Put all routines in a module—this ensures that argument lists are
  checked.

* Use the ``only`` clause in module ``use`` statements to explicitly
  make clear what is being accessed.

* Fortran routines that are called from C++ should being with ``ca_``
