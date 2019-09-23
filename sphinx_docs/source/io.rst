.. _ch:io:

************
Parallel I/O
************

Both checkpoint files and plotfiles are really directories containing
subdirectories: one subdirectory for each level of the AMR hierarchy.
The fundamental data structure we read/write to disk is a ``MultiFab``,
which is made up of multiple FAB’s, one FAB per grid. Multiple
``MultiFab`` s may be written to each directory in a checkpoint file.
``MultiFab`` s of course are shared across CPUs; a single ``MultiFab`` may be
shared across thousands of CPUs. Each CPU writes the part of the
``MultiFab`` that it owns to disk, but they don’t each write to their own
distinct file. Instead each MultiFab is written to a runtime
configurable number of files :math:`N` (:math:`N` can be set in the inputs file as the
parameter ``amr.checkpoint_nfiles`` and ``amr.plot_nfiles``; the
default is 64). That is to say, each ``MultiFab`` is written to disk
across at most :math:`N` files, plus a small amount of data that gets written
to a header file describing how the file is laid out in those :math:`N` files.

What happens is :math:`N` CPUs each opens a unique one of the :math:`N` files into
which the ``MultiFab`` is being written, seeks to the end, and writes
their data. The other CPUs are waiting at a barrier for those :math:`N`
writing CPUs to finish. This repeats for another :math:`N` CPUs until all the
data in the ``MultiFab`` is written to disk. All CPUs then pass some data
to CPU 0 which writes a header file describing how the ``MultiFab`` is
laid out on disk.

We also read ``MultiFabs`` from disk in a “chunky” manner, opening only :math:`N`
files for reading at a time. The number :math:`N`, when the ``MultiFab`` s were
written, does not have to match the number :math:`N` when the ``MultiFab`` s are
being read from disk. Nor does the number of CPUs running while
reading in the ``MultiFab`` need to match the number of CPUs running when
the ``MultiFab`` was written to disk.

Think of the number :math:`N` as the number of independent I/O pathways in
your underlying parallel filesystem. Of course a “real” parallel
filesytem should be able to handle any reasonable value of :math:`N`. The
value -1 forces :math:`N` to the number of CPUs on which you’re
running, which means that each CPU writes to a unique file, which can
create a very large number of files, which can lead to inode issues.
