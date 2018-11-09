Castro requires both a C compiler (with support for C 11) and
a Fortran compiler with 2003+ support.

General Compilers
=================

GCC
---

The GCC compilers are the preferred compilers on Linux machines.
Castro runs without issue with GCC 7.x and GCC 8.x.

Intel
-----

Intel compilers 17.x, 18.x, and 19.x produce internal compiler errors
and should not be used.

PGI
---

The PGI compilers (16.4–16.10) are also tested with Castro and have
no known issues. Since PGI uses the GCC header files, a compatible
version of GCC is needed. This is particularly problematic with C 11
support. At the moment, with PGI 16.10 compilers the latest GCC compiler
that has compatible C 11 support is GCC 4.9.4.

Working at OLCF (ORNL)
======================

The build system needs python 2.7+ or 3.5+. The default python at OLCF
is too old, so you will need to manually load a newer python as:

::

    module load python

Cray
----

Note, you will need to swap the default PGI environment to use Cray:

::

    module swap PrgEnv-pgi PrgEnv-cray

cce/8.5.7 through 8.6.3 fail to compile AMReX. The default version at OLCF as of May 2018,
cce/8.6.4, resolves this issue.

.. _pgi-1:

PGI
---

If you are using OpenACC on the GPUs, then you need to use the PGI
compilers instead. It is best to use the latest available PGI
compilers, which are 17.7 at this writing.

::

    module swap pgi pgi/17.7

These seem to build the code without problems.

This is perhaps out of date:
Finally, to use OpenACC, you need to make the CUDA libraries available to the PGI compiler:

::

    module load craype-accel-nvidia35

Automatic Restarting and Archiving of Data
------------------------------------------

See the Maestro User’s Guide

Working at NERSC
================

General issues
--------------

Edison compilers
----------------

-  Intel compilers

   The default compilers on edison are the Intel compilers.

   The Intel 17.0.2 compilers produce code that crashes at runtime when self-gravity
   is used (the error is in the F_MG directory. This appears to be an issue since
   16.0.2.

   At the moment, avoid Intel compilers on Edison.

-  Cray compilers

   -  The Cray compilers from versions 8.5.7 to 8.6.3 have issues compiling AMReX. This
      can be resolved by loading version 8.6.4 or later. As of May 2018, you can do

      ::

                module swap cce cce/8.6.5
              

-  GNU

   The GNU compilers (6.3) seem to work fine.

Cori (KNL) compilers
--------------------

To compile with Cray compilers on Cori, we need to swap the compiler
wrappers to use the AVX-512 instruction set supported on the Intel Phi
processors instead of the AVX-2 extensions used by the Intel Haswell
architecture. This is done as:

::

    module swap craype-{haswell,mic-knl}

It could happen that even when the various verbosities are set to 0, when using several nodes (more than 64) in a run compiled with Intel, the follwing error shows:

::

    "forrtl: severe (40): recursive I/O operation, unit -1, file unknown"

Seems like the error is due to all threads printing to stdout. Adding the following to the inputs file, prevents this error to occur:

::

    castro.print_fortran_warnings = 0

Hypre and radiation
-------------------

On edison, the Cray *Third Party Scientific Libraries* provide
hypre in a form that works directly with the compiler wrappers
used on that machine (CC, ftn, :math:`\ldots`). To use this,
simply do:

::

    module load cray-tpsl

There is no need to set HYPRE_DIR, but note however that the
dependency checker script (BoxLib/Tools/C_scripts/mkdep) will
complain about:

::

    /path/to/Hypre--with-openmp/include does not exist

This can be ignored an compilation will finish. If you do wish to
silence it, you can set HYPRE_DIR to the path shown by

::

    module show cray-tpsl

as

::

    export HYPRE_DIR=${CRAY_TPSL_PREFIX_DIR}

This path will change dynamically to reflect which compiler programming
environment you have loaded. (You can also see that this is the path
sent to the compilation by doing ftn -craype-verbose).

Running jobs
------------

edison is configured with 24 cores per node split between two Intel
IvyBridge 12-core processors. Each processor connects to 1/2 of the
node’s memory and is called a NUMA node, so there are 2 NUMA nodes per
edison node. Best performance is seen when running with 6 or 12 threads.

Jobs should be run in your $SCRATCH or $CSCATCH directory.
By default, SLURM will change directory into the submission directory.

A sample job submission script, edison.MPI.OMP.slurm is in
Castro/Util/job_scripts/edison/, and includes logic to
automatically add the correct restart options to the run to continue a
simulation from the last checkpoint file in the submission directory.

To chain jobs, such that one queues up after the previous job finished,
use the chainslurm.sh script in that same directory:

::

    chainslurm.sh jobid number script

where jobid is the existing job you want to start you chain
from, number is the number of new jobs to chain from this
starting job, and script is the job submission script to use
(the same one you used originally most likely). You can view the job
dependency using:

::

    squeue -l -j job-id                                                             

where job-id is the number of the job.

Jobs are submitted with sbatch. A job can be canceled using
scancel, and the status can be checked using squeue -u
*username*.

Archiving data to HPSS
----------------------

The script edison.xfer.slurm in
Castro/Util/job_scripts/edison/ can be used to archive data to
HPSS automatically. This is submitted to the xfer queue and
runs the script process.xrb which continually looks for output
and stores it to HPSS.

To use the scripts, first create a directory in HPSS that has the same
name as the directory on lustre you are running in (just the directory
name, not the full path). E.g. if you are running in a directory
call wdconvect_run, then do:

::

    hsi                                                                             
    mkdir wdconvect_run                                                             

(Note: if the hsi command prompts you for your password, you will need
to talk to the NERSC help desk to ask for password-less access to
HPSS).

The script process.xrb is called from the xfer job and will
run in the background and continually wait until checkpoint or
plotfiles are created (actually, it always leaves the most recent one
alone, since data may still be written to it, so it waits until there
are more than 1 in the directory).

Then the script will use htar to archive the plotfiles and
checkpoints to HPSS. If the htar command was successful, then
the plotfiles are copied into a plotfile/ subdirectory. This is
actually important, since you don’t want to try archiving the data a
second time and overwriting the stored copy, especially if a purge
took place. The same is done with checkpoint files.

Additionally, if the ftime executable is in your path (
ftime.f90 lives in BoxLib/Tools/Postprocessing/F_src/), then
the script will create a file called ftime.out that lists the
name of the plotfile and the corresponding simulation time.

Finally, right when the job is submitted, the script will tar up all
of the diagnostic files, ftime.out, submission script, inputs
and probin, and archive them on HPSS. The .tar file is given a
name that contains the date-string to allow multiple archives to
co-exist.

When process.xrb is running, it creates a lockfile (called
process.pid) that ensures that only one instance of the script
is running at any one time. Sometimes if the machine crashes, the
process.pid file will be left behind, in which case, the script
aborts. Just delete that if you know the script is not running.

Jobs in the xfer queue start up quickly. The best approach is
to start one as you start your main job (or make it dependent on the
main job). The sample process.xrb script will wait for output
and then archive it as it is produced, using the techniques described
for titan above.

To check the status of a job in the xfer queue, use:

::

    squeue -u username -M all                                                       

Working at LANL
===============

For the following LANL systems, which all have access to a joint file system,

-  Cielito

-  Conejo

-  Lightshow

-  Moonlight

-  Pinto

-  Wolf

-  Mustang

-  Trinitite

the following steps are needed to get Castro compiling (reported by Platon Karpov, 6/13/2016):

| 1) Compile on the login node (thus host is lanl.gov)
| 2) Do *not* set env variable ``BOXLIB_USE_MPI_WRAPPERS`` at all
| 3) Execute ``module load python-epd`` at the command line, but do *not* load anaconda

Scaling
=======

Data from scaling studies is archived in Castro/Docs/ManagingJobs/scaling/

Needs to be updated

Gotyas
======

#. 3/4/16: The default version loaded of Python on Mira is not
   recent enough to support the Python scripts in our build system. Add
   +python to your .soft to fix this.

#. 2/18/16: The default version loaded of Python on Titan (2.6.9)
   is not recent enough to support the Python scripts in our build
   system. At the terminal, do module load python to fix this.

GPUs
====

Bender
======

Compile as:

::

    make CUDA_VERSION=cc60 COMPILE_CUDA_PATH=/usr/local/cuda-9.2 USE_CUDA=TRUE COMP=PGI -j 4

To run the CUDA code path without device launching, do:

::

    make -j4 COMP=PGI USE_CUDA=TRUE USE_MPI=FALSE DEBUG=TRUE NO_DEVICE_LAUNCH=TRUE CUDA_VERSION=cc60 COMPILE_CUDA_PATH=/usr/local/cuda-9.2
