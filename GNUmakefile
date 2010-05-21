PBOXLIB_HOME= ../../..
FBOXLIB_HOME= ../../../../fParallel

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 1

COMP	   = g++
FCOMP	   = gfortran

USE_MPI    = FALSE

USE_REACT  = TRUE


# These set the EOS and reaction network directories
EOS_dir     := $(FBOXLIB_HOME)/extern/EOS/helmeos
Network_dir := $(FBOXLIB_HOME)/extern/networks/ignition_simple


Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
