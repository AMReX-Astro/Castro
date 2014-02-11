PRECISION  = DOUBLE

DEBUG      = FALSE
DEBUG      = TRUE

PROFILE    = FALSE

USE_MPI    = TRUE

USE_HYPRE  = TRUE
USE_RAD    = TRUE

USE_NEUTRINO = TRUE

# uncomment these three to build without radiation
#USE_HYPRE  = FALSE
#USE_RAD    = FALSE
#USE_NEUTRINO = FALSE

# use radiation interpbndry
RAD_INTERP = TRUE

USE_GRAV       = TRUE
USE_AUX_UPDATE = FALSE

COMP	   = g++
FCOMP      = gfortran
#COMP	   = Intel
#FCOMP	   = Intel
DIM        = 1
#DIM        = 2
#DIM        = 3

EOS_dir     := collapse
Network_dir := collapse
Opacity_dir := burrows

ifeq ($(EOS_dir), cott)
   USE_HDF5 = TRUE
endif

ifeq ($(Opacity_dir), cott)
   USE_HDF5 = TRUE
endif

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro

