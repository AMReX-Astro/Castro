PRECISION  = DOUBLE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = TRUE

USE_GRAV      = TRUE
USE_ROTATION  = TRUE
USE_POINTMASS = TRUE

DIMENSION_AGNOSTIC = TRUE

USE_HYBRID_MOMENTUM = TRUE

# define the location of the CASTRO top directory
CASTRO_DIR  := ../..

ifdef MICROPHYSICS_DIR

# This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
EOS_dir     := polytrope

# This sets the EOS directory in $(MICROPHYSICS_DIR)/networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(MICROPHYSICS_DIR)/networks/$(Network_dir)/ignition.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)

endif

Bpack	:= ./Make.package
Blocs	:= .

include ../Make.Castro
