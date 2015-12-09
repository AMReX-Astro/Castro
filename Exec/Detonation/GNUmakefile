PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 1

COMP	   = g++
FCOMP	   = gfortran

USE_MPI    = FALSE

USE_REACT  = TRUE

CASTRO_DIR = ../..

ifdef MICROPHYSICS_DIR

# This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
EOS_dir     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_DIR)/networks
Network_dir := aprox13

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
