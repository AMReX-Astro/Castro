PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_GRAV   = TRUE
USE_MAESTRO_INIT = FALSE
USE_OLDPLOTPER = TRUE
USE_REACT = TRUE

CASTRO_DIR = ../..

ifdef MICROPHYSICS_DIR

  # This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
  EOS_dir     := gamma_law_general

  # This sets the network directory in $(MICROPHYSICS_DIR)/networks
  Network_dir := kpp

else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
