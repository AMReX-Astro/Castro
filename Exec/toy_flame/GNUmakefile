PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 1

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_REACT = TRUE

USE_DIFFUSION = TRUE

CASTRO_DIR = ../..

ifdef MICROPHYSICS_DIR

  # This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
  EOS_dir     := gamma_law_general

  # This sets the network directory in $(MICROPHYSICS_DIR)/networks
  Network_dir := powerlaw

else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)

endif

Conductivity_dir := constant

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
