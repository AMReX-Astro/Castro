PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT = TRUE
USE_ROTATION = TRUE
USE_DIFFUSION = TRUE

CASTRO_DIR = ../..

ifdef MICROPHYSICS_DIR

  # This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
  EOS_dir     := gamma_law_general

  # This sets the network directory in $(MICROPHYSICS_DIR)/networks
  Network_dir := kpp

else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)

endif

Conductivity_dir := constant

# This sets the network directory in $(MICROPHYSICS_DIR)/Networks
Network_dir := powerlaw

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
