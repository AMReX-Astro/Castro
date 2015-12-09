CASTRO_DIR := ../..

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

ifdef MICROPHYSICS_DIR
  # This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
  EOS_dir     := gamma_law_general
else
  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)
endif

# This sets the network directory in $(CASTRO_DIR)/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_DIR)/Networks/$(Network_dir)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
