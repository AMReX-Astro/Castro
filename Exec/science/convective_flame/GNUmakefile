PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT = TRUE

USE_ROTATION = TRUE
USE_DIFFUSION = TRUE

CASTRO_HOME = ../../..

ifdef MICROPHYSICS_HOME

  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_dir     := gamma_law_general

  # This sets the network directory in $(MICROPHYSICS_HOME)/networks
  Network_dir := powerlaw

else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Conductivity_dir := constant

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
