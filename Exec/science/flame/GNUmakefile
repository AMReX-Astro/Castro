PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 1

COMP	   = gnu

USE_MPI    = TRUE

USE_REACT = TRUE

USE_DIFFUSION = TRUE

CASTRO_HOME = ../../..

ifdef MICROPHYSICS_HOME

  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_DIR     := helmholtz

  # This sets the network directory in $(MICROPHYSICS_HOME)/networks
  NETWORK_DIR := triple_alpha_plus_cago

else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

#Conductivity_dir := constant_opacity
CONDUCTIVITY_DIR := stellar

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
