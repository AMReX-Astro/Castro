PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT = TRUE

USE_ROTATION = TRUE
USE_DIFFUSION = TRUE

USE_PROB_PARAMS = TRUE

CASTRO_HOME = ../../..

ifdef MICROPHYSICS_HOME

  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_DIR     := gamma_law_general

  # This sets the network directory in $(MICROPHYSICS_HOME)/networks
  NETWORK_DIR := powerlaw

else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

CONDUCTIVITY_DIR := constant

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
