PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = TRUE

USE_ROTATION  = TRUE
USE_DIFFUSION = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../../..

GPU_COMPATIBLE_PROBLEM = TRUE

USE_PROB_PARAMS = TRUE

ifdef MICROPHYSICS_HOME

  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_DIR     := helmholtz

  # This sets the network directory in $(MICROPHYSICS_HOME)/networks
  NETWORK_DIR := aprox13

  INTEGRATOR_DIR := VODE
else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

CONDUCTIVITY_DIR := stellar

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
