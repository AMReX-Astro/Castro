PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 1

COMP	   = gnu

USE_MPI    = FALSE

USE_REACT  = TRUE

USE_PROB_PARAMS = TRUE

CASTRO_HOME = ../../..

GPU_COMPATIBLE_PROBLEM = TRUE

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := aprox19


else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
