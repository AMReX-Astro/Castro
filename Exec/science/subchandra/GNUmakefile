CASTRO_HOME := ../../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = TRUE

USE_CXX_MODEL_PARSER = TRUE

USE_SIMPLIFIED_SDC=TRUE

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := subch_approx

INTEGRATOR_DIR := VODE

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
