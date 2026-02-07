CASTRO_HOME ?= ../../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = TRUE
USE_SHOCK_VAR = TRUE

USE_MODEL_PARSER = TRUE

USE_SIMPLIFIED_SDC = TRUE

# for the networks we tend to use, we can still use Jac caching on GPUs
USE_JACOBIAN_CACHING = TRUE

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := he-burn/he-burn-19am

INTEGRATOR_DIR := VODE

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
