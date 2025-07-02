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

USE_MODEL_PARSER = TRUE
NUM_MODELS := 2

CASTRO_HOME ?= ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := gamma_law

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := powerlaw

CONDUCTIVITY_DIR := constant

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
