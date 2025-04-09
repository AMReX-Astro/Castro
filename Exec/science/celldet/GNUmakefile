PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE

USE_REACT  = TRUE

USE_GRAV   = FALSE
USE_PARTICLES = FALSE

USE_MHD = FALSE

#USE_PROB_PARAMS = TRUE

USE_MODEL_PARSER = TRUE

CASTRO_HOME ?= ../../..

#GPU_COMPATIBLE_PROBLEM = TRUE

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := aprox13

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
