PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 2

COMP	         = gnu

USE_MPI          = TRUE
USE_GRAV         = TRUE
USE_REACT        = TRUE

USE_MODEL_PARSER = TRUE

USE_MAESTRO_INIT = FALSE

CASTRO_HOME ?= ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := triple_alpha_plus_cago

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
