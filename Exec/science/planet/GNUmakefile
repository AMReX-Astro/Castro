COMP	         = gnu

DIM              = 3

PRECISION        = DOUBLE

USE_MPI          = TRUE
USE_OMP          = FALSE

DEBUG            = FALSE

PROFILE          = FALSE


USE_RAD = TRUE

USE_GRAV         = TRUE
USE_REACT        = FALSE
USE_MODELPARSER = TRUE
USE_PARTICLE     = FALSE
USE_SHOCK_VAR    = TRUE
USE_ROTATION = FALSE

USE_MODEL_PARSER = TRUE

CASTRO_HOME ?= ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := gamma_law

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null
NETWORK_INPUTS := ./taeho.net

# power-law opacity
OPACITY_DIR := rad_power_law

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
