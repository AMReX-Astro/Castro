PRECISION  = DOUBLE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_GRAV      = TRUE
USE_ROTATION  = TRUE

USE_HYBRID_MOMENTUM = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  ?= ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := polytrope

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null
NETWORK_INPUTS = ignition.net

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
