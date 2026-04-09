PRECISION ?= DOUBLE
PROFILE ?= FALSE

DEBUG ?= FALSE
BACKTRACE ?= FALSE
TEST ?= FALSE

USE_MPI ?= TRUE
USE_OMP ?= FALSE

DIM ?= 3

USE_GRAV ?= TRUE
USE_REACT ?= TRUE
USE_ROTATION ?= TRUE

USE_MODEL_PARSER = TRUE
NUM_MODELS = 2

# Disable burning in shocksa
USE_SHOCK_VAR = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  ?= ../../..

# Turn on hybrid momentum? but note that it only is used in 3D at present.

ifeq ($(DIM), 3)
  USE_HYBRID_MOMENTUM ?= FALSE
endif

# Set the default integrator for the burner.

INTEGRATOR_DIR ?= VODE

# Default compiler.

COMP ?= gnu

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS

EOS_DIR ?= helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks

NETWORK_DIR ?= aprox13

# If using general_null, make sure we have a reasonable set of isotopes.

ifeq ($(NETWORK_DIR), general_null)
  NETWORK_INPUTS = aprox13.net
endif

# Location where the source files are coming from.

PROBLEM_DIR ?= ./

Blocs = $(PROBLEM_DIR)
Bpack = $(PROBLEM_DIR)/Make.package

include $(CASTRO_HOME)/Exec/Make.Castro
