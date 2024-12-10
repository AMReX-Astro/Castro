COMP	   = gnu

DIM        = 1

USE_MPI    = TRUE
DEBUG      = FALSE

PROFILE    = FALSE

USE_GRAV   = TRUE

USE_RAD = TRUE

USE_MODEL_PARSER = TRUE

CASTRO_HOME ?= ../../..

EOS_DIR     = breakout
NETWORK_DIR := general_null
NETWORK_INPUTS = breakout.net
OPACITY_DIR = breakout

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
