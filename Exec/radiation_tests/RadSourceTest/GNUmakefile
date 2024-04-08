COMP	   = gnu

DIM        = 1

USE_MPI    = TRUE
DEBUG      = FALSE

PROFILE    = FALSE

USE_RAD = TRUE


CASTRO_HOME ?= ../../..

# has the f90 routine ca_initgroups and group structure in probdata_module
#HAS_INITGROUPS = TRUE

EOS_DIR     := rad_power_law

NETWORK_DIR := general_null
NETWORK_INPUTS := gammalaw.net

Opacity_dir := rad_power_law

# use radiation interpbndry
RAD_INTERP = TRUE

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
