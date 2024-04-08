PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP       = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

# define CASTRO top directory
CASTRO_HOME ?= ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := gamma_law

# This sets the network directory in $(MICROPHYSICS_HOME)/Networks
NETWORK_DIR := general_null
NETWORK_INPUTS = gammalaw.net

PROBLEM_DIR ?= ./

Bpack := $(PROBLEM_DIR)/Make.package
Blocs := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
