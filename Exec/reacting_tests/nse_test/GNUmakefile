PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP       = gnu

USE_MPI    = TRUE

USE_REACT  = TRUE

USE_SIMPLIFIED_SDC := TRUE

CASTRO_HOME ?= ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
ifeq ($(USE_NSE_NET), TRUE)
	NETWORK_DIR := ase
	SCREEN_METHOD := chabrier1998
else
	USE_NSE_TABLE := TRUE
	NETWORK_DIR := aprox19
endif

INTEGRATOR_DIR := VODE

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
