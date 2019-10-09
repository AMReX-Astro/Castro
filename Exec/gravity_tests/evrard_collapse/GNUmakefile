# Define the location of the CASTRO top directory,
# if not already defined by an environment variable.

CASTRO_HOME := ../../..

# Location of this directory. Useful if
# you're trying to compile this from another location.

TEST_DIR = $(CASTRO_HOME)/Exec/gravity_tests/evrard_collapse

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_GRAV   = TRUE

USE_PROB_PARAMS = TRUE

ifdef MICROPHYSICS_HOME
  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_DIR     := gamma_law_general
else
  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)
endif

# This sets the network directory in Castro/Networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(NETWORK_DIR)/gammalaw.net

Bpack   += $(TEST_DIR)/Make.package
Blocs   += $(TEST_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
