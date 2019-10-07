CASTRO_HOME := ../../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_GRAV   = TRUE

USE_PROB_PARAMS = TRUE

USE_MAESTRO_INIT = FALSE

ifdef MICROPHYSICS_HOME
  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_DIR     := gamma_law_general
else
  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)
endif

# This sets the network directory in $(CASTRO_HOME)/Networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(NETWORK_DIR)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
