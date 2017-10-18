
PRECISION  = DOUBLE
PROFILE    = TRUE
COMM_PROFILE    = TRUE

DEBUG      = TRUE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../../..

ifdef MICROPHYSICS_HOME
  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_DIR     := gamma_law_general
else
  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)
endif

# This sets the EOS directory in Castro/EOS
EOS_DIR     := gamma_law_general

# This sets the network directory in Castro/Networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(NETWORK_DIR)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
