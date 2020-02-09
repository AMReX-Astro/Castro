
PRECISION  = DOUBLE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_PROB_PARAMS = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := gamma_law_general

# This sets the network directory in Castro/Networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(NETWORK_DIR)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
