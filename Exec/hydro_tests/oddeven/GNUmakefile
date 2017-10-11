PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_SHOCK_VAR = TRUE

CASTRO_HOME = ../../..

# This sets the EOS directory in $(CASTRO_HOME)/EOS
EOS_DIR     := gamma_law

# This sets the network directory in $(CASTRO_HOME)/Networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(NETWORK_DIR)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
