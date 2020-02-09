CASTRO_HOME := ../../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = FALSE

GPU_COMPATIBLE_PROBLEM = TRUE

USE_PROB_PARAMS = TRUE

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := gamma_law_general

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS = simple.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
