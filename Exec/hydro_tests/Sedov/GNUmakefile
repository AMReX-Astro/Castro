PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

GPU_COMPATIBLE_PROBLEM = TRUE
USE_PROB_PARAMS = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := gamma_law

# This sets the network directory in $(MICROPHYSICS_HOME)/Networks
NETWORK_DIR := general_null
NETWORK_INPUTS = gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
