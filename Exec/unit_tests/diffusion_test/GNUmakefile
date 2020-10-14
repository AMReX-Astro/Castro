PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_DIFFUSION = TRUE

USE_POST_SIM = TRUE


GPU_COMPATIBLE_PROBLEM = TRUE

CASTRO_HOME := ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := gamma_law_general

CONDUCTIVITY_DIR = constant

# This sets the network directory in $(MICROPHYSICS_HOME)/Networks
NETWORK_DIR := general_null
NETWORK_INPUTS = gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
