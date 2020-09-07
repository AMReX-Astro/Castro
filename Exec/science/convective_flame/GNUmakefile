PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT = TRUE

USE_ROTATION = TRUE
USE_DIFFUSION = TRUE

USE_PROB_PARAMS = TRUE

USE_INTERPOLATE = TRUE

CASTRO_HOME = ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := gamma_law_general

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := powerlaw

CONDUCTIVITY_DIR := constant

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
