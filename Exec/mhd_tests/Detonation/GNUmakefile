PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE

USE_REACT  = TRUE

USE_GRAV   = FALSE
USE_PARTICLES = FALSE

USE_MHD = TRUE

USE_PROB_PARAMS = TRUE

CASTRO_HOME = ../../..

GPU_COMPATIBLE_PROBLEM = TRUE

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the network directory in Castro/Networks
NETWORK_DIR := aprox13

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
