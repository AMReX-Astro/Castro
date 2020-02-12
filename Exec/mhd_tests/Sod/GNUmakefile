PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_MHD    = TRUE
USE_PROB_PARAMS = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../../..

# This sets the EOS directory in Castro/EOS
EOS_DIR     := gamma_law

# This sets the network directory in Castro/Networks
NETWORK_DIR := general_null
NETWORK_INPUTS = gammalaw.net


Blocs   := .
Bpack   := ./Make.package

include $(CASTRO_HOME)/Exec/Make.Castro
