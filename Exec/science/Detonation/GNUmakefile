PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 1

COMP	   = gnu

USE_MPI    = FALSE

USE_REACT  = TRUE

CASTRO_HOME = ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
ifeq ($(USE_NSE_NET), TRUE)
	NETWORK_DIR := ase
	SCREEN_METHOD := chabrier1998
else
	NETWORK_DIR := aprox19
endif

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
