PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_THORNADO = TRUE

DEFINES += -DMOMENT_CLOSURE_MINERBO

# define the location of the CASTRO top directory
CASTRO_HOME := ../../..

MICROPHYSICS_HOME := ../../../../Microphysics

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(NETWORK_DIR)/ignition.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

# This sets the network directory in Castro/Networks

Bpack   := ./Make.package 

include $(CASTRO_HOME)/Exec/Make.Castro
