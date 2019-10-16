PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 1

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE

USE_PROB_PARAMS = TRUE

CASTRO_HOME = ../../..

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(NETWORK_DIR)/ignition.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack	:= ./Make.package
Blocs	:= .

include $(CASTRO_HOME)/Exec/Make.Castro
