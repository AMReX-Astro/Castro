CASTRO_HOME := ../../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = Intel
FCOMP	   = Intel

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = TRUE

USE_MODELPARSER  = TRUE

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
EOS_dir     := gamma_law_general

# This sets the network directory in $(MICROPHYSICS_DIR)/networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(Network_dir)/H_He.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
