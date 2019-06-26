PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 2

COMP	         = gnu

USE_MPI          = TRUE
USE_GRAV         = TRUE
USE_REACT        = TRUE

USE_MAESTRO_INIT = FALSE

CASTRO_HOME = ../../..

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := triple_alpha_plus_cago

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
