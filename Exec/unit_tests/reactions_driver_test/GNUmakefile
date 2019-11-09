PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 3

COMP	         = gnu

USE_MPI          = FALSE
USE_OMP          = FALSE
USE_ACC          = FALSE

USE_DIFFUSION    = FALSE
USE_GRAV         = FALSE
USE_RAD          = FALSE
USE_PARTICLES    = FALSE
USE_ROTATION     = FALSE

USE_REACT        = TRUE

USE_MAESTRO_INIT = FALSE

CASTRO_HOME = ../../..

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the Network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := ignition_simple

# This sets the integrator directory in $(MICROPHYSICS_HOME)/integration
INTEGRATOR_DIR := VODE90

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
