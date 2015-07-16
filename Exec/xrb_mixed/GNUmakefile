PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 2
COMP	         = g++
FCOMP	         = gfortran

USE_MPI          = TRUE
USE_OMP          = FALSE

USE_GRAV         = TRUE
USE_REACT        = TRUE
USE_MODELPARSER  = TRUE

USE_OLDPLOTPER   = TRUE


# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := helmeos

# This sets the network directory in $(NETWORK_HOME)
NETWORK_HOME := $(MAESTRO_HOME)/Microphysics/networks/
Network_dir := rprox

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
