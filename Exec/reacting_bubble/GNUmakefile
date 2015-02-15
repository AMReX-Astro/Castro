PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 2

COMP	         = g++
FCOMP	         = gfortran

USE_MPI          = TRUE
USE_GRAV         = TRUE
USE_REACT        = TRUE

USE_MODELPARSER  = TRUE
USE_MAESTRO_INIT = FALSE
USE_OLDPLOTPER   = FALSE


CASTRO_DIR = ../..

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := helmeos

# This sets the network directory in $(CASTRO_DIR)/Networks
Network_dir := ignition_simple

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
