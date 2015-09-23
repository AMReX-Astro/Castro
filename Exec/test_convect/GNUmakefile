CASTRO_DIR ?= ../..

PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 2
COMP	         = g++
FCOMP	         = gfortran

USE_MPI          = TRUE
USE_OMP          = FALSE

USE_GRAV         = TRUE
USE_REACT        = FALSE
USE_MODELPARSER  = TRUE

USE_OLDPLOTPER   = TRUE

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := helmeos

# This sets the network directory in $(NETWORK_HOME)
Network_dir := general_null
GENERAL_NET_INPUTS := $(CASTRO_DIR)/Networks/$(Network_dir)/ignition.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
