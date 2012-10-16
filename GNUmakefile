CASTRO_DIR ?= /path/to/Castro

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 1

COMP	   = g++
FCOMP	   = gfortran

USE_MPI    = FALSE

USE_REACT  = TRUE

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := helmeos

# This sets the EOS directory in $(CASTRO_DIR)/Networks
NETWORK_HOME := $(ASTRODEV_DIR)/networks
Network_dir := approx8

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
