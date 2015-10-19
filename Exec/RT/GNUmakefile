PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE

USE_GRAV   = TRUE

# define the location of the CASTRO top directory
CASTRO_DIR  := ../..

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := gamma_law

# This sets the network directory in $(CASTRO_DIR)/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_DIR)/Networks/$(Network_dir)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
