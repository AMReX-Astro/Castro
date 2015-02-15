PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

CASTRO_DIR = ../..

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := gamma_law_general

# This sets the network directory in $(CASTRO_DIR)/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_DIR)/Networks/$(Network_dir)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
