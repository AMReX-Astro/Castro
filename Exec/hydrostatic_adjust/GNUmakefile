PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 1

COMP	   = Intel
FCOMP	   = Intel

COMP	   = g++
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_GRAV   = TRUE

CASTRO_DIR = ../..

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := helmeos

# This sets the EOS directory in $(CASTRO_DIR)/Networks
Network_dir := general_null
#GENERAL_NET_INPUTS := $(CASTRO_DIR)/Networks/general_null/triple_alpha_plus_o.net
GENERAL_NET_INPUTS := $(CASTRO_DIR)/Networks/general_null/ignition.net

Bpack	:= ./Make.package
Blocs	:= .

include $(CASTRO_DIR)/Exec/Make.Castro
