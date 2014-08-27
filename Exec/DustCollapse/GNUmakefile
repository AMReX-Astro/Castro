PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 1
#DIM        = 2
DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_GRAV   = TRUE

# define the location of the CASTRO top directory
CASTRO_DIR  := ../..

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := gamma_law_general

# This sets the EOS directory in $(CASTRO_DIR)/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_DIR)/Networks/$(Network_dir)/ignition.net

Bpack	:= ./Make.package
Blocs	:= .

include ../Make.Castro
