
PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

#DIM        = 1
#DIM        = 2
DIM        = 1

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

# define the location of the CASTRO top directory
CASTRO_DIR  := ../..

# This sets the EOS directory in Castro/EOS
#EOS_dir     := helmeos
EOS_dir     := gamma_law_general

# This sets the network directory in Castro/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_DIR)/Networks/$(Network_dir)/ignition.net

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
