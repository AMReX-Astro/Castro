PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

# define the location of the CASTRO top directory
CASTRO_DIR  := ../..

# This sets the EOS directory in Castro/EOS
EOS_dir     := gamma_law

# This sets the network directory in Castro/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_DIR)/Microphysics/networks/$(Network_dir)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
