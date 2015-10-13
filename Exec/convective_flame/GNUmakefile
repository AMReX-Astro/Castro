PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT = TRUE
USE_ROTATION = TRUE
USE_DIFFUSION = TRUE

CASTRO_DIR = ../..


# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := gamma_law_general

Conductivity_dir := constant

# This sets the network directory in $(CASTRO_DIR)/Networks
Network_dir := powerlaw

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
