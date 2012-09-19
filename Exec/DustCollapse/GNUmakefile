CASTRO_DIR ?= /path/to/Castro

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_GRAV   = TRUE

# This sets the EOS directory in $(CASTRO_DIR)/EOS
EOS_dir     := gamma_law_general

# This sets the EOS directory in $(CASTRO_DIR)/Networks
Network_dir := null

Bpack	:= ./Make.package
Blocs	:= .

include $(CASTRO_DIR)/Exec/Make.Castro
