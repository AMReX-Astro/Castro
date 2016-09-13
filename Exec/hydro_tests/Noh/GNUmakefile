# Define the location of the CASTRO top directory,
# if not already defined by an environment variable.

CASTRO_HOME := ../../..

# Location of this directory. Useful if
# you're trying to compile this from another location.

TEST_DIR = $(CASTRO_HOME)/Exec/hydro_tests/Noh

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3
DIMENSION_AGNOSTIC = TRUE

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = TRUE
USE_OMP    = FALSE

EOS_dir   := gamma_law

# This sets the network directory in Castro/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(Network_dir)/gammalaw.net

Bpack   += $(TEST_DIR)/Make.package
Blocs   += $(TEST_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
