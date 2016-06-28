
PRECISION  = DOUBLE
PROFILE    = TRUE
COMM_PROFILE    = TRUE

DEBUG      = TRUE

DIM        = 3

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../..

ifdef MICROPHYSICS_HOME
  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_dir     := gamma_law_general
else
  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)
endif

# This sets the EOS directory in Castro/EOS
EOS_dir     := gamma_law_general

# This sets the network directory in Castro/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(Network_dir)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
