PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

CASTRO_HOME = ../..

# This sets the EOS directory in $(CASTRO_HOME)/EOS
EOS_dir     := gamma_law

# This sets the network directory in $(CASTRO_HOME)/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(Network_dir)/gammalaw.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
