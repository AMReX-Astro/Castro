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
CASTRO_HOME  := ../..

# This sets the EOS directory in $(CASTRO_HOME)/EOS
EOS_dir     := gamma_law

# This sets the EOS directory in $(CASTRO_HOME)/Networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(CASTRO_HOME)/Microphysics/networks/$(Network_dir)/ignition.net

Bpack	:= ./Make.package
Blocs	:= .

include ../Make.Castro
