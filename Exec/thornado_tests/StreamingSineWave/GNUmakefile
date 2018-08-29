PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3
DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_HDF5   = TRUE

USE_THORNADO = TRUE

DEFINES += -DMOMENT_CLOSURE_MINERBO

# define the location of the CASTRO top directory
CASTRO_HOME := ../../..
THORNADO_HOME := ../../../../thornado

MICROPHYSICS_HOME := ../../../../Microphysics

#HDF5_HOME   = /home/almgren/bin/hdf5-1.8.3/hdf5
#HDF5_HOME   = /home/almgren/bin/hdf5
HDF5_HOME   = /home/almgren/bin/miniconda2

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null

GENERAL_NET_INPUTS = ./ignition.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

# This sets the network directory in Castro/Networks

Bpack   := ./Make.package 

include $(CASTRO_HOME)/Exec/Make.Castro
